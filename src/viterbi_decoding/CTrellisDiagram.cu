
#include "CTrellisDiagram.cuh"

#include <limits.h>
#include <stddef.h>
#include <algorithm>
#include <iostream>
#include <cmath>

#ifdef _USE_MPI_
#include <mpi/mpi.h> // For MPI API
#endif

#ifdef _USE_OPENMP_
#include <omp.h> // For OPENMP API
#endif

#ifdef _USE_CUDA_
#include <cuda_runtime.h> // For CUDA API
#include <vector>
#endif

#include <list>

#include <stdio.h>

#include "../lookup_tables/CNextStateLookupTable.cuh"
#include "../lookup_tables/CNextEncodedSymbolLookupTable.h"

const unsigned int CTrellisDiagram::s_cnInfinity = UINT_MAX / 2;

#ifdef _USE_CUDA_
__constant__ unsigned int g_cpnSetBitsInNibble[] =
{ 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
__constant__ unsigned int g_cnInfinity = UINT_MAX / 2;
#endif

const unsigned int CTrellisDiagram::s_pcnSetBitsInNibble[] =
{ 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };


CTrellisDiagram::CTrellisDiagram(unsigned int const p_cnNumberOfSymbols, unsigned int const p_cnK, CNextStateLookupTable const& p_rccNextStateLookupTable,
		CNextEncodedSymbolLookupTable const& p_rccNextEncodedSymbolLookupTable, unsigned int const p_cnSizeOfDecodedDataType) :
		m_nRowsCount(0), m_nColumnsCount(0), m_nNumberOfSymbols(p_cnNumberOfSymbols), m_cnMemory(p_cnK), m_rccNextStateLookupTable(
				p_rccNextStateLookupTable), m_rccNextEncodedSymbolLookupTable(p_rccNextEncodedSymbolLookupTable), m_cnNibblesInDataType(
				(p_cnSizeOfDecodedDataType * 8) / 4), m_nJump(0)
{

#ifdef _USE_MPI_
	// Setup API configurations and needed variables
	InitializeMpi();
#endif

	m_nMaxState = CStateLookupTable::GetMaxState(m_cnMemory);

	// Initialize this diagram only in construction
	// In new trellis calculation with the same construction properties,
	// only ClearTrellisDiagram() is needed
	InitializeTrellisDiagram();

#ifdef _USE_CUDA_
	SetGpuInfo(); // Get GPU properties

	// copy initial data
	AllocateAndSetDeviceMemory();
#endif


}

CTrellisDiagram::~CTrellisDiagram(void)
{
#ifdef _USE_MPI_
	// Finalize MPI if needed at this object destruction
	int l_nMpiFinilized;
	MPI_Finalized(&l_nMpiFinilized);
	if(!l_nMpiFinilized)
	{
		MPI_Finalize();
	}
#endif

#ifdef _USE_CUDA_
	FreeMemoryAndResetDevice();
#endif
}

std::list<std::vector<unsigned int> > const& CTrellisDiagram::GetDistancesMatrix(void) const
{
	return m_cDistancesMatrix;
}

void CTrellisDiagram::InitializeTrellisDiagram(void)
{

#ifdef _USE_MPI_
	// In case of MPI, change the number of symbols according to the current number of processes
	// Each process gets an equal amount of symbols to decode

	unsigned int const l_cnDataTypeSizeInBits = m_cnNibblesInDataType * 4;
	unsigned int const l_cnNumberOfElementsInDecodedSymbols = (m_nNumberOfSymbols - m_cnMemory) / l_cnDataTypeSizeInBits;
	unsigned int const l_cnElementsPerProcess = l_cnNumberOfElementsInDecodedSymbols / m_nNumberOfProcesses;
	m_nNumberOfSymbols = l_cnElementsPerProcess * l_cnDataTypeSizeInBits;
	if(m_nMyRank == m_nNumberOfProcesses - 1)
	{
		unsigned int const l_cnRemainElements = l_cnNumberOfElementsInDecodedSymbols % m_nNumberOfProcesses;
		if(l_cnRemainElements != 0)
		{
			m_nNumberOfSymbols += l_cnRemainElements;
		}
	}
	m_nNumberOfSymbols += m_cnMemory; // Complete spares after cut them in the calculation above
#endif

	m_nColumnsCount = m_nNumberOfSymbols + 1;
	m_nRowsCount = m_nMaxState + 1;

	for(unsigned int l_nI = 0; l_nI < m_nColumnsCount; l_nI++)
	{
		std::vector<unsigned int> l_cColumn(m_nRowsCount);
		m_cDistancesMatrix.push_back(l_cColumn);
	}

}

void CTrellisDiagram::ClearTrellisDiagram(void)
{
	m_nJump = std::pow(2, m_cnMemory);

	ClearDistancesMatrix(); // Fill distances matrix with infinity
}

void CTrellisDiagram::ClearDistancesMatrix(void)
{
	for(std::list<std::vector<unsigned int> >::iterator l_cIt = m_cDistancesMatrix.begin(); l_cIt != m_cDistancesMatrix.end(); ++l_cIt)
	{
		std::fill(&(*l_cIt)[0], &(*l_cIt)[0] + m_nRowsCount, s_cnInfinity);
	}
}

void CTrellisDiagram::PerformStages(std::vector<std::list<unsigned int>* > const& p_rccSymbols)
{


	// Prepare diagram for calculation
	ClearTrellisDiagram();

	std::list<std::vector<unsigned int> >::iterator l_cTrellisIt = m_cDistancesMatrix.begin();

#ifdef _USE_MPI_

	// Calculate the position of the current process buffer according to the parallel division
	std::list<unsigned int>* const* l_nStartPtr = &p_rccSymbols[0] + (m_nNumberOfSymbols - m_cnMemory) * m_nMyRank;

	// Perform first iteration for first symbol
	PerformFirstStages(**l_nStartPtr, l_cTrellisIt);
	l_cTrellisIt++;

	// Perform the rest of the iteration for the rest of the symbols
	std::list<unsigned int>* const* l_pcpcSymbolsPtrLast = l_nStartPtr + m_nNumberOfSymbols;
	for (std::list<unsigned int>** l_pnSymbolsPtr = const_cast<std::list<unsigned int>**>(l_nStartPtr) + 1; l_pnSymbolsPtr < l_pcpcSymbolsPtrLast; l_pnSymbolsPtr++)
	{
		PerformSequentialStage(**l_pnSymbolsPtr, l_cTrellisIt);
		l_cTrellisIt++;

		if(m_nJump != 1)
		{
			m_nJump /= 2;
		}
	}


#else

	// Perform first iteration for first symbol
	PerformFirstStages(*(p_rccSymbols[0]), l_cTrellisIt);
	l_cTrellisIt++;
	if(m_nJump != 1)
	{
		m_nJump /= 2;
	}

	// Perform the rest of the iteration for the rest of the symbols
	std::list<unsigned int>* const* l_pcpcSymbolsPtrLast = &p_rccSymbols[0] + m_nNumberOfSymbols;
	for (std::list<unsigned int>** l_ppcSymbolsPtr = const_cast<std::list<unsigned int>**>(&p_rccSymbols[0] + 1); l_ppcSymbolsPtr < l_pcpcSymbolsPtrLast; l_ppcSymbolsPtr++)
	{
		PerformSequentialStage(**l_ppcSymbolsPtr, l_cTrellisIt);
		l_cTrellisIt++;

		if(m_nJump != 1)
		{
			m_nJump /= 2;
		}
	}

#endif

}

void CTrellisDiagram::PerformFirstStages(std::list<unsigned int> const& p_rccSymbol, std::list<std::vector<unsigned int> >::iterator& p_rcTrellisItPtr)
{
	// This is the trellis calculation for the first iteration


	// Start from the possible states and symbols, from the first state in the table
	unsigned int const* l_pcnNextPossibleStates = m_rccNextStateLookupTable.GetSequentialStates(0);
	std::list<unsigned int> const* l_pccNextPossibleSymbols = m_rccNextEncodedSymbolLookupTable.GetNextEncodedSymbols(0);

	std::vector<unsigned int>& l_cCurrentColumn = *(p_rcTrellisItPtr++);
	std::vector<unsigned int>& l_cNextColumn = *(p_rcTrellisItPtr--);

	// Perform calculation for the possible destination states
	for (unsigned int l_nInput = 0; l_nInput < 2; l_nInput++)
	{
		// 1. Calculate the hamming distance between the received symbol and the predicted symbol from the table
		// 2. Get the initial distance "already calculated" from the distances matrix
		// 3. Take the minimum between those two values
		unsigned int l_nDistance = std::min(l_cCurrentColumn[0],
				GetHammingDistance(p_rccSymbol, l_pccNextPossibleSymbols[l_nInput], m_cnNibblesInDataType));

		l_cNextColumn[l_pcnNextPossibleStates[l_nInput]] = l_nDistance;
	}

}

void CTrellisDiagram::PerformSequentialStage(std::list<unsigned int> const& p_rccSymbol, std::list<std::vector<unsigned int> >::iterator& p_rcTrellisItPtr)
{
	// This is the trellis calculation for the rest of the iterations

	std::vector<unsigned int>& l_cCurrentColumn = *(p_rcTrellisItPtr++);
	std::vector<unsigned int>& l_cNextColumn = *(p_rcTrellisItPtr--);

#ifdef _USE_CUDA_


	// Take the symbol as vector instead of linked list in order to reduce the number of calls to cudaMemcpy
	std::vector<unsigned int> l_cSymbolAsVector(m_rccNextEncodedSymbolLookupTable.GetNumberOfClusters());
	unsigned int l_nCurrentCluster = 0;
	for(std::list<unsigned int>::const_iterator l_cIt = p_rccSymbol.begin(); l_cIt != p_rccSymbol.end(); ++l_cIt)
	{
		l_cSymbolAsVector[l_nCurrentCluster++] = *l_cIt;
	}

	// Copy the symbol to the device global memory
	cudaMemcpy(m_pnIterationSymbolDevicePtr, &l_cSymbolAsVector[0], sizeof(unsigned int) * m_rccNextEncodedSymbolLookupTable.GetNumberOfClusters(), cudaMemcpyHostToDevice);

	// Get maximum number of threads per block and calculate the number of blocks
	// that we need
	unsigned int const l_cnNumberOfThreadsPerBlock = m_sGpuProps.maxThreadsPerBlock / 2;
	unsigned int l_nNumberOfBlocks = (m_nRowsCount / l_cnNumberOfThreadsPerBlock == 0) ? 1 : m_nRowsCount / l_cnNumberOfThreadsPerBlock;

	// Get maximum number of blocks according to the GPU properties (only one dimension needed)
	const unsigned int l_cnMaxNumberOfBlocks = m_sGpuProps.maxThreadsDim[0];

	// copy two columns to device memory
	cudaMemcpy(m_pnCurrentAndNextColumnDevicePtr, &l_cCurrentColumn[0], sizeof(unsigned int) * m_nRowsCount, cudaMemcpyHostToDevice);
	cudaMemcpy(m_pnCurrentAndNextColumnDevicePtr + m_nRowsCount, &l_cNextColumn[0], sizeof(unsigned int) * m_nRowsCount, cudaMemcpyHostToDevice);

	// If all the threads can be run at once
	if(l_nNumberOfBlocks <= l_cnMaxNumberOfBlocks)
	{
		// call appropriate kernel ( <<< >>> )
		TrellisIteration_Kernel<<<l_nNumberOfBlocks, l_cnNumberOfThreadsPerBlock>>>(m_pnIterationSymbolDevicePtr, m_pnCurrentAndNextColumnDevicePtr, m_pnNextSymbolTableDevicePtr, m_pnNextStateTableDevicePtr, m_nRowsCount, m_cnNibblesInDataType, m_rccNextEncodedSymbolLookupTable.GetNumberOfClusters(), 0, m_nJump);
	}
	else
	{

		unsigned int l_nCurrentState = 0;
		while(l_nNumberOfBlocks >= l_cnMaxNumberOfBlocks)
		{
			// call kernel ( <<< >>> )
			TrellisIteration_Kernel<<<l_cnMaxNumberOfBlocks, l_cnNumberOfThreadsPerBlock>>>(m_pnIterationSymbolDevicePtr, m_pnCurrentAndNextColumnDevicePtr, m_pnNextSymbolTableDevicePtr, m_pnNextStateTableDevicePtr, m_nRowsCount, m_cnNibblesInDataType, m_rccNextEncodedSymbolLookupTable.GetNumberOfClusters(), l_nCurrentState, m_nJump);

			l_nNumberOfBlocks -= l_cnMaxNumberOfBlocks;
			l_nCurrentState += l_cnMaxNumberOfBlocks * l_cnNumberOfThreadsPerBlock;
		}
		if(l_nNumberOfBlocks > 0)
		{
			// call kernel ( <<< >>> )
			TrellisIteration_Kernel<<<l_nNumberOfBlocks, l_cnNumberOfThreadsPerBlock>>>(m_pnIterationSymbolDevicePtr, m_pnCurrentAndNextColumnDevicePtr, m_pnNextSymbolTableDevicePtr, m_pnNextStateTableDevicePtr, m_nRowsCount, m_cnNibblesInDataType, m_rccNextEncodedSymbolLookupTable.GetNumberOfClusters(), l_nCurrentState, m_nJump);
		}


	}

	// copy results from device to host ( cudaMemcpy ( ... cudaMemcpyDeviceToHost ) )
	cudaMemcpy(&l_cNextColumn[0], m_pnCurrentAndNextColumnDevicePtr + m_nRowsCount, sizeof(unsigned int) * m_nRowsCount, cudaMemcpyDeviceToHost);

#else
#ifdef _USE_OPENMP_
	// Divide the loop work between the machine cores
	// In fact, every core gets an equal amount of states and the load balancing for
	// each state is identical
#pragma omp parallel for
#endif
	// Loop on all active states
	for (unsigned int l_nCurrentState = 0; l_nCurrentState < m_nRowsCount; l_nCurrentState += m_nJump)
	{
		// Get the possible states and symbols, from the current active state
		unsigned int const* l_pcnNextPossibleStates = m_rccNextStateLookupTable.GetSequentialStates(l_nCurrentState);
		std::list<unsigned int> const* l_pccNextPossibleSymbols = m_rccNextEncodedSymbolLookupTable.GetNextEncodedSymbols(l_nCurrentState);

		// Get the distance of the state from which we got to the current one
		unsigned int l_nPreviousDistanceOfInterest = l_cCurrentColumn[l_nCurrentState];

		// Check the possible transition states (for inputs 0 and 1)
		for (unsigned int l_nInput = 0; l_nInput < 2; l_nInput++)
		{
			// Get current possible state for the given input
			unsigned int l_nNextPossibleState = l_pcnNextPossibleStates[l_nInput];

			// 1. Calculate the hamming distance between the current symbol and the predicted symbol,
			//	  plus the distance calculated in the previous iteration
			// 2. Get the distance from the next possible state (in the next iteration. It could already
			//	  calculated from other state in this iteration)
			// 3. Take the minimum between those two values above
			unsigned int l_nDistance = std::min(l_cNextColumn[l_nNextPossibleState],
					GetHammingDistance(p_rccSymbol, l_pccNextPossibleSymbols[l_nInput], m_cnNibblesInDataType)
					+ l_nPreviousDistanceOfInterest);


			// Update the distances matrix with the new distance
			l_cNextColumn[l_nNextPossibleState] = l_nDistance;
		}
	}

#endif

	/****** DEBUG START ******/
//	std::cout << "Trellis iteration performed for symbol : ";
//	for(std::list<unsigned int>::const_iterator it = p_rccSymbol.begin(); it != p_rccSymbol.end(); ++it) {
//		std::cout << "-> " << CPolynomialGenerator::GetBinaryString(*it, 32) << " ";
//	}
//	std::cout << std::endl;
//	std::cout << "and the diagram now is :" << std::endl;
//	std::cout << *this;
//	std::cout << "##########################################################################################################" << std::endl;
	/*************************/

}

unsigned int CTrellisDiagram::GetHammingDistance(unsigned int const p_cnFirstSymbol, unsigned int const p_cnSecondSymbol,
		unsigned int const p_cnNibblesInDataType)
{
	unsigned int l_nXor = p_cnFirstSymbol ^ p_cnSecondSymbol;
	unsigned int l_nSetBitsCount = 0;
	for (unsigned int l_nNibblesCount = 0; l_nNibblesCount < p_cnNibblesInDataType; l_nNibblesCount++)
	{
		l_nSetBitsCount += s_pcnSetBitsInNibble[(l_nXor >> (l_nNibblesCount * 4)) & 0xf];
	}

	return l_nSetBitsCount;
}

unsigned int CTrellisDiagram::GetHammingDistance(std::list<unsigned int> const& p_rccFirstSymbol, std::list<unsigned int> const& p_rccSecondSymbol,
		unsigned int const p_cnNibblesInDataType)
{
	unsigned int l_nSetBitsCount = 0;

	if(p_rccFirstSymbol.size() == p_rccSecondSymbol.size())
	{

		std::list<unsigned int>::const_iterator l_cFIt = p_rccFirstSymbol.begin();
		std::list<unsigned int>::const_iterator l_cSIt = p_rccSecondSymbol.begin();

		for(
				;
				l_cFIt != p_rccFirstSymbol.end() && l_cSIt != p_rccSecondSymbol.end();
				++l_cFIt, ++l_cSIt)
		{
			l_nSetBitsCount += GetHammingDistance(*l_cFIt, *l_cSIt, p_cnNibblesInDataType);
		}
	}

	return l_nSetBitsCount;
}

unsigned int CTrellisDiagram::GetRowsCount(void) const
{
	return m_nRowsCount;
}

unsigned int CTrellisDiagram::GetColumnsCount(void) const
{
	return m_nColumnsCount;
}

std::ostream& operator<<(std::ostream& p_rcOs, CTrellisDiagram const& p_rccTrellis)
{
	p_rcOs << "trellis diagram distances table (" << p_rccTrellis.m_nRowsCount << "x" << p_rccTrellis.m_nColumnsCount << "):" << std::endl;

	for(unsigned int l_nRI = 0; l_nRI < p_rccTrellis.m_nRowsCount; l_nRI++)
	{
		p_rcOs << CPolynomialGenerator::GetBinaryString(l_nRI, p_rccTrellis.m_cnMemory) << "\t";

		for(std::list<std::vector<unsigned int> >::const_iterator l_cIt = p_rccTrellis.m_cDistancesMatrix.begin(); l_cIt != p_rccTrellis.m_cDistancesMatrix.end(); ++l_cIt)
		{
			unsigned int const* l_pnCurrentColumn = &(*l_cIt)[0];

			if(l_pnCurrentColumn[l_nRI] == CTrellisDiagram::s_cnInfinity)
			{
				p_rcOs << "Inf";
			}
			else
			{
				p_rcOs << l_pnCurrentColumn[l_nRI];
			}
			p_rcOs << "\t";
		}

		p_rcOs << std::endl;
	}

    return p_rcOs;
}

#ifdef _USE_MPI_
void CTrellisDiagram::InitializeMpi(void)
{
	// Before initializing MPI, check if it already initialized
	int l_nMpiInitialized;
	MPI_Initialized(&l_nMpiInitialized);
	if(!l_nMpiInitialized)
	{
		MPI_Init(NULL, NULL); // Initialize
	}

	// Get current process rank and number of runnning processes
	MPI_Comm_rank(MPI_COMM_WORLD, &m_nMyRank);
	m_bImMaster = (0 == m_nMyRank) ? true : false; // Set master identification for current process

	MPI_Comm_size(MPI_COMM_WORLD, &m_nNumberOfProcesses);

}
#endif

#ifdef _USE_CUDA_

void CTrellisDiagram::SetGpuInfo(void)
{
	// Get the device properties into the class member
	cudaGetDeviceProperties((cudaDeviceProp*)&m_sGpuProps, 0);
}

void CTrellisDiagram::AllocateAndSetDeviceMemory(void)
{
	m_pnCurrentAndNextColumnDevicePtr = NULL;
	m_pnNextSymbolTableDevicePtr = NULL;
	m_pnNextStateTableDevicePtr = NULL;
	m_pnIterationSymbolDevicePtr = NULL;
	// allocate memory on the device ( cudaMalloc ) - can be done outside the function for all the iterations
	cudaMalloc((void**)&m_pnCurrentAndNextColumnDevicePtr, sizeof(unsigned int) * m_nRowsCount * 2);
	cudaMalloc((void**)&m_pnNextSymbolTableDevicePtr, sizeof(unsigned int) * m_rccNextEncodedSymbolLookupTable.GetSize() * m_rccNextEncodedSymbolLookupTable.GetNumberOfClusters());
	cudaMalloc((void**)&m_pnNextStateTableDevicePtr, sizeof(unsigned int) * m_rccNextStateLookupTable.GetSize());
	cudaMalloc((void**)&m_pnIterationSymbolDevicePtr, sizeof(unsigned int) * m_rccNextEncodedSymbolLookupTable.GetNumberOfClusters());

	// Copy the next state lookup table
	cudaMemcpy(m_pnNextStateTableDevicePtr, &(m_rccNextStateLookupTable.AsVector()[0]), sizeof(unsigned int) * m_rccNextStateLookupTable.GetSize(), cudaMemcpyHostToDevice);

	// Copy the next symbol lookup table
	std::vector<std::list<unsigned int> > const& l_rccNextSymbolTable = m_rccNextEncodedSymbolLookupTable.AsVector();

	// Take the symbol as vector instead of linked list in order to reduce the number of calls to cudaMemcpy
	std::vector<unsigned int> l_cSymbolAsVector(m_rccNextEncodedSymbolLookupTable.GetNumberOfClusters());

	for(unsigned int l_nTI = 0; l_nTI < l_rccNextSymbolTable.size(); l_nTI++)
	{
		unsigned int l_nCurrentCluster = 0;
		std::list<unsigned int> const& l_rccSymbol = l_rccNextSymbolTable[l_nTI];
		for(std::list<unsigned int>::const_iterator l_cIt = l_rccSymbol.begin(); l_cIt != l_rccSymbol.end(); ++l_cIt)
		{
			l_cSymbolAsVector[l_nCurrentCluster++] = *l_cIt;
		}
		cudaMemcpy(m_pnNextSymbolTableDevicePtr + l_nTI * m_rccNextEncodedSymbolLookupTable.GetNumberOfClusters(), &l_cSymbolAsVector[0], sizeof(unsigned int) * m_rccNextEncodedSymbolLookupTable.GetNumberOfClusters(), cudaMemcpyHostToDevice);
	}
}

void CTrellisDiagram::FreeMemoryAndResetDevice(void)
{
	cudaFree(m_pnCurrentAndNextColumnDevicePtr);
	cudaFree(m_pnNextSymbolTableDevicePtr);
	cudaFree(m_pnNextStateTableDevicePtr);
	cudaFree(m_pnIterationSymbolDevicePtr);
	cudaDeviceReset(); // Reset the device at object destruction
}

__device__ unsigned int CudaHammingDistance(unsigned int const* p_pcnFirstSymbol, unsigned int const* p_pcnSecondSymbol, unsigned int const p_cnNibblesInDataType, unsigned int const p_cnNumberOfClusters)
{
	unsigned int l_nSetBitsCount = 0;

	for(unsigned int l_nI = 0; l_nI < p_cnNumberOfClusters; l_nI++)
	{
		unsigned int l_nXor = p_pcnFirstSymbol[l_nI] ^ p_pcnSecondSymbol[l_nI];
		for (unsigned int l_nNibblesCount = 0; l_nNibblesCount < p_cnNibblesInDataType; l_nNibblesCount++)
		{
			l_nSetBitsCount += g_cpnSetBitsInNibble[(l_nXor >> (l_nNibblesCount * 4)) & 0xf];
		}
	}

	return l_nSetBitsCount;
}

__global__ void TrellisIteration_Kernel(unsigned int const* p_pcnCurrentSymbol, unsigned int* p_pnCurrentAndNextColumnPtr, unsigned int const* p_pcnNextSymbolTablePtr, unsigned int const* p_pcnNextStateTablePtr, unsigned int const p_cnColumnSize, unsigned int const p_cnNibblesInDataType, unsigned int const p_cnSymbolClusterSize, unsigned int const p_cnStartState, unsigned int const p_cnJump)
{
	// First thing, we take the current cell according to the current thread
	unsigned int l_nI = blockDim.x * blockIdx.x + threadIdx.x + p_cnStartState;


	if(l_nI < p_cnColumnSize) // If the cell index doesn't exit from table limits
	{

		if(l_nI / p_cnJump == l_nI / (double) p_cnJump)
		{
			// ANOTHER FILTER FOR SYNCHRONIZATION
			unsigned int l_nInput = (l_nI % 2 == 1) ? 0 : 1;

			// Get the possible states and symbols, from the current active state
			unsigned int const* l_pcnNextPossibleStates = p_pcnNextStateTablePtr + l_nI * 2;
			unsigned int const* l_pccNextPossibleSymbols = p_pcnNextSymbolTablePtr + l_nI * 2 * p_cnSymbolClusterSize;

			// Get the distance of the state from which we got to the current one
			unsigned int l_nPreviousDistanceOfInterest = p_pnCurrentAndNextColumnPtr[l_nI];

			// Get current possible state for the given input
			unsigned int l_nNextPossibleState = l_pcnNextPossibleStates[l_nInput];

			// 1. Calculate the hamming distance between the current symbol and the predicted symbol,
			//	  plus the distance calculated in the previous iteration
			// 2. Get the distance from the next possible state (in the next iteration. It could already
			//	  calculated from other state in this iteration)
			// 3. Take the minimum between those two values above

			unsigned int l_nDistance = p_pnCurrentAndNextColumnPtr[l_nNextPossibleState + p_cnColumnSize];
			unsigned int l_nCalculatedDistance = CudaHammingDistance(p_pcnCurrentSymbol, l_pccNextPossibleSymbols + l_nInput * p_cnSymbolClusterSize, p_cnNibblesInDataType, p_cnSymbolClusterSize) + l_nPreviousDistanceOfInterest;
			if(l_nCalculatedDistance < l_nDistance)
			{
				l_nDistance = l_nCalculatedDistance;
			}



			// Update the distances matrix with the new distance
			p_pnCurrentAndNextColumnPtr[l_nNextPossibleState + p_cnColumnSize] = l_nDistance;

			l_nInput = (l_nInput == 0) ? 1 : 0;

//			__syncthreads();

			// Get current possible state for the given input
			l_nNextPossibleState = l_pcnNextPossibleStates[l_nInput];

			// 1. Calculate the hamming distance between the current symbol and the predicted symbol,
			//	  plus the distance calculated in the previous iteration
			// 2. Get the distance from the next possible state (in the next iteration. It could already
			//	  calculated from other state in this iteration)
			// 3. Take the minimum between those two values above

			l_nDistance = p_pnCurrentAndNextColumnPtr[l_nNextPossibleState + p_cnColumnSize];
			l_nCalculatedDistance = CudaHammingDistance(p_pcnCurrentSymbol, l_pccNextPossibleSymbols + l_nInput * p_cnSymbolClusterSize, p_cnNibblesInDataType, p_cnSymbolClusterSize) + l_nPreviousDistanceOfInterest;
			if(l_nCalculatedDistance < l_nDistance)
			{
				l_nDistance = l_nCalculatedDistance;
			}



			// Update the distances matrix with the new distance
			p_pnCurrentAndNextColumnPtr[l_nNextPossibleState + p_cnColumnSize] = l_nDistance;
		}
	}
}

#endif
