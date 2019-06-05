
#include "CStateLookupTable.cuh"

#include <stddef.h>
#include <math.h>
#include <iostream>

#include <iomanip>
#include <algorithm>
#include <vector>

#ifdef _USE_CUDA_
#include <cuda_runtime.h> // For CUDA API
#include <stdio.h>
#endif

#include "CPolynomialGenerator.h"

CStateLookupTable::CStateLookupTable(unsigned int const p_nK
#ifdef _USE_CUDA_
		, LookupTableType const p_eTableType // Table type parameter
#endif
		) :
		m_cnMemory(p_nK)
#ifdef _USE_CUDA_
, m_ceTableType(p_eTableType) // Set table type
#endif
{
#ifdef _USE_CUDA_
	SetGpuInfo(); // Get GPU properties
#endif
}

CStateLookupTable::~CStateLookupTable(void)
{
#ifdef _USE_CUDA_
	cudaDeviceReset(); // Reset the device at object destruction
#endif
}

unsigned int CStateLookupTable::GetSequentialState(unsigned int const p_cnCurrentState, unsigned int const p_cnInput) const
{
	return GetSequentialStates(p_cnCurrentState)[p_cnInput];
}

unsigned int const* CStateLookupTable::GetSequentialStates(unsigned int const p_cnCurrentState) const
{
	// Get sequential state according to given state
	return &m_cTable[p_cnCurrentState * 2];
}

unsigned int CStateLookupTable::GetMaxState(unsigned int const p_cnMemory)
{
	// Return the maximal state possible ( 2^K - 1)
	return pow(2, p_cnMemory) - 1;
}

unsigned int CStateLookupTable::GetSize(void) const
{
	return m_cTable.size();
}

std::vector<unsigned int> const& CStateLookupTable::AsVector(void) const
{
	return m_cTable;
}

void CStateLookupTable::InitializeStateLookupTable(void)
{
#ifdef _USE_CUDA_

	// First, allocate the memory for the host table, and set its size
	unsigned int const l_cnMaxState = GetMaxState(m_cnMemory);
	m_cTable.resize((l_cnMaxState + 1) * 2);

	// Get maximum number of threads per block and calculate the number of blocks
	// that we need
	unsigned int const l_cnNumberOfThreadsPerBlock = m_sGpuProps.maxThreadsPerBlock;
	unsigned int l_nNumberOfBlocks = (m_cTable.size() / l_cnNumberOfThreadsPerBlock == 0) ? 1 : m_cTable.size() / l_cnNumberOfThreadsPerBlock;

	// Get maximum number of blocks according to the GPU properties (only one dimension needed)
	const unsigned int l_cnMaxNumberOfBlocks = m_sGpuProps.maxThreadsDim[0];

	// Pointer to device memory
	unsigned int* l_pnDeviceBuffer = NULL;

	// If all the threads can be run at once
	if(l_nNumberOfBlocks <= l_cnMaxNumberOfBlocks)
	{
		cudaMalloc((void**)&l_pnDeviceBuffer, m_cTable.size() * sizeof(unsigned int));

		Kernel_Calc<<<l_nNumberOfBlocks, l_cnNumberOfThreadsPerBlock>>>(m_ceTableType, m_cnMemory, l_cnMaxState, l_pnDeviceBuffer, m_cTable.size(), 0);

		cudaMemcpy(&m_cTable[0], l_pnDeviceBuffer, m_cTable.size() * sizeof(unsigned int), cudaMemcpyDeviceToHost);

		cudaFree(l_pnDeviceBuffer);
	}
	else
	{
		cudaMalloc((void**)&l_pnDeviceBuffer, m_cTable.size() * sizeof(unsigned int));

		unsigned int l_nCurrentState = 0;
		while(l_nNumberOfBlocks >= l_cnMaxNumberOfBlocks)
		{
			Kernel_Calc<<<l_cnMaxNumberOfBlocks, l_cnNumberOfThreadsPerBlock>>>(m_ceTableType, m_cnMemory, l_cnMaxState, l_pnDeviceBuffer, m_cTable.size(), l_nCurrentState);
			cudaMemcpy(&m_cTable[l_nCurrentState], l_pnDeviceBuffer, l_cnMaxNumberOfBlocks * l_cnNumberOfThreadsPerBlock * sizeof(unsigned int), cudaMemcpyDeviceToHost);
			l_nNumberOfBlocks -= l_cnMaxNumberOfBlocks;
			l_nCurrentState += l_cnMaxNumberOfBlocks * l_cnNumberOfThreadsPerBlock;
		}
		if(l_nNumberOfBlocks > 0)
		{
			Kernel_Calc<<<l_nNumberOfBlocks, l_cnNumberOfThreadsPerBlock>>>(m_ceTableType, m_cnMemory, l_cnMaxState, l_pnDeviceBuffer, m_cTable.size(), l_nCurrentState);
			cudaMemcpy(&m_cTable[l_nCurrentState], l_pnDeviceBuffer, (m_cTable.size() - l_nCurrentState) * sizeof(unsigned int), cudaMemcpyDeviceToHost);
		}

		cudaFree(l_pnDeviceBuffer);
	}
	cudaDeviceReset();

#else
	m_cTable.resize((GetMaxState(m_cnMemory) + 1) * 2);
	unsigned int* l_pnTablePtr = &m_cTable[0];
	for (unsigned int l_nCurrentState = 0; l_nCurrentState < m_cTable.size() / 2; l_nCurrentState++)
	{
		// For each state possible set its sequential states according
		// to shifting with 0 or 1
		l_pnTablePtr[0] = Shift(l_nCurrentState, 0);
		l_pnTablePtr[1] = Shift(l_nCurrentState, 1);
		l_pnTablePtr += 2;
	}
#endif
}

std::ostream& operator<<(std::ostream& p_rcOs, CStateLookupTable const& p_rccTable)
{

	p_rcOs << "state lookup table:" << std::endl;

	unsigned int const l_cnStateWidth = p_rccTable.m_cnMemory;

	unsigned int const l_cnStateColumnWidth = std::max(5, (int)l_cnStateWidth);
	unsigned int const l_cnInputColumnWidth = std::max(8, (int)l_cnStateWidth);

	unsigned int const l_cnDelimeterWidth = std::max(23, (int)(l_cnInputColumnWidth * 2 + l_cnStateColumnWidth + 2));

	p_rcOs << '+';
	for(unsigned int l_nI = 0; l_nI < l_cnDelimeterWidth; l_nI++)
	{
		p_rcOs << '-';
	}
	p_rcOs << '+' << std::endl;

	p_rcOs << '|' << std::setw(l_cnStateColumnWidth) << std::left << "state" << '|'
			<< std::setw(l_cnInputColumnWidth) << std::left << "input: 0" << '|'
			<< std::setw(l_cnInputColumnWidth) << std::left << "input: 1" << '|'
			<< std::endl;

	p_rcOs << '|';
	for(unsigned int l_nI = 0; l_nI < l_cnDelimeterWidth; l_nI++)
	{
		p_rcOs << '-';
	}
	p_rcOs << '|' << std::endl;

	unsigned int l_nCurrentState = 0;
	for(unsigned int l_nRI = 0; l_nRI < p_rccTable.m_cTable.size(); l_nRI+=2)
	{
		p_rcOs << '|' << std::setw(l_cnStateColumnWidth) << std::left << CPolynomialGenerator::GetBinaryString(l_nCurrentState++, p_rccTable.m_cnMemory) << '|'
				<< std::setw(l_cnInputColumnWidth) << std::left << CPolynomialGenerator::GetBinaryString(p_rccTable.m_cTable[l_nRI + 0], p_rccTable.m_cnMemory) << '|'
				<< std::setw(l_cnInputColumnWidth) << std::left << CPolynomialGenerator::GetBinaryString(p_rccTable.m_cTable[l_nRI + 1], p_rccTable.m_cnMemory) << '|'
				<< std::endl;
	}

	p_rcOs << '+';
	for(unsigned int l_nI = 0; l_nI < l_cnDelimeterWidth; l_nI++)
	{
		p_rcOs << '-';
	}
	p_rcOs << '+' << std::endl;

    return p_rcOs;
}

#ifdef _USE_CUDA_

void CStateLookupTable::SetGpuInfo(void)
{
	// Get the device properties into the class member
	cudaGetDeviceProperties((cudaDeviceProp*)&m_sGpuProps, 0);
}

__device__ unsigned int ShiftNext(unsigned int const p_cnState, unsigned int const p_cnInput, unsigned int const p_cnMemory)
{
	return (p_cnState >> 1) + ((p_cnInput == 0) ? (0) : (1 << (p_cnMemory - 1)));
}

__device__ unsigned int ShiftPrev(unsigned int const p_cnState, unsigned int const p_cnInput, unsigned int const p_cnMaxState)
{
	return ((p_cnState << 1) & p_cnMaxState) + p_cnInput;
}

__global__ void Kernel_Calc(int const p_ceTableType, unsigned int const p_cnMemory, unsigned int const p_cnMaxState, unsigned int* p_pnTable, unsigned int const p_cnTableSize, unsigned int const p_cnStartState)
{
	// In case of next state lookup table
	if(p_ceTableType == 0)
	{
		CalculateCell(ShiftNext, p_pnTable, p_cnTableSize, p_cnStartState, p_cnMemory);
	}
	// In case of previous state lookup table
	else if(p_ceTableType == 1)
	{
		CalculateCell(ShiftPrev, p_pnTable, p_cnTableSize, p_cnStartState, p_cnMaxState);
	}
}

__device__ void CalculateCell(unsigned int (*p_fpShiftFunc)(unsigned int const, unsigned int const, unsigned int const), unsigned int* p_pnTable, unsigned int const p_cnTableSize, unsigned int const p_cnStartState, unsigned int const p_cnAdditionalInfo)
{
	// First thing, we take the current cell according to the current thread
	int l_nI = blockDim.x * blockIdx.x + threadIdx.x;

	if(l_nI < p_cnTableSize) // If the cell index doesn't exit from table limits
	{
		int l_nCurrentState = p_cnStartState + l_nI / 2; // Get the current state
		int l_nInput = (l_nI % 2 == 0) ? 0 : 1; // Decide what is the input

		p_pnTable[l_nI] = (*p_fpShiftFunc)(l_nCurrentState, l_nInput, p_cnAdditionalInfo); // Perform the shifting
	}
}

#endif
