
#include "CViterbiDecoder.h"

#include <stddef.h>
#include <iostream>

#include <list>
#include <vector>

#ifdef _USE_MPI_
#include <mpi/mpi.h>
#include <string.h>
#endif

#include "../lookup_tables/CPreviousStateLookupTable.cuh"
#include "CTrellisDiagram.cuh"

template<typename T>
CViterbiDecoder<T>::CViterbiDecoder(unsigned int const p_cnK, CTrellisDiagram const& p_rccTrellisDiagram,
		CPreviousStateLookupTable const& p_rccPreviousStateLookupTable) :
		m_rccTrellisDiagram(p_rccTrellisDiagram), m_rccPreviousStateLookupTable(p_rccPreviousStateLookupTable),
		m_nErrorBitsCount(0), m_cnSparesNumber(p_cnK)
{

#ifdef _USE_MPI_
	// Setup MPI
	InitializeMpi();
#endif

	// initialize buffers for backtrace, decoded data, etc.
	InitializeBuffers();
}

template<typename T>
CViterbiDecoder<T>::~CViterbiDecoder(void)
{
#ifdef _USE_MPI_
	int l_nMpiFinilized;
	MPI_Finalized(&l_nMpiFinilized);
	if(!l_nMpiFinilized)
	{
		MPI_Finalize();
	}
#endif
}

template<typename T>
void CViterbiDecoder<T>::InitializeBuffers(void)
{
	// Get the number of data bits
	unsigned int const l_cnNumberOfDataBits = m_rccTrellisDiagram.GetColumnsCount() - 1 - m_cnSparesNumber;

	// Initialize the backtrace size and memory
	m_cBacktraceRoad.resize(l_cnNumberOfDataBits + m_cnSparesNumber + 1);

	// Initialize the decoded data buffer
	m_cDecodedData.resize(l_cnNumberOfDataBits / (sizeof(T) * 8)); // total_number_of_bits / data_size_in_bits
}

template<typename T>
std::vector<T> const& CViterbiDecoder<T>::GetDecodedData(void) const
{
	return m_cDecodedData;
}

template<typename T>
unsigned int CViterbiDecoder<T>::GetDecodedDataSize(void) const
{
	return m_cDecodedData.size();
}

template<typename T>
std::vector<unsigned int> const& CViterbiDecoder<T>::GetBacktraceRoad(void) const
{
	return m_cBacktraceRoad;
}

template<typename T>
unsigned int CViterbiDecoder<T>::GetBacktraceRoadSize(void) const
{
	return m_cBacktraceRoad.size();
}

template<typename T>
unsigned int CViterbiDecoder<T>::GetErrorBitsCount(void) const
{
	return m_nErrorBitsCount;
}

template<typename T>
void CViterbiDecoder<T>::Backtrace(void)
{
	m_nErrorBitsCount = 0;

	unsigned int l_nBitsInDataItemCount = 0;
	T l_nDecodedDataItem = 0;
	unsigned int l_nTotalBitsCount = 0;

	unsigned int l_nDataCount = 0;
	unsigned int l_nTransitionsCount = 0;

	std::list<std::vector<unsigned int> > const l_ppnDistancesMatrix = m_rccTrellisDiagram.GetDistancesMatrix();

	// Find the minimum distance from the last column of the trellis diagram,
	// it is also the number of error bits received after decoding and the first
	// value of the backtrace
	unsigned int l_nLastStageIndex = m_rccTrellisDiagram.GetColumnsCount() - 1;

	std::list<std::vector<unsigned int> >::const_reverse_iterator l_cLastColumnIt = l_ppnDistancesMatrix.rbegin();
	unsigned int const* l_pnLastColumn = &(*l_cLastColumnIt)[0];
	unsigned int l_nMinDistance = l_pnLastColumn[0];

	unsigned int l_nSuitableState = 0;
	for (unsigned int l_nCurrentState = 1; l_nCurrentState < m_rccTrellisDiagram.GetRowsCount(); l_nCurrentState++)
	{
		if (l_pnLastColumn[l_nCurrentState] < l_nMinDistance)
		{
			l_nMinDistance = l_pnLastColumn[l_nCurrentState];
			l_nSuitableState = l_nCurrentState;
		}
	}
	m_nErrorBitsCount = l_nMinDistance; // number of error bits
	m_cBacktraceRoad[l_nTransitionsCount++] = l_nSuitableState; // first backtrace value

	l_cLastColumnIt++;

	// get the possible previous states according to the state of the minimum found above
	unsigned int const* l_pnPossiblePredecessorsStates = m_rccPreviousStateLookupTable.GetSequentialStates(l_nSuitableState);

	// loop over all the columns of the trellis diagram, except from the spares (represents the extra
	// zero values from the encoder)
	for (unsigned int l_nCurrentState = l_nLastStageIndex - 1; // one column before the last
			// don't consider the spare bits
			l_nCurrentState >= 2 && l_nTotalBitsCount <= (l_nLastStageIndex + 1) - m_cnSparesNumber;
			l_nCurrentState--)
	{
		unsigned int const* l_pnCurrentColumn =  &(*l_cLastColumnIt)[0];

		// take the predecessor with the minimal value and decode its bit value to the data item
		if (l_pnCurrentColumn[l_pnPossiblePredecessorsStates[0]]
				> l_pnCurrentColumn[l_pnPossiblePredecessorsStates[1]])
		{
			l_nSuitableState = l_pnPossiblePredecessorsStates[1]; // save state for next iteration
			l_nDecodedDataItem += (1 << ((sizeof(T) * 8 - l_nBitsInDataItemCount) - 1)); // decode bit
		}
		else
		{
			l_nSuitableState = l_pnPossiblePredecessorsStates[0]; // save state for next iteration
		}

		// take previous possible states from the updated current state, and update the backtrace
		l_pnPossiblePredecessorsStates = m_rccPreviousStateLookupTable.GetSequentialStates(l_nSuitableState);
		m_cBacktraceRoad[l_nTransitionsCount++] = l_nSuitableState;

		l_nBitsInDataItemCount++;
		l_nTotalBitsCount++;

		// when one item decoding finished, put it in the buffer
		if(l_nBitsInDataItemCount == (sizeof(T) * 8))
		{
			l_nBitsInDataItemCount = 0;
			m_cDecodedData[l_nDataCount++] = l_nDecodedDataItem;
			l_nDecodedDataItem = 0;
		}

		l_cLastColumnIt++;
	}
}

template<typename T>
void CViterbiDecoder<T>::Decode(void)
{
	// Perform backtrace and also decode the data on the way
	Backtrace();

#ifdef _USE_MPI_

	if(m_bImMaster)
	{
		// Reallocate the space for the decoded data
		unsigned int l_nPrevDataSize = m_cDecodedData.size();

		// Reallocate the memory of the master and copy to this memory its decoded data
		m_cDecodedData.resize(m_cDecodedData.size() * m_nNumberOfProcesses); // Set master size according the number of processes
		memcpy(&m_cDecodedData[0] + (m_nNumberOfProcesses - 1) * l_nPrevDataSize, &m_cDecodedData[0], l_nPrevDataSize * sizeof(T));

		// Receive the rest of the decoded data from the rest of the processes
		// and copy it to the master memory
		int l_nRecvCtr = 0;
		MPI_Status l_sStatus;
		T* l_pnRecvBuffer = new T[l_nPrevDataSize];
		while(l_nRecvCtr++ < m_nNumberOfProcesses - 1)
		{
			// Receive data from some other process
			MPI_Recv(l_pnRecvBuffer, l_nPrevDataSize * sizeof(T), MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &l_sStatus);
			memcpy(&m_cDecodedData[0] + (m_nNumberOfProcesses - 1 - l_sStatus.MPI_SOURCE) * l_nPrevDataSize, l_pnRecvBuffer, l_nPrevDataSize * sizeof(T));
		}


		delete[] l_pnRecvBuffer; // Release reception memory area
	}
	// The slaves send their decoded data to the master
	else
	{
		// Send all the data as bytes because the template variable
		MPI_Send(&m_cDecodedData[0], m_cDecodedData.size() * sizeof(T), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	}

#endif
}

template<typename T>
std::ostream& operator<<(std::ostream& p_rcOs, CViterbiDecoder<T> const& p_rccDecoder)
{
	p_rcOs << "viterbi decoded data:" << std::endl;

	p_rcOs << "backtrace: ";
	for(unsigned int l_nI = 0; l_nI < p_rccDecoder.m_cBacktraceRoad.size(); l_nI++)
	{
		p_rcOs << p_rccDecoder.m_cBacktraceRoad[l_nI] << " ";
	}

	p_rcOs << std::endl << "decoded data: ";
	for(unsigned int l_nI = 0; l_nI < p_rccDecoder.m_cDecodedData.size(); l_nI++)
	{
		p_rcOs << p_rccDecoder.m_cDecodedData[l_nI] << " ";
	}

	p_rcOs << std::endl << "number of bit errors: " << p_rccDecoder.m_nErrorBitsCount << std::endl;

    return p_rcOs;
}

#ifdef _USE_MPI_
template<typename T>
void CViterbiDecoder<T>::InitializeMpi(void)
{
	// Before initializing MPI, check if it already initialized
	int l_nMpiInitialized;
	MPI_Initialized(&l_nMpiInitialized);
	if(!l_nMpiInitialized)
	{
		MPI_Init(NULL, NULL);
	}

	// Get current process rank and number of runnning processes
	MPI_Comm_rank(MPI_COMM_WORLD, &m_nMyRank);
	m_bImMaster = (0 == m_nMyRank) ? true : false;

	MPI_Comm_size(MPI_COMM_WORLD, &m_nNumberOfProcesses);

}
#endif
