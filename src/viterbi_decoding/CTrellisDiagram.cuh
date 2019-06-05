/**
 * @file CTrellisDiagram.h
 * @brief Class definition of trellis diagram for convolutional encoding
 * @author Oren Cohen
 * @bug No known bugs
 */

#ifndef VITERBI_DECODING_CTRELLISDIAGRAM_H_
#define VITERBI_DECODING_CTRELLISDIAGRAM_H_

#include <iostream>

#include <list>
#include <vector>

#ifdef _USE_CUDA_
#include <cuda_runtime.h> // For CUDA API
#endif

#include "../lookup_tables/CNextStateLookupTable.cuh"
#include "../lookup_tables/CNextEncodedSymbolLookupTable.h"



/**
 * Calculation of trellis diagram for convolutional encoding
 */
class CTrellisDiagram
{
public:

	/**
	 * Parameterized constructor
	 * @param p_nNumberOfSymbols Number of encoded symbols (column for each symbol)
	 * @param p_nK Encoding memory size
	 * @param p_pcNextStateLookupTable Lookup table for next predicted state
	 * @param p_pcNextEncodedSymbolLookupTable Next possible symbol lookup table
	 * @param p_nSizeOfDecodedDataType The size in bytes of the decoded data type (for example: sizeof operator of the type)
	 */
	CTrellisDiagram(unsigned int const p_cnNumberOfSymbols, unsigned int const p_cnK, CNextStateLookupTable const& p_rccNextStateLookupTable,
			CNextEncodedSymbolLookupTable const& p_rccNextEncodedSymbolLookupTable, unsigned int const p_cnSizeOfDecodedDataType);
	/**
	 * Destructor
	 */
	virtual ~CTrellisDiagram(void);

	/**
	 * Get the distances matrix
	 * @return The distances matrix of the trellis diagram
	 */
	std::list<std::vector<unsigned int> > const& GetDistancesMatrix(void) const;

	/**
	 * Calculate trellis diagram
	 * @param p_pnSymbols The decoded symbol
	 */
	void PerformStages(std::vector<std::list<unsigned int>* > const& p_rccSymbols);

	/**
	 * Get number of rows
	 * @return Number of rows in the diagram
	 */
	unsigned int GetRowsCount(void) const;
	/**
	 * Get number of columns
	 * @return Number of columns in the diagram
	 */
	unsigned int GetColumnsCount(void) const;

	/**
	 * For Debug purposes - overload the << operator for this class
	 * @param p_rcOs The output stream
	 * @param p_rcTrellis The diagram for output
	 * @return The current output stream
	 */
	friend std::ostream& operator<<(std::ostream& p_rcOs, CTrellisDiagram const& p_rccTrellis);

private:

	/**
	 * Initialize the diagram parameters
	 */
	void InitializeTrellisDiagram(void);
	/**
	 * Clear the diagram (for new decoding)
	 */
	void ClearTrellisDiagram(void);

	/**
	 * Initialize the distances matrix of the diagram (initial value = INFINITY)
	 */
	void ClearDistancesMatrix(void);

	/**
	 * Perform the first iteration of the diagram (it slight different from \n
	 * the next calculations)
	 * @param p_nSymbol Given encoded symbol
	 */
	void PerformFirstStages(std::list<unsigned int> const& p_rccSymbol, std::list<std::vector<unsigned int> >::iterator& p_rcTrellisItPtr);
	/**
	 * Perform one sequential iteration of the diagram
	 * @param p_nSymbol Given encoded symbol
	 */
	void PerformSequentialStage(std::list<unsigned int> const& p_rccSymbol, std::list<std::vector<unsigned int> >::iterator& p_rcTrellisItPtr);

	/**
	 * Get the hamming distance between two symbols \n
	 * Hamming distance is the number of different parallel bits between two binary strings
	 * @param p_nFirstSymbol The first symbol
	 * @param p_nSecondSymbol The second symbol
	 * @param p_nNibblesInDataType The number of nibbles in the given data type (for example: 2 for char)
	 * @return The hamming distance between the two given numbers
	 */
	static unsigned int GetHammingDistance(unsigned int const p_cnFirstSymbol, unsigned int const p_cnSecondSymbol,
			unsigned int const p_cnNibblesInDataType);
	static unsigned int GetHammingDistance(std::list<unsigned int> const& p_rccFirstSymbol, std::list<unsigned int> const& p_rccSecondSymbol,
				unsigned int const p_cnNibblesInDataType);

#ifdef _USE_MPI_
	/**
	 * Get the relevant parameters and configurations for MPI compilation
	 */
	void InitializeMpi(void);
#endif

#ifdef _USE_CUDA_

	/**
	 * Get the GPU information
	 */
	void SetGpuInfo(void);

	/**
	 * Copy the needed lookup tables to the device in the initial stage of this object construction
	 */
	void AllocateAndSetDeviceMemory(void);

	/**
	 * Reset the GPU and free the memory allocated on it
	 */
	void FreeMemoryAndResetDevice(void);

#endif

	/**
	 * Constant representation for used infinity
	 */
	static const unsigned int s_cnInfinity;
	/**
	 * Map that tells the number of set bits in every nibble for the numbers 0...15
	 */
	static const unsigned int s_pcnSetBitsInNibble[];

	/**
	 * The distances matrix
	 */
	std::list<std::vector<unsigned int> > m_cDistancesMatrix;
	/**
	 * Number of rows in distances matrix
	 */
	unsigned int m_nRowsCount;
	/**
	 * Number of columns in distances matrix
	 */
	unsigned int m_nColumnsCount;

	/**
	 * Number of encoded symbols
	 */
	unsigned int m_nNumberOfSymbols;

	/**
	 * Encoder memory size
	 */
	unsigned int const m_cnMemory;

	/**
	 * Value of maximal state possible
	 */
	unsigned int m_nMaxState;

	/**
	 * Next possible states lookup table
	 */
	CNextStateLookupTable const& m_rccNextStateLookupTable;
	/**
	 * Next possible symbols lookup table
	 */
	CNextEncodedSymbolLookupTable const& m_rccNextEncodedSymbolLookupTable;

	/**
	 * Number of nibbles in data type (for no returned calculations)
	 */
	unsigned int const m_cnNibblesInDataType;

	unsigned int m_nJump;

#ifdef _USE_MPI_

	/**
	 * Parallel : current process rank
	 */
	int m_nMyRank;

	/**
	 * Parallel: number of processes
	 */
	int m_nNumberOfProcesses;

	/**
	 * Parallel: indication if the current process is the master
	 */
	bool m_bImMaster;
#endif

#ifdef _USE_CUDA_

	/**
	 * CUDA: Holds the properties of the machine GPU
	 */
	cudaDeviceProp m_sGpuProps;

	/**
	 * Pointer to the memory allocated on the device for two columns of the diagram
	 */
	unsigned int* m_pnCurrentAndNextColumnDevicePtr;

	/**
	 * Pointer to the memory allocated on the device for the next symbol lookup table
	 */
	unsigned int* m_pnNextSymbolTableDevicePtr;

	/**
	 * Pointer to next state lookup table on the device
	 */
	unsigned int* m_pnNextStateTableDevicePtr;

	/**
	 * Pointer to the symbol of the current trellis iteration
	 */
	unsigned int* m_pnIterationSymbolDevicePtr;

#endif

};

#ifdef _USE_CUDA_

__device__ unsigned int CudaHammingDistance(unsigned int const* p_pcnFirstSymbol, unsigned int const* p_pcnSecondSymbol, unsigned int const p_cnNibblesInDataType, unsigned int const p_cnNumberOfClusters);

__global__ void TrellisIteration_Kernel(unsigned int const* p_pcnCurrentSymbol, unsigned int* p_pnCurrentAndNextColumnPtr, unsigned int const* p_pcnNextSymbolTablePtr, unsigned int const* p_pcnNextStateTablePtr, unsigned int const p_cnColumnSize, unsigned int const p_cnNibblesInDataType, unsigned int const p_cnSymbolClusterSize, unsigned int const p_cnStartState, unsigned int const p_cnJump);

#endif

#endif /* VITERBI_DECODING_CTRELLISDIAGRAM_H_ */
