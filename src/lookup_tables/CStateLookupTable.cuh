/**
 * @file CStateLookupTable.h
 * @brief Class definition for abstract states lookup table
 * @author Oren Cohen
 * @bug No known bugs
 */

#ifndef LOOKUP_TABLES_CSTATELOOKUPTABLE_H_
#define LOOKUP_TABLES_CSTATELOOKUPTABLE_H_

#include <iostream>
#include <vector>

#ifdef _USE_CUDA_
#include <cuda_runtime.h> // For CUDA API
#endif

/**
 * Abstract class for getting some predicted state from lookup table
 *
 * Usage example:
 * @code
 * See inherited classes
 * @endcode
 */
class CStateLookupTable
{
public:
	/**
	 * Destructor
	 */
	virtual ~CStateLookupTable(void);

	/**
	 * Get sequential state form table
	 * @param p_nCurrentState The state before the sequential one
	 * @param p_nInput The bit value for the sequential state
	 * @return The sequential state
	 */
	unsigned int GetSequentialState(unsigned int const p_cnCurrentState, unsigned int const p_cnInput) const;
	/**
	 * Get sequential states (possible two) form table
	 * @param p_nCurrentState The state before the sequential states
	 * @return The sequential states
	 */
	unsigned int const* GetSequentialStates(unsigned int const p_cnCurrentState) const;

	/**
	 * Get the maximal state possible according to memory size
	 * @param p_nMemory Given memory size
	 * @return The maximal state
	 */
	static unsigned int GetMaxState(unsigned int const p_cnMemory);

	/**
	 * Get the table size
	 * @return The table size
	 */
	unsigned int GetSize(void) const;

	/**
	 * Get the current table as STL vector
	 * @return The current lookup table
	 */
	std::vector<unsigned int> const& AsVector(void) const;

	/**
	 * For Debug purposes - overload the << operator for this class
	 * @param p_rcOs The output stream
	 * @param p_rcTable The states lookup table for output
	 * @return The current output stream
	 */
	friend std::ostream& operator<<(std::ostream& p_rcOs, CStateLookupTable const& p_rccTable);

	/**
	 * Abstract function for shifting to next state
	 * @param p_nState The current state
	 * @param p_nInput The given bit according to it the shifting performed
	 * @return The new (shifted) state
	 */
	virtual unsigned int Shift(unsigned int const p_cnState, unsigned int const p_cnInput) const = 0;

protected:

#ifdef _USE_CUDA_

	/**
	 * When using CUDA, lookup table type identifier is needed
	 * This member meets the need
	 */
	enum LookupTableType
	{
		LookupTableType_Next,
		LookupTableType_Prev
	} const m_ceTableType;

#endif

	/**
	 * Initialize (fill) the states lookup table
	 */
	void InitializeStateLookupTable(void);
	/**
	 * Protected parameterized constructor that creates the table \n
	 * The constructor is protected because of design purposes \n
	 * The actual creation is trough the Create method in the inherited classes \n
	 * @param p_nK The memory size of the table
	 * @param p_eTableType If using CUDA : the table type
	 */
	CStateLookupTable(unsigned int const p_cnK
#ifdef _USE_CUDA_
			, LookupTableType const p_ceTableType
#endif
			);

	/**
	 * The memory size of the lookup table
	 */
	unsigned int const m_cnMemory;

private:

#ifdef _USE_CUDA_

	/**
	 * Get the GPU information
	 */
	void SetGpuInfo(void);

#endif

	/**
	 * The lookup table values (a matrix)
	 */
	std::vector<unsigned int> m_cTable;

#ifdef _USE_CUDA_

	/**
	 * CUDA: Holds the properties of the machine GPU
	 */
	cudaDeviceProp m_sGpuProps;



#endif

};

#ifdef _USE_CUDA_

/**
 * Kernel function that calculates on row in the lookup table \n
 * In case that the table is too big for the GPU global memory, it will be done by clusters
 * @param p_pnTable Pointer to the lookup table that copied to the GPU
 * @param p_nTableSize The size of the given table
 * @param p_nStartState The first state which the calculation needs to start from
 */

/**
 * Device function that implements next state shifting for next state lookup table
 * @param p_nState The current state
 * @param p_nInput The given bit input (0 or 1)
 * @param p_nMemory The shift register memory size
 * @return The next suitable state according to the given parameters
 */
__device__ unsigned int ShiftNext(unsigned int const p_cnState, unsigned int const p_cnInput, unsigned int const p_cnMemory);
/**
 * Device function that implements previous state shifting for previous state lookup table
 * @param p_nState The current state
 * @param p_nInput The given bit input (0 or 1)
 * @param p_nMemory The shift maximal state possible
 * @return The previous suitable state according to the given parameters
 */
__device__ unsigned int ShiftPrev(unsigned int const p_cnState, unsigned int const p_cnInput, unsigned int const p_cnMaxState);

/**
 * Device function that calculates one cell in some states lookup table
 * @param p_fpShiftFunc Function pointer to the appropriate shifting function
 * @param p_pnTable Pointer to the table allocated memory area in the GPU global memory
 * @param p_nTableSize The table allocated memory size
 * @param p_nStartState In case that the calculation is done by clusters, the first state of the cluster is needed
 * @param p_nAdditionalInfo Additional info that the shifting function needs (memory or maximal state, depends on the table type)
 */
__device__ void CalculateCell(unsigned int (*p_fpShiftFunc)(unsigned int const, unsigned int const, unsigned int const), unsigned int* p_pnTable, unsigned int const p_cnTableSize, unsigned int const p_cnStartState, unsigned int const p_cnAdditionalInfo);

/**
 * Kernel function that calculates one states lookup table cell
 * @param p_eTableType The table type (next/previous state)
 * @param p_nMemory The table memory size
 * @param p_nMaxState The table maximal state
 * @param p_pnTable Pointer to the table on the global memory
 * @param p_nTableSize The table size
 * @param p_nStartState The cluster first state
 */
__global__ void Kernel_Calc(int const p_ceTableType, unsigned int const p_cnMemory, unsigned int const p_cnMaxState, unsigned int* p_pnTable, unsigned int const p_cnTableSize, unsigned int const p_cnStartState);

#endif

#endif /* LOOKUP_TABLES_CSTATELOOKUPTABLE_H_ */
