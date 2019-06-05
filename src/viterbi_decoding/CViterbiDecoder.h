/**
 * @file CViterbiDecoder.h
 * @brief Class definition of viterbi decoder
 * @author Oren Cohen
 * @bug No known bugs
 */

#ifndef VITERBI_DECODING_CVITERBIDECODER_H_
#define VITERBI_DECODING_CVITERBIDECODER_H_

#include <vector>

#include "../lookup_tables/CPreviousStateLookupTable.cuh"
#include "CTrellisDiagram.cuh"

/**
 * Viterbi decode using trellis diagram
 */
template<typename T>
class CViterbiDecoder
{
public:
	/**
	 * Parameterized constructor for decoder initialization
	 * @param p_nK The decoder memory size
	 * @param p_pcTrellisDiagram Suitable trellis diagram for the decoder
	 * @param p_pcPreviousStateLookupTable Lookup table for traceback procedure
	 */
	CViterbiDecoder(unsigned int const p_cnK, CTrellisDiagram const& p_rccTrellisDiagram,
			CPreviousStateLookupTable const& p_rccPreviousStateLookupTable);
	/**
	 * Destructor
	 */
	virtual ~CViterbiDecoder(void);

	/**
	 * Get the decoded data
	 * @return The decoded data of some data type
	 */
	std::vector<T> const& GetDecodedData(void) const;
	/**
	 * Get decoded data size
	 * @return The number of decoded data elements
	 */
	unsigned int GetDecodedDataSize(void) const;
	/**
	 * Get backtrace
	 * @return The backtrace road of the decoding process
	 */
	std::vector<unsigned int> const& GetBacktraceRoad(void) const;
	/**
	 * Get backtrace size
	 * @return The number of steps in the backtrace
	 */
	unsigned int GetBacktraceRoadSize(void) const;
	/**
	 * Get error bits count
	 * @return Number of error bits discovered between the decoded data and the received data
	 */
	unsigned int GetErrorBitsCount(void) const;

	/**
	 * Decode the received symbols
	 */
	void Decode(void);

	/**
	 * For Debug purposes - overload the << operator for this class
	 * @param p_rcOs The output stream
	 * @param p_rcDecoder The decoder for output
	 * @return The current output stream
	 */
	template<typename SAME_T>
	friend std::ostream& operator<<(std::ostream& p_rcOs, CViterbiDecoder<SAME_T> const& p_rccDecoder);

private:
	/**
	 * Initialize the decoder buffers
	 */
	void InitializeBuffers(void);

	/**
	 * Perform backtrace on the trellis diagram
	 */
	void Backtrace(void);

#ifdef _USE_MPI_
	/**
	 * Configuration of MPI
	 */
	void InitializeMpi(void);
#endif

	/**
	 * The decoder trellis diagram
	 */
	CTrellisDiagram const& m_rccTrellisDiagram;

	/**
	 * The decoder states lookup table for backtrace
	 */
	CPreviousStateLookupTable const& m_rccPreviousStateLookupTable;

	/**
	 * Decoded data
	 */
	std::vector<T> m_cDecodedData;

	/**
	 * Backtrace information
	 */
	std::vector<unsigned int> m_cBacktraceRoad;

	/**
	 * Number of error bits
	 */
	unsigned int m_nErrorBitsCount;

	/**
	 * Number of spares in the decoded data (not really a data)
	 */
	unsigned int const m_cnSparesNumber;

#ifdef _USE_MPI_

	/**
	 * Parallel: number of processes
	 */
	int m_nNumberOfProcesses;

	/**
	 * Parallel : current process rank
	 */
	int m_nMyRank;

	/**
	 * Parallel: indication if the current process is the master
	 */
	bool m_bImMaster;

#endif

};

#include "CViterbiDecoder.hpp"

#endif /* VITERBI_DECODING_CVITERBIDECODER_H_ */
