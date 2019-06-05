/**
 * @file CConvolutionEncoderKToN.h
 * @brief Class definition for convolutional encoder
 * @author Oren Cohen
 * @bug No known bugs
 */

#ifndef CONVOLUTION_ENCODING_CCONVOLUTIONENCODERKTON_H_
#define CONVOLUTION_ENCODING_CCONVOLUTIONENCODERKTON_H_

#include <iostream>
#include <vector>

#include "../lookup_tables/CNextStateLookupTable.cuh"
#include "../lookup_tables/CNextEncodedSymbolLookupTable.h"
#include "CShiftRegister.h"

/**
 * Convolutional encoding for some data type
 */
template<typename T>
class CConvolutionEncoderKToN
{
public:
	/**
	 * Parameterized constructor which initialize the encoder by specific parameters
	 * @param p_pcNextStateLookupTable Lookup table for getting next state
	 * @param p_pcNextEncodedSymbolLookupTable Lookup table for getting next possible symbol by current state
	 * @param p_nK The memory size of the encoder
	 * @param p_nN The symbol size
	 */
	CConvolutionEncoderKToN(CNextStateLookupTable const& p_rccNextStateLookupTable,
			CNextEncodedSymbolLookupTable const& p_rccNextEncodedSymbolLookupTable, unsigned int const p_cnK, unsigned int const p_cnN);
	/**
	 * Destructor
	 */
	virtual ~CConvolutionEncoderKToN(void);

	/**
	 * Encode the given data
	 * @param p_pnData A given data buffer of some type
	 * @param p_nSize Number of elements in the given buffer
	 */
	void Convolove(std::vector<T> const& p_pnData);
	/**
	 * Get the encoded symbols
	 * @return Pointer to encoded symbols buffer
	 */
	std::vector<std::list<unsigned int>* > const& GetSymbols(void) const;
	/**
	 * Get number of encoded symbols
	 * @return Number of encoded symbols
	 */
	unsigned int GetNumberOfSymbols(void) const;

	/**
	 * For Debug purposes - overload the << operator for this class
	 * @param p_rcOs The output stream
	 * @param p_rcEncoder The encoder for output
	 * @return The current output stream
	 */
	template<typename SAME_T>
	friend std::ostream& operator<<(std::ostream& p_rcOs, CConvolutionEncoderKToN<SAME_T> const& p_rccEncoder);

private:
	/**
	 * The shift register that participates with the encoding
	 */
	CShiftRegister m_cShiftRegister;

	/**
	 * Encoded symbols buffer
	 */
	std::vector<std::list<unsigned int>* > m_cSymbols;

	/**
	 * Encoder memory size
	 */
	unsigned int m_nMemory;
	/**
	 * Length in bit on encoded symbol
	 */
	unsigned int m_nSymbolLength;
};

#include "CConvolutionEncoderKToN.hpp"

#endif /* CONVOLUTION_ENCODING_CCONVOLUTIONENCODERKTON_H_ */
