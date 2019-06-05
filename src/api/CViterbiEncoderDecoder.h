/*
 * CViterbiEncoderDecoder.h
 *
 *  Created on: Feb 6, 2015
 *      Author: orenc
 */

#ifndef CVITERBIENCODERDECODER_H_
#define CVITERBIENCODERDECODER_H_

#include <list>
#include <vector>

#include <iostream>

#include "../lookup_tables/CNextEncodedSymbolLookupTable.h"
#include "../lookup_tables/CNextStateLookupTable.cuh"
#include "../lookup_tables/CPolynomialGenerator.h"
#include "../lookup_tables/CPreviousStateLookupTable.cuh"

#include "../convolution_encoding/CConvolutionEncoderKToN.h"

#include "../viterbi_decoding/CTrellisDiagram.cuh"
#include "../viterbi_decoding/CViterbiDecoder.h"

template<typename T>
class CViterbiEncoderDecoder
{
public:
	CViterbiEncoderDecoder(unsigned int const p_cnK, unsigned int const p_cnN, unsigned int const p_cnSize);
	virtual ~CViterbiEncoderDecoder(void);

	std::vector<std::list<unsigned int>* > const& Encode(std::vector<T> const& p_rccData);
	std::vector<T> const& Decode(std::vector<std::list<unsigned int>* > const& p_rccSymbols);

	template<typename SAME_T>
	friend std::ostream& operator<<(std::ostream& p_rcOs, CViterbiEncoderDecoder<SAME_T> const& p_rccEncoderDecoder);

private:

	unsigned int const m_cnSize;

	CPolynomialGenerator m_cPolynomialGenerator;

	CNextStateLookupTable* m_pcNextStateLookupTable;
	CPreviousStateLookupTable* m_pcPrevStateLookupTable;

	CNextEncodedSymbolLookupTable m_cNextSymbolLookupTable;

	CConvolutionEncoderKToN<T> m_cConvolutionEncoder;

	CTrellisDiagram m_cTrellisDiagram;
	CViterbiDecoder<T> m_cViterbiDecoder;
};

#include "CViterbiEncoderDecoder.hpp"

#endif /* CVITERBIENCODERDECODER_H_ */
