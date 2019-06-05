/*
 * CViterbiEncoderDecoder.cpp
 *
 *  Created on: Feb 6, 2015
 *      Author: orenc
 */

#include "CViterbiEncoderDecoder.h"

#include <stddef.h>
#include <iostream>

#include "../lookup_tables/CNextEncodedSymbolLookupTable.h"
#include "../lookup_tables/CNextStateLookupTable.cuh"
#include "../lookup_tables/CPolynomialGenerator.h"
#include "../lookup_tables/CPreviousStateLookupTable.cuh"

#include "../convolution_encoding/CConvolutionEncoderKToN.h"

#include "../viterbi_decoding/CTrellisDiagram.cuh"
#include "../viterbi_decoding/CViterbiDecoder.h"

template<typename T>
CViterbiEncoderDecoder<T>::CViterbiEncoderDecoder(unsigned int const p_cnK,
		unsigned int const p_cnN, unsigned int const p_cnSize) :
		m_cnSize(p_cnSize),
		m_cPolynomialGenerator(p_cnN, p_cnK),
		m_pcNextStateLookupTable(CNextStateLookupTable::Create(p_cnK)),
		m_pcPrevStateLookupTable(CPreviousStateLookupTable::Create(p_cnK)),
		m_cNextSymbolLookupTable(m_cPolynomialGenerator),
		m_cConvolutionEncoder(*m_pcNextStateLookupTable, m_cNextSymbolLookupTable, p_cnK, p_cnN),
		m_cTrellisDiagram(p_cnSize * sizeof(T) * 8 + p_cnK, p_cnK, *m_pcNextStateLookupTable, m_cNextSymbolLookupTable, sizeof(T)),
		m_cViterbiDecoder(p_cnK, m_cTrellisDiagram, *m_pcPrevStateLookupTable)
{
}

template<typename T>
CViterbiEncoderDecoder<T>::~CViterbiEncoderDecoder(void)
{
	if (NULL != m_pcNextStateLookupTable)
	{
		delete m_pcNextStateLookupTable;
	}
	if (NULL != m_pcPrevStateLookupTable)
	{
		delete m_pcPrevStateLookupTable;
	}
}

template<typename T>
std::vector<std::list<unsigned int>* > const& CViterbiEncoderDecoder<T>::Encode(std::vector<T> const& p_rccData)
{
	m_cConvolutionEncoder.Convolove(p_rccData);
	return m_cConvolutionEncoder.GetSymbols();
}

template<typename T>
std::vector<T> const& CViterbiEncoderDecoder<T>::Decode(std::vector<std::list<unsigned int>* > const& p_rccSymbols)
{
	m_cTrellisDiagram.PerformStages(p_rccSymbols);
	m_cViterbiDecoder.Decode();
	return m_cViterbiDecoder.GetDecodedData();
}

template<typename T>
std::ostream& operator<<(std::ostream& p_rcOs, CViterbiEncoderDecoder<T> const& p_rccEncoderDecoder)
{
	p_rcOs << "cluster size: " << p_rccEncoderDecoder.m_cnSize << std::endl;
	p_rcOs << p_rccEncoderDecoder.m_cPolynomialGenerator;
	p_rcOs << *(p_rccEncoderDecoder.m_pcNextStateLookupTable);
	p_rcOs << *(p_rccEncoderDecoder.m_pcPrevStateLookupTable);
	p_rcOs << p_rccEncoderDecoder.m_cNextSymbolLookupTable;
	p_rcOs << p_rccEncoderDecoder.m_cConvolutionEncoder;
	p_rcOs << p_rccEncoderDecoder.m_cTrellisDiagram;
	p_rcOs << p_rccEncoderDecoder.m_cViterbiDecoder;

	return p_rcOs;
}

