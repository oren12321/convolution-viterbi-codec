

#include "CConvolutionEncoderKToN.h"

#include <stddef.h>
#include <iostream>
#include <string.h>

#include <list>
#include <vector>

#include "CShiftRegister.h"

template<typename T>
CConvolutionEncoderKToN<T>::CConvolutionEncoderKToN(CNextStateLookupTable const& p_rccNextStateLookupTable,
		CNextEncodedSymbolLookupTable const& p_rccNextEncodedSymbolLookupTable, unsigned int const p_cnK, unsigned int const p_cnN) :
		m_cShiftRegister(p_rccNextStateLookupTable, p_rccNextEncodedSymbolLookupTable), m_nMemory(p_cnK), m_nSymbolLength(p_cnN)
{

}

template<typename T>
CConvolutionEncoderKToN<T>::~CConvolutionEncoderKToN(void)
{
}

template<typename T>
void CConvolutionEncoderKToN<T>::Convolove(std::vector<T> const& p_pnData)
{
	// First check what is the number of the encoded symbols
	// In case that the needed length is different from the last time,
	// release the current buffer and allocate a new one
	// If the number of symbols is identical to the last encoding, only
	// clear the buffer
	// Because we are encode each bit of the data separately, the size of the symbols buffer
	// is equal to : data_size * bits_in_one_data_element + memory_size (needed for decoding)
	//
	// Example:
	// 		data = 0xA => order = 0101b => symbol for each bit
	unsigned int const l_cnNumberOfSymbols = p_pnData.size() * sizeof(T) * 8 + m_nMemory;
	if(m_cSymbols.size() != l_cnNumberOfSymbols)
	{
		m_cSymbols.resize(l_cnNumberOfSymbols);
	}

	// Get pointer to the start of the buffer
	std::list<unsigned int> const** l_pccSymbolsPtr = const_cast<std::list<unsigned int> const**>(&m_cSymbols[0]);

	// Fill the symbols buffer from the last data bit to the first one
	for (T const* l_pnDataPtr = &p_pnData[0] + (p_pnData.size() - 1); l_pnDataPtr >= &p_pnData[0]; l_pnDataPtr--)
	{
		// For each data element, encode from its LSB to its MSB
		unsigned int l_nBitsInTCount = sizeof(T) * 8;
		T l_nDataItem = *l_pnDataPtr;
		while (l_nBitsInTCount-- > 0)
		{

			// Activate the shift register for the current bit
			*l_pccSymbolsPtr = &m_cShiftRegister.ProcessSymbol(l_nDataItem & 1);
			l_pccSymbolsPtr++;
			l_nDataItem >>= 1;

		}
	}

	// After filling the data, add the spares (The decoding will performed in the
	// opposite way)
	for (unsigned int l_nSIndex = 0; l_nSIndex < m_nMemory; l_nSIndex++)
	{
		*l_pccSymbolsPtr = &m_cShiftRegister.ProcessSymbol(0);
		l_pccSymbolsPtr++;
	}

}

template<typename T>
std::vector<std::list<unsigned int>* > const& CConvolutionEncoderKToN<T>::GetSymbols(void) const
{
	return m_cSymbols;
}

template<typename T>
unsigned int CConvolutionEncoderKToN<T>::GetNumberOfSymbols(void) const
{
	return m_cSymbols.size();
}

template<typename T>
std::ostream& operator<<(std::ostream& p_rcOs, CConvolutionEncoderKToN<T> const& p_rccEncoder)
{
	p_rcOs << "symbols (count = " << p_rccEncoder.m_cSymbols.size() << "):" << std::endl << "\t";
	for(unsigned int l_nI = 0; l_nI < p_rccEncoder.m_cSymbols.size(); l_nI++)
	{
		p_rcOs << CPolynomialGenerator::GetBinaryString(*(p_rccEncoder.m_cSymbols[l_nI]), p_rccEncoder.m_nSymbolLength) << " ";
	}

	p_rcOs << std::endl;

    return p_rcOs;
}
