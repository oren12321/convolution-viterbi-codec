
#include "CNextEncodedSymbolLookupTable.h"

#include <stddef.h>
#include <iostream>

#include <iomanip>
#include <algorithm>
#include <cmath>

#include <list>

#include "CStateLookupTable.cuh"
#include "CPolynomialGenerator.h"

CNextEncodedSymbolLookupTable::CNextEncodedSymbolLookupTable(CPolynomialGenerator const& p_rccPolynomialGenerator) :
		m_rccPolynomialGenerator(p_rccPolynomialGenerator), m_cnEntrySize(2)
{
	// Build the lookup table
	InitializeNextEncodedSymbolLookupTable();
}

CNextEncodedSymbolLookupTable::~CNextEncodedSymbolLookupTable(void)
{
}

std::list<unsigned int> const& CNextEncodedSymbolLookupTable::GetNextEncodedSymbol(unsigned int const p_cnCurrentState, unsigned int const p_cnInput) const
{
	return GetNextEncodedSymbols(p_cnCurrentState)[p_cnInput];
}

std::list<unsigned int> const* CNextEncodedSymbolLookupTable::GetNextEncodedSymbols(unsigned int const p_cnCurrentState) const
{
	// Possible symbols according to current state
	return &m_cTable[p_cnCurrentState * 2];
}

unsigned int CNextEncodedSymbolLookupTable::CalculateSymbolBit(unsigned int const p_cnGIndex, unsigned int const p_cnCurrentState,
		unsigned int const p_cnInput, unsigned int const p_cnMemory)
{
	unsigned int l_nSymbolBit = 0; // Holds one bit of the generated symbol
	bool l_bFirstXor = true;

	// Get the appropriate polynomial
	// Each bit in the symbol calculated by other polynomial function
	unsigned int l_nGenerator = m_rccPolynomialGenerator.GetPolynomial(p_cnGIndex);

	// Get next possible symbol by given input (push left 0 or 1)
	unsigned int l_nRegisters = (p_cnCurrentState >> 1) + (p_cnInput << (p_cnMemory - 1));

	// Move over all the bits of the polynomial
	// If the current bit of the polynomial is set then perform accumulated XOR
	// with between the current symbol bit and the parallel bit to the polynomial
	// in the shift register
	// For example :
	// 		polynomial = 1001101
	//                   |  || |	=>	symbol_bit = 1^1^0^1
	//		register   = 1001011
	for (unsigned int l_nBIndex = 0; l_nBIndex < p_cnMemory; l_nBIndex++)
	{
		if ((l_nGenerator & 1) == 1)
		{
			if (l_bFirstXor) // For the first XOR operation accumulation doesn't needed
			{
				l_nSymbolBit = (l_nRegisters & 1);
				l_bFirstXor = false;
			}
			else
			{
				l_nSymbolBit ^= (l_nRegisters & 1);
			}
		}
		// Go to the next bits of the polynomial and the shift register
		l_nGenerator >>= 1;
		l_nRegisters >>= 1;
	}

	return l_nSymbolBit;
}

void CNextEncodedSymbolLookupTable::InitializeNextEncodedSymbolLookupTable(void)
{
	// Get number of XORs needed for symbol (its the polynomial generator size)
	unsigned int const l_cnNumberOfXors = m_rccPolynomialGenerator.GetSize();
	// Get the memory size, means bits size of the states
	unsigned int const l_cnMemory = m_rccPolynomialGenerator.GetDegree() + 1;

	// Allocate memort for the lookup table and its rows size (= max state size)
	m_cTable.resize((CStateLookupTable::GetMaxState(l_cnMemory) + 1) * 2);
	std::list<unsigned int>* l_ppnTablePtr = &m_cTable[0];

	const unsigned int l_cnMaxBitsInCluster = sizeof(unsigned int) * 8;

	// Loop over the rows
	for (unsigned int l_nCurrentState = 0; l_nCurrentState < m_cTable.size() / 2; l_nCurrentState++)
	{
		// Reveal the next possible symbol for every possible input (0 and 1)
		for (unsigned int l_nInput = 0; l_nInput < m_cnEntrySize; l_nInput++)
		{
			unsigned int l_nSymbol = 0;

			// Perform needed number of XORs
			for (unsigned int l_nBIndex = 0; l_nBIndex < l_cnNumberOfXors; l_nBIndex++)
			{
				// Before calculating the next bit, shift the symbol one time left
				l_nSymbol <<= 1;
				l_nSymbol += CalculateSymbolBit(l_nBIndex, l_nCurrentState, l_nInput, l_cnMemory);

				// When reaching bitsize(unsigned int) or symbol size, store the bits in the list
				if((l_nBIndex + 1) % l_cnMaxBitsInCluster == 0 || l_nBIndex == l_cnNumberOfXors - 1)
				{
					l_ppnTablePtr[l_nInput].push_back(l_nSymbol);
				}
			}

		}

		// Insert row record to lookup table
		l_ppnTablePtr += 2;

	}
}

unsigned int CNextEncodedSymbolLookupTable::GetSize(void) const
{
	return m_cTable.size();
}

unsigned int CNextEncodedSymbolLookupTable::GetNumberOfClusters(void) const
{
	return std::ceil(m_rccPolynomialGenerator.GetSize() / 32.0);
}

std::vector<std::list<unsigned int> > const& CNextEncodedSymbolLookupTable::AsVector(void) const
{
	return m_cTable;
}

std::ostream& operator<<(std::ostream& p_rcOs, CNextEncodedSymbolLookupTable const& p_rccTable)
{

	p_rcOs << "next symbol lookup table:" << std::endl;

	unsigned int l_nMemory = p_rccTable.m_rccPolynomialGenerator.GetDegree() + 1;
	unsigned int l_nNumberOfXors = p_rccTable.m_rccPolynomialGenerator.GetSize();

	unsigned int l_nStateWidth = l_nMemory;

	unsigned int l_nStateColumnWidth = std::max(5, (int)l_nStateWidth);
	unsigned int l_nInputColumnWidth = std::max(8, (int)l_nNumberOfXors);

	unsigned int l_nDelimeterWidth = std::max(23, (int)(l_nInputColumnWidth * 2 + l_nStateColumnWidth + 2));

	p_rcOs << '+';
	for(unsigned int l_nI = 0; l_nI < l_nDelimeterWidth; l_nI++)
	{
		p_rcOs << '-';
	}
	p_rcOs << '+' << std::endl;

	p_rcOs << '|' << std::setw(l_nStateColumnWidth) << std::left << "state" << '|'
			<< std::setw(l_nInputColumnWidth) << std::left << "input: 0" << '|'
			<< std::setw(l_nInputColumnWidth) << std::left << "input: 1" << '|'
			<< std::endl;

	p_rcOs << '|';
	for(unsigned int l_nI = 0; l_nI < l_nDelimeterWidth; l_nI++)
	{
		p_rcOs << '-';
	}
	p_rcOs << '|' << std::endl;


	for(unsigned int l_nRI = 0; l_nRI < p_rccTable.m_cTable.size(); l_nRI+=2)
	{
		p_rcOs << '|' << std::setw(l_nStateColumnWidth) << std::left << CPolynomialGenerator::GetBinaryString(l_nRI / 2, l_nMemory) << '|'
				<< std::setw(l_nInputColumnWidth) << std::left << CPolynomialGenerator::GetBinaryString(p_rccTable.m_cTable[l_nRI + 0], l_nNumberOfXors) << '|'
				<< std::setw(l_nInputColumnWidth) << std::left << CPolynomialGenerator::GetBinaryString(p_rccTable.m_cTable[l_nRI + 1], l_nNumberOfXors) << '|'
				<< std::endl;
	}

	p_rcOs << '+';
	for(unsigned int l_nI = 0; l_nI < l_nDelimeterWidth; l_nI++)
	{
		p_rcOs << '-';
	}
	p_rcOs << '+' << std::endl;

    return p_rcOs;
}
