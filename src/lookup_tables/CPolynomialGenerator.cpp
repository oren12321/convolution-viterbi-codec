/**
 * @file CPolynomialGenerator.cpp
 * @brief Class implementation for polynomials generator
 * @author Oren Cohen
 * @bug No known bugs
 */

#include "CPolynomialGenerator.h"

#include <stddef.h>

#include <string>
#include <iostream>
#include <bitset>

#include <list>

CPolynomialGenerator::CPolynomialGenerator(unsigned int const p_cnN, unsigned int const p_cnK) :
		m_cnDegree(p_cnK - 1), m_cPolynomials(p_cnN)
{
	// Create the polynomials
	CreatePolynomials();
}

CPolynomialGenerator::~CPolynomialGenerator(void)
{
}

unsigned int CPolynomialGenerator::GetDegree(void) const
{
	return m_cnDegree;
}

unsigned int CPolynomialGenerator::GetSize(void) const
{
	return m_cPolynomials.size();
}

unsigned int CPolynomialGenerator::GetPolynomial(unsigned int const p_cnIndex) const
{
	return m_cPolynomials[p_cnIndex];
}

void CPolynomialGenerator::CreatePolynomials(void)
{
	// Calculate the basic polynomial (10....01)
	// It needs to always be in this pattern to ensure that all its
	// bits will participate in the convolutional encoding
	unsigned int const l_cnBoundries = 1 + (1 << m_cnDegree);

	// Now for the middle part just put the numbers from 0 to m_nSize
	unsigned int l_nBase = 0;
	for(std::vector<unsigned int>::iterator l_cIt = m_cPolynomials.begin(); l_cIt != m_cPolynomials.end(); ++l_cIt)
	{
		*l_cIt = l_cnBoundries + /* counter alignment */((l_nBase++) << 1);
	}
}


std::string const CPolynomialGenerator::GetBinaryString(unsigned int const p_cnNumber, unsigned int const p_cnLength)
{
	// Convert the number to binary string with maximum size and the cut it
	// according to given length
	std::string l_c64BitEncoding(std::bitset<sizeof(p_cnNumber) * 8>(p_cnNumber).to_string());
	return l_c64BitEncoding.substr(l_c64BitEncoding.size() - p_cnLength);
}

std::string const CPolynomialGenerator::GetBinaryString(std::list<unsigned int> const& p_rccNumber, unsigned int const p_cnLength)
{
	std::string l_cEncoding;
	int l_nLength = p_cnLength;
	int const l_cnMaxClusterLength = sizeof(unsigned int) * 8;
	// Iterate on the collection
	for(std::list<unsigned int>::const_iterator l_cIt = p_rccNumber.begin(); l_cIt != p_rccNumber.end(); ++l_cIt)
	{
		int l_nReducedLength;
		if(l_nLength - l_cnMaxClusterLength > 0)
		{
			l_cEncoding += GetBinaryString(*l_cIt, l_cnMaxClusterLength);
			l_nReducedLength = l_cnMaxClusterLength;
		}
		else
		{
			l_cEncoding += GetBinaryString(*l_cIt, l_nLength);
			l_nReducedLength = l_nLength;
		}
		l_nLength -= l_nReducedLength;
	}
	return l_cEncoding;
}

std::ostream& operator<<(std::ostream& p_rcOs, CPolynomialGenerator const& p_rccGenerator)
{
	p_rcOs << "polynomial generator of degree " << p_rccGenerator.m_cnDegree << ": " << std::endl << "\t";

	for(std::vector<unsigned int>::const_iterator l_cIt = p_rccGenerator.m_cPolynomials.begin(); l_cIt != p_rccGenerator.m_cPolynomials.end(); ++l_cIt)
	{
		p_rcOs << CPolynomialGenerator::GetBinaryString(*l_cIt, p_rccGenerator.m_cnDegree + 1) << " ";
	}

	p_rcOs << std::endl;

    return p_rcOs;
}
