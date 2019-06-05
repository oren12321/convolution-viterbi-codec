/**
 * @file CNextEncodedSymbolLookupTable.h
 * @brief Class definition for next predicted symbol lookup table
 * @author Oren Cohen
 * @bug No known bugs
 */

#ifndef LOOKUP_TABLES_CNEXTENCODEDSYMBOLLOOKUPTABLE_H_
#define LOOKUP_TABLES_CNEXTENCODEDSYMBOLLOOKUPTABLE_H_

#include <stddef.h>
#include <iostream>

#include <list>

#include "CPolynomialGenerator.h"

/**
 * Lookup table of possible encoded symbols (by convolution) \n
 * according to state
 *
 * Usage example:
 * @code
 * CNextEncodedSymbolLookupTable table(CPolynomialGenerator(3,2));
 * This will generates :
 * state	input:0	input:1
 * 000	00	11
 * 001	11	00
 * 010	01	10
 * 011	10	01
 * 100	11	00
 * 101	00	11
 * 110	10	01
 * 111	01	10
 * @endcode
 */
class CNextEncodedSymbolLookupTable
{
public:

	/**
	 * Parameterized constructor that build the lookup table according to \n
	 * polynomials generator
	 * @param p_pcPolynomialGenerator The given polynomials generator
	 */
	CNextEncodedSymbolLookupTable(CPolynomialGenerator const& p_rccPolynomialGenerator);
	/**
	 * Destructor
	 */
	virtual ~CNextEncodedSymbolLookupTable(void);

	/**
	 * Get next encoded symbol according to state and input
	 * @param p_nCurrentState The given current state
	 * @param p_nInput The given input (0 or 1)
	 * @return The next encoded symbol
	 */
	std::list<unsigned int> const& GetNextEncodedSymbol(unsigned int const p_cnCurrentState, unsigned int const p_cnInput) const;
	/**
	 * Get next possible encoded symbols according to the current state
	 * @param p_nCurrentState The current state
	 * @return Possible symbols
	 */
	std::list<unsigned int> const* GetNextEncodedSymbols(unsigned int p_cnCurrentState) const;

	/**
	 * Get table size
	 * @return The table size (number of states)
	 */
	unsigned int GetSize(void) const;

	/**
	 * Get number of symbol clusters
	 * @return The number of clusters
	 */
	unsigned int GetNumberOfClusters(void) const;

	/**
	 * Get the table as STL vector
	 * @return The table as vector
	 */
	std::vector<std::list<unsigned int> > const& AsVector(void) const;

	/**
	 * For Debug purposes - overload the << operator for this class
	 * @param p_rcOs The output stream
	 * @param p_rcTable The symbols lookup table for output
	 * @return The current output stream
	 */
	friend std::ostream& operator<<(std::ostream& p_rcOs, CNextEncodedSymbolLookupTable const& p_rccTable);

private:

	/**
	 * Calculate one bit of the next symbol according to given parameters
	 * @param p_nGIndex The index of the appropriate polynomial
	 * @param p_nCurrentState The current state
	 * @param p_nInput The current input for the state (0 or 1)
	 * @param p_nMemory The memory size of the symbol
	 */
	unsigned int CalculateSymbolBit(unsigned int const p_cnGIndex, unsigned int const p_cnCurrentState, unsigned int const p_cnInput, unsigned int const p_cnMemory);
	/**
	 * Initialize the symbols lookup table
	 */
	void InitializeNextEncodedSymbolLookupTable(void);

	/**
	 * Given polynomials generator
	 */
	CPolynomialGenerator const& m_rccPolynomialGenerator;

	/**
	 * The data of the lookup table (a matrix)
	 */
	std::vector<std::list<unsigned int> > m_cTable;
	/*
	 * Table entry size (number of columns), it always equals to 2
	 */
	unsigned int const m_cnEntrySize;
};

#endif /* LOOKUP_TABLES_CNEXTENCODEDSYMBOLLOOKUPTABLE_H_ */
