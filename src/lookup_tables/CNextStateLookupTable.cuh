/**symbols lookup table
 * @file CNextStateLookupTable.h
 * @brief Class definition for next states lookup table
 * @author Oren Cohen
 * @bug No known bugs
 */

#ifndef LOOKUP_TABLES_CNEXTSTATELOOKUPTABLE_H_
#define LOOKUP_TABLES_CNEXTSTATELOOKUPTABLE_H_

#include "CStateLookupTable.cuh"

/**
 * Lookup table of possible next states
 *
 * Usage example:
 * @code
 * CNextStateLookupTable table(3);
 * This will generates :
 * state	input:0	input:1
 * 000	000	100
 * 001	000	100
 * 010	001	101
 * 011	001	101
 * 100	010	110
 * 101	010	110
 * 110	011	111
 * 111	011	111
 * @endcode
 */
class CNextStateLookupTable: public CStateLookupTable
{
public:
	static CNextStateLookupTable* Create(unsigned int const p_cnK);

	virtual ~CNextStateLookupTable(void);

private:
	CNextStateLookupTable(unsigned int const p_cnK);

	unsigned int Shift(unsigned int const p_cnState, unsigned int const p_cnInput) const;
};

#endif /* LOOKUP_TABLES_CNEXTSTATELOOKUPTABLE_H_ */
