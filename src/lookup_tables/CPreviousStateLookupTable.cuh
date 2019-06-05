/**
 * @file CPreviousStateLookupTable.h
 * @brief Class definition for previous states lookup table
 * @author Oren Cohen
 * @bug No known bugs
 */

#ifndef LOOKUP_TABLES_CPREVIOUSSTATELOOKUPTABLE_H_
#define LOOKUP_TABLES_CPREVIOUSSTATELOOKUPTABLE_H_

#include "CStateLookupTable.cuh"

/**
 * Lookup table of possible previous states
 *
 * Usage example:
 * @code
 * CPreviousStateLookupTable table(3);
 * This will generates :
 * state	input:0	input:1
 * 000	000	001
 * 001	010	011
 * 010	100	101
 * 011	110	111
 * 100	000	001
 * 101	010	011
 * 110	100	101
 * 111	110	111
 * @endcode
 */
class CPreviousStateLookupTable: public CStateLookupTable
{
public:
	static CPreviousStateLookupTable* Create(unsigned int const p_cnK);

	virtual ~CPreviousStateLookupTable(void);

private:
	CPreviousStateLookupTable(unsigned int const p_cnK);

	unsigned int Shift(unsigned int const p_cnState, unsigned int const p_cnInput) const;

	unsigned int const m_cnMaxState;
};

#endif /* LOOKUP_TABLES_CPREVIOUSSTATELOOKUPTABLE_H_ */
