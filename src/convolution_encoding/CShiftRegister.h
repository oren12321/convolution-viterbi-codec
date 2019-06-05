/**
 * @file CShiftRegister.h
 * @brief Class definition for shift register
 * @author Oren Cohen
 * @bug No known bugs
 */

#ifndef CONVOLUTION_ENCODING_CSHIFTREGISTER_H_
#define CONVOLUTION_ENCODING_CSHIFTREGISTER_H_

#include <list>

#include "../lookup_tables/CNextStateLookupTable.cuh"
#include "../lookup_tables/CNextEncodedSymbolLookupTable.h"

/**
 * Shift register for convolutional encoding
 *
 * Usage example:
 * @code
 * CShiftRegister shifter(stateLookupTable, symbolsLookupTable);
 * @endcode
 */
class CShiftRegister
{
public:
	/**
	 * Parameterized constructor for shift register
	 * @param p_pcNextStateLookupTable Lookup table for reviling next state of the register
	 * @param p_pcNextEncodedSymbolLookupTable Lookup table for taking the possible symbols according to state
	 */
	CShiftRegister(CNextStateLookupTable const& p_rccNextStateLookupTable,
			CNextEncodedSymbolLookupTable const& p_rccNextEncodedSymbolLookupTable);

	/**
	 * Destructor
	 */
	virtual ~CShiftRegister(void);

	/**
	 * Get encoded symbol
	 * @param p_nInput Given input (0 or 1) for symbol generation
	 * @return Encoded symbol
	 */
	std::list<unsigned int> const& ProcessSymbol(unsigned int const p_cnInput);

private:

	/**
	 * The current value of the shift register
	 */
	unsigned int m_nRegisters;

	/**
	 * The shift register lookup tables
	 */
	CNextStateLookupTable const& m_rccNextStateLookupTable;
	CNextEncodedSymbolLookupTable const& m_rccNextEncodedSymbolLookupTable;
};

#endif /* CONVOLUTION_ENCODING_CSHIFTREGISTER_H_ */
