
#include "CShiftRegister.h"

#include <list>

#include "../lookup_tables/CNextStateLookupTable.cuh"
#include "../lookup_tables/CNextEncodedSymbolLookupTable.h"

CShiftRegister::CShiftRegister(CNextStateLookupTable const& p_rccNextStateLookupTable,
		CNextEncodedSymbolLookupTable const& p_rccNextEncodedSymbolLookupTable) :
		m_nRegisters(0), m_rccNextStateLookupTable(p_rccNextStateLookupTable), m_rccNextEncodedSymbolLookupTable(
				p_rccNextEncodedSymbolLookupTable)
{
}

CShiftRegister::~CShiftRegister(void)
{
}

std::list<unsigned int> const& CShiftRegister::ProcessSymbol(unsigned int const p_cnInput)
{
	// Take symbol according to current state of the shift register (its state)
	std::list<unsigned int> const& l_nSymbol = m_rccNextEncodedSymbolLookupTable.GetNextEncodedSymbol(m_nRegisters, p_cnInput);

	// Update the shift register value to the next state
	m_nRegisters = m_rccNextStateLookupTable.GetSequentialState(m_nRegisters, p_cnInput);

	// Return the symbol
	return l_nSymbol;
}
