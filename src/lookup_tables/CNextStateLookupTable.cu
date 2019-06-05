
#include "CNextStateLookupTable.cuh"

#include "CStateLookupTable.cuh"

CNextStateLookupTable::CNextStateLookupTable(unsigned int const p_cnK) :
		CStateLookupTable(p_cnK
#ifdef _USE_CUDA_
				, LookupTableType_Next // Set table type
#endif
				)
{

}

CNextStateLookupTable::~CNextStateLookupTable(void)
{
}

unsigned int CNextStateLookupTable::Shift(unsigned int const p_cnState, unsigned int const p_cnInput) const
{
	// Shift right once and put the input in the MSB
	return (p_cnState >> 1) + ((p_cnInput == 0) ? (0) : (1 << (m_cnMemory - 1)));
}

CNextStateLookupTable* CNextStateLookupTable::Create(unsigned int const p_cnK)
{
	CNextStateLookupTable* l_pcInstance = new CNextStateLookupTable(p_cnK);
	l_pcInstance->InitializeStateLookupTable();
	return l_pcInstance;
}
