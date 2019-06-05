
#include "CPreviousStateLookupTable.cuh"

#include "CStateLookupTable.cuh"

CPreviousStateLookupTable::CPreviousStateLookupTable(unsigned int const p_cnK) :
		CStateLookupTable(p_cnK
#ifdef _USE_CUDA_
				, LookupTableType_Prev // Set table type
#endif
				), m_cnMaxState(GetMaxState(m_cnMemory))
{
}

CPreviousStateLookupTable::~CPreviousStateLookupTable(void)
{
}

unsigned int CPreviousStateLookupTable::Shift(unsigned int const p_cnState, unsigned int const p_cnInput) const
{
	// Shift left once and put the input in the LSB
	return ((p_cnState << 1) & m_cnMaxState) + p_cnInput;
}

CPreviousStateLookupTable* CPreviousStateLookupTable::Create(unsigned int const p_cnK)
{
	CPreviousStateLookupTable* l_pcInstance = new CPreviousStateLookupTable(p_cnK);
	l_pcInstance->InitializeStateLookupTable();
	return l_pcInstance;
}

