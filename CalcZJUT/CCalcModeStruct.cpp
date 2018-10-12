#include "CCalcModeStruct.h"
#include "CPeriodicFramework.h"

namespace CALCZJUT{

CCalcModeStruct::CCalcModeStruct(CParameter* temParameter)
:m_pParameter(temParameter)
{
    m_IsNeedRandomInit=true;
}

CCalcModeStruct::~CCalcModeStruct()
{
    //dtor
}
void CCalcModeStruct::init()
{

}
CATAZJUT::CPeriodicFramework* CCalcModeStruct::periodicFramework()
{
    return m_pPeriodicFramework;
}
CATAZJUT::CPeriodicFramework* CCalcModeStruct::periodicFramework(size_t index_int)
{
    if( index_int > this->m_PopuPeriodicFramework.size()-1 && index_int<0 )
        return nullptr;

    return m_PopuPeriodicFramework[index_int];
}
void CCalcModeStruct::setPeriodicFramekwork(CATAZJUT::CPeriodicFramework* mbf)
{
    m_pPeriodicFramework = mbf;
}
void CCalcModeStruct::createSupport(Bitset & mt)
{
    // no doing
}
void CCalcModeStruct::createMoleAdsorb(Bitset & mt)
{
    // no doing
}
void CCalcModeStruct::createStructureAtGene()
{

}
void CCalcModeStruct::setRandomInitState(const bool& mht)
{
    this->m_IsNeedRandomInit=mht;
}
bool CCalcModeStruct::RandomInitState()
{
    return this->m_IsNeedRandomInit;
}
void CCalcModeStruct::removeStructureOfGene()
{
   for(size_t i=0;i<m_PopuPeriodicFramework.size();i++)
       delete m_PopuPeriodicFramework[i];
   m_PopuPeriodicFramework.clear();
}
void CCalcModeStruct::init()
{
    // nothing to do
}
std::vector<std::pair<std::string,size_t>>& CCalcModeStruct::chemicalFormula()
{

}
void CCalcModeStruct::setChemicalFormula(const std::vector<std::pair<std::string,size_t>>& mth)
{

}



}
