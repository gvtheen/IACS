#include "CCalcModeStruct.h"
#include "CPeriodicFramework.h"

namespace CALCZJUT{

CCalcModeStruct::CCalcModeStruct(CParameter* temParameter)
:m_pParameter(temParameter)
{

}

CCalcModeStruct::~CCalcModeStruct()
{
    //dtor
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
void CCalcModeStruct::removeStructureOfGene()
{
   for(size_t i=0;i<m_PopuPeriodicFramework.size();i++)
       delete m_PopuPeriodicFramework[i];
   m_PopuPeriodicFramework.clear();
}


}
