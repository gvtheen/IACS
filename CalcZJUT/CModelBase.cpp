#include "CModelBase.h"
#include "../CataZJUT/CPeriodicFramework.h"

namespace CALCZJUT{

CModelBase::CModelBase(CParameter* temParameter)
:m_pParameter(temParameter)
{
    m_IsNeedRandomInit=true;
}

CModelBase::~CModelBase()
{
    //dtor
}
void CModelBase::init()
{

}
CATAZJUT::CPeriodicFramework* CModelBase::periodicFramework()
{
    return m_pPeriodicFramework;
}
void CModelBase::setPeriodicFramekwork(CATAZJUT::CPeriodicFramework* mbf)
{
    m_pPeriodicFramework = mbf;
}
void CModelBase::createSupport(Bitset & mt)
{
    // no doing
}
void CModelBase::createMoleAdsorb(Bitset & mt)
{
    // no doing
}
void CModelBase::createStructureAtGene()
{

}
void CModelBase::setRandomInitState(const bool& mht)
{
    this->m_IsNeedRandomInit=mht;
}
bool CModelBase::RandomInitState()
{
    return this->m_IsNeedRandomInit;
}
void CModelBase::removeStructureOfGene()
{
   for(size_t i=0;i<m_PopuPeriodicFramework.size();i++)
       delete m_PopuPeriodicFramework[i];
   m_PopuPeriodicFramework.clear();
}
void CModelBase::init()
{
    // nothing to do
}
std::vector<std::pair<std::string,size_t>>& CModelBase::chemicalFormula()
{

}
void CModelBase::setChemicalFormula(const std::vector<std::pair<std::string,size_t>>& mth)
{

}



}
