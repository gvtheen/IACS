#include "CCalcClusterSupport.h"
namespace CALCZJUT{

CCalcClusterSupport::CCalcClusterSupport(CParameter* mPara)
:CCalcModeStruct(mPara)
{
    m_pPeriodicFramework = new CATAZJUT::CPeriodicFramework(mPara);
    m_backupPeriodicFramework = m_pPeriodicFramework;
    m_pSupport = nullptr;
    m_pAdsorbMolecule = nullptr;
    m_pCrystalPlanes = nullptr;
}
CCalcModeStruct* CCalcClusterSupport::clone()
{
   CCalcClusterSupport* res = new CCalcClusterSupport(this->m_pParameter);
   res->m_pGeneVAR->assign(this->m_pGeneVAR->begin(),this->m_pGeneVAR->end());
   res->createMoleAdsorb(this->MoleAdsorbBit());
   res->createSupport(this->SupportBit());
   res->setRandomInitState(this->RandomInitState());
   res->setCrystalPlanes(this->crystalPlanes());
   return res;

}
CCalcClusterSupport::~CCalcClusterSupport()
{
    //dtor
}

void CCalcClusterSupport::createSupport(const Bitset& mht)
{
   assert(m_pPeriodicFramework);
   m_BitbackupSupport = mht;
   m_pSupport= new CCalcSupportBase(this->m_pPeriodicFramework,mht);
}
void CCalcClusterSupport::createMoleAdsorb(const Bitset& mht)
{
   assert(m_pPeriodicFramework);
   m_BitbackupAdsorbMolecule = mht;
   m_pAdsorbMolecule = new CCalcMoleculeAdsorbent(this->m_pPeriodicFramework,mht);
}
void CCalcClusterSupport::createStructureAtGene()
{
   CCalcModeStruct::createStructureAtGene();
   m_pSupport->setConfiguration(this->periodicFramework());
   m_pAdsorbMolecule->setConfiguration(this->periodicFramework());
}
// crystal plane is set and get.
CATAZJUT::CCrystalPlanes* CCalcClusterSupport::crystalPlanes()
{
   return this->m_pCrystalPlanes;
}
void CCalcClusterSupport::setCrystalPlanes(CATAZJUT::CCrystalPlanes* mth)
{
   if(this->m_pCrystalPlanes!=nullptr)
       delete this->m_pCrystalPlanes;
   this->m_pCrystalPlanes = mth->Clone();
}
void CCalcClusterSupport::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{

}
void CCalcClusterSupport::getGeneValuefromStruct(std::vector<double>& currentGeneRealValue)
{

}
void CCalcClusterSupport::GeneVARRange(std::vector<GeneVAR>& currentGeneVARible)
{

}

}
