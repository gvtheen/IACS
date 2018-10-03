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
void CCalcClusterSupport::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{

}
std::vector<double>* CCalcClusterSupport::getGeneValuefromStruct()const
{

}
std::vector<GeneVAR>* CCalcClusterSupport::GeneVARRange()
{

}

}
