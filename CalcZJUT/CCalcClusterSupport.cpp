#include "CCalcClusterSupport.h"
#include "CCalcSupportBase.h"
#include "CCalcMoleculeAdsorbent.h"
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CCrystalPlanes.h"
#include "../CataZJUT/CCrystalPlanes.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../Util/Bitset.h"
#include "../Util/foreach.h"
#include "../Util/Point-Vector.h"

using util::Bitset;
using util::Log;
using util::Matrix;

namespace CALCZJUT{

CCalcClusterSupport::CCalcClusterSupport(CParameter* mPara,CATAZJUT::CPeriodicFramework** copy_ppPeriodicFramework)
:CCalcModeStruct(mPara)
{
    m_pPeriodicFramework = new CATAZJUT::CPeriodicFramework(mPara);
    m_ppBackupPeriodicFramework = copy_ppPeriodicFramework;
    m_pSupport = nullptr;
    m_pAdsorbMolecule = nullptr;
    m_pCrystalPlanes = nullptr;
    m_IsCrystalPlanePerceived = false;
}
CCalcModeStruct* CCalcClusterSupport::clone()
{
   CCalcClusterSupport* res = new CCalcClusterSupport(this->m_pParameter,this->m_ppBackupPeriodicFramework);
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
   m_pSupport= new CCalcSupportBase(this->m_pPeriodicFramework,m_BitbackupSupport);
}
void CCalcClusterSupport::createMoleAdsorb(const Bitset& mht)
{
   assert(m_pPeriodicFramework);
   m_BitbackupAdsorbMolecule = mht;
   m_pAdsorbMolecule = new CCalcMoleculeAdsorbent(this->m_pPeriodicFramework,m_BitbackupAdsorbMolecule);
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
   if(currentGeneVARible.size()!=0)
       currentGeneVARible.clear();

   if( !this->m_IsCrystalPlanePerceived )
        perceiveCrystalPlane();
    // index of crystal planes in the cluster.                    1st gene descriptor.
    // max is m_pCrystalPlanes->crystalPlaneNum() - 1; substracting 1 is necessary!
    currentGeneVARible.push_back({0,m_pCrystalPlanes->crystalPlaneNum()-1,0.5});

    // Height the adsorbing molecule on the index crystal plane.  2nd gene descriptor.

    // ratio

    // angle


}
void CCalcClusterSupport::perceiveCrystalPlane()
{
    Eigen::MatrixXd *atom_Coordinate =  new (Eigen::MatrixXd)(m_pSupport->atomCount(),3);
    size_t index=0;
    Point3 tempP;
    foreach(CATAZJUT::CAtom* atom_s,this->m_pSupport->atoms()){
        tempP = atom_s->position();
        for(size_t i=0;i<3;i++){
            (*atom_Coordinate)(index,i) =tempP[i];
        }
        index++;
    }
    this->m_pCrystalPlanes =  new CATAZJUT::CCrystalPlanes(atom_Coordinate);
    this->m_pCrystalPlanes->CreateCrystalPlane();
    this->m_IsCrystalPlanePerceived = true;
}


}
