#ifndef CCALCCLUSTERSUPPORT_H
#define CCALCCLUSTERSUPPORT_H

#include <Eigen/Dense>
#include "CCalcModeStruct.h"
#include "../GaZJUT/GaDeclaration.h"
#include "../Util/Bitset.h"

using GAZJUT::GeneVAR;
using util::Bitset;

// import crystal plane object from namespace CATAZJUT.
namespace CATAZJUT{
   class CCrystalPlanes;
}

namespace CALCZJUT{

class CCalcSupportBase;
class CCalcMoleculeAdsorbent;
class CParameter;

class CCalcClusterSupport:public CCalcModeStruct
{
    public:
        CCalcClusterSupport(CParameter*,CATAZJUT::CPeriodicFramework**);
        virtual ~CCalcClusterSupport();

        CCalcModeStruct* clone();

        //virtual function from CCalcModeStruct
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        void getGeneValuefromStruct(std::vector<double>&);
        void GeneVARRange(std::vector<GeneVAR>&);

        void createSupport(const Bitset &);
        void createMoleAdsorb(const Bitset &);

        void perceiveCrystalPlane();

        Bitset SupportBit();
        Bitset MoleAdsorbBit();

        CATAZJUT::CCrystalPlanes* crystalPlanes();
        void setCrystalPlanes(CATAZJUT::CCrystalPlanes*);

    protected:

    private:
          CCalcSupportBase*      m_pSupport;
    CCalcMoleculeAdsorbent*      m_pAdsorbMolecule;
  CATAZJUT::CCrystalPlanes*      m_pCrystalPlanes;
                     Bitset      m_BitbackupSupport;
                     Bitset      m_BitbackupAdsorbMolecule;
                       bool      m_IsCrystalPlanePerceived;



};



}
#endif // CCalcClusterSupport_H
