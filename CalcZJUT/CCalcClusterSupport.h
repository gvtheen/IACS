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
   class CSphere;
}

namespace CALCZJUT{

class CCalcSupportBase;
class CCalcMoleculeAdsorbent;
class CParameter;

class CCalcClusterSupport:public CCalcModeStruct
{
    public:
        typedef enum EA_001286{
                     POLYHEDRON = 0x98126,
                         SPHERE = 0x98127}CLUSTER_MODEL_TYPE;

        CCalcClusterSupport(CParameter*,CATAZJUT::CPeriodicFramework**);
        virtual ~CCalcClusterSupport();

        CCalcModeStruct* clone();

        //virtual function from CCalcModeStruct
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        void getGeneValuefromStruct(std::vector<double>&);
        void GeneVARRange(std::vector<GeneVAR>&);

        void createSupport(const Bitset &);
        void createMoleAdsorb(const Bitset &);

        void perceiveClusterModel();

        Bitset SupportBit();
        Bitset MoleAdsorbBit();

        CATAZJUT::CCrystalPlanes* crystalPlanes();
        void setCrystalPlanes(CATAZJUT::CCrystalPlanes*);

    protected:
        void eliminateCloseContacts(double distanceCutOff=1.0);
    private:
          CCalcSupportBase*      m_pSupport;
    CCalcMoleculeAdsorbent*      m_pAdsorbMolecule;

  CATAZJUT::CCrystalPlanes*      m_pCrystalPlanes;
         CATAZJUT::CSphere*      m_pSphere;
                     // define the index of support and adsorbmolecule in whole configuration.
                     Bitset      m_BitbackupSupport;
                     Bitset      m_BitbackupAdsorbMolecule;
                     // define the model of current cluster
         CLUSTER_MODEL_TYPE      m_ClusterModelType;
                       bool      m_ClusterModelPerceived;



};



}
#endif // CCalcClusterSupport_H
