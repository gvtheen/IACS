#ifndef CModelClusterSUPPORT_H
#define CModelClusterSUPPORT_H

#include <Eigen/Dense>
#include "CModelBase.h"
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

class CModelSupport;
class CModelMoleculeAdsorbent;
class CParameter;

class CModelClusterSupport:public CModelBase
{
    public:
        typedef enum EA_001286{
                     POLYHEDRON = 0x98126,
                         SPHERE = 0x98127}CLUSTER_MODEL_TYPE;

        CModelClusterSupport(CParameter*,CATAZJUT::CPeriodicFramework**,size_t);
        virtual ~CModelClusterSupport();

        CModelBase* clone();

        //virtual function from CModelBase
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
          CModelSupport*      m_pSupport;
    CModelMoleculeAdsorbent*      m_pAdsorbMolecule;

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
#endif // CModelClusterSupport_H
