#ifndef CMODELCLUSTERLOADED2DSUPPORT_H
#define CMODELCLUSTERLOADED2DSUPPORT_H

#include <vector>
#include "../Util/Bitset.h"
#include "../GaZJUT/GaDeclaration.h"
#include "CModelBase.h"
#include "../Util/Point-Vector.h"
#include "../IACS.h"

using util::Point3;
using util::Vector3;

namespace CATAZJUT{
  class CConfigurationBase;
  class CPlane;
  class CSphere;

}
namespace CALCZJUT{

class CModelMoleculeAdsorbent;
class CModelSupport;
class CParameter;

class CModelClusterLoaded2DSupport:public CModelBase
{
    public:
        typedef enum LA_001{A_AXIS=0x40AB,
                            B_AXIS=0x40AC,
                            C_AXIS=0x40AD,
                          NONE_DIR=0x40AE}LATT_DIRECTION;

        CModelClusterLoaded2DSupport(CParameter*,CATAZJUT::CConfigurationBase**,size_t);
        virtual ~CModelClusterLoaded2DSupport();

        CModelBase* clone();   //clone function

        //virtual function from CModelBase
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        void getGeneValuefromStruct(std::vector<double>&);
        void VarRangeStructRange(std::vector<VarRangeStruct>&);

        void createSupport(const Bitset &);
        void createMoleAdsorb(const Bitset &);

        void createSupportedCluster(const Bitset&);
        void createSupportSurface(const Bitset&);

        void createClusterSphere();
        void createSupportPlane();

        Bitset SupportBit();
        Bitset MoleAdsorbBit();
        Bitset ClusterBit();
        Bitset SupportSurfaceBit();
        //overload the father class's function
        void setPeriodicFramekwork(CATAZJUT::CConfigurationBase* mbf);

        LATT_DIRECTION latticeDirection();
        void setLatticeDirection(LATT_DIRECTION);

    protected:
        void IdentifyvacuumLayerDirection();
        void perceiveSupportSurface();
        void eliminateCloseContacts(double distanceCutOff=1.0);

    private:
        CModelSupport            *m_pSupport;
        CModelMoleculeAdsorbent  *m_pAdsorbMolecule;

        CModelMoleculeAdsorbent  *m_pClusterOnSupport;
        CModelSupport            *m_pPureSupportSurface;

        LATT_DIRECTION            m_latticeDirection;
        CATAZJUT::CPlane         *m_support_surface;
        CATAZJUT::CSphere        *m_pSphereCluster;

                   Bitset         m_BitbackupSupport;
                   Bitset         m_BitbackupAdsorbMolecule;
                   Bitset         m_BitbackupCluster;
                   Bitset         m_BitBackupSupportNoCluster;
};



}


#endif // CMODELCLUSTERLOADED2DSUPPORT_H
