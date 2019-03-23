#ifndef CModel2DSupport_H
#define CModel2DSupport_H

#include <vector>
#include "../Util/Bitset.h"
#include "CModelBase.h"
#include "../Util/Point-Vector.h"
#include "../IACS.h"

using IACSZJUT::VarRangeStruct;
using util::Point3;
using util::Vector3;

namespace CATAZJUT{
  class CConfigurationBase;
  class CPlane;
  class CAtom;
}

namespace CALCZJUT{

class CModelMoleculeAdsorbent;
class CModelSupport;
class CParameter;

class CModel2DSupport:public CModelBase
{
    public:
        typedef enum LA_001{A_AXIS=0x30AB,
                            B_AXIS=0x30AC,
                            C_AXIS=0x30AD,
                          NONE_DIR=0x30AE}LATT_DIRECTION;

        CModel2DSupport(CParameter*,
                        CATAZJUT::CConfigurationBase**,
                        size_t index);
        virtual ~CModel2DSupport();

        CModelBase* clone();

        //virtual function from CModelBase
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        void getGeneValuefromStruct(std::vector<double>&);
        void VarRangeStructRange(std::vector<VarRangeStruct>&);

        void createSupport(const Bitset &);
        void createMoleAdsorb(const Bitset &);

        Bitset SupportBit();
        Bitset MoleAdsorbBit();
        //overload the father class's function
        void setPeriodicFramekwork(CATAZJUT::CConfigurationBase* mbf);
        void outputStructureToFile();

        CATAZJUT::CPlane* supportSurfacePlane();
        void setSupportSurfacePlane(CATAZJUT::CPlane*);

        LATT_DIRECTION latticeDirection();
        void setLatticeDirection(LATT_DIRECTION);

        void IdentifyvacuumLayerDirection();
        void perceiveSupportSurface();
        void eliminateCloseContacts(double distanceCutOff=1.0);

    protected:

    private:
        CModelSupport            *m_pSupport;
        CModelMoleculeAdsorbent  *m_pAdsorbMolecule;
        LATT_DIRECTION            m_latticeDirection;
        CATAZJUT::CPlane         *m_support_surface;

                   Bitset         m_BitbackupSupport;
                   Bitset         m_BitbackupAdsorbMolecule;

};

}
#endif // CModel2DSupport_H
