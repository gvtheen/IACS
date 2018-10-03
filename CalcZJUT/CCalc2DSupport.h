#ifndef CCALC2DSUPPORT_H
#define CCALC2DSUPPORT_H

#include <vector>
#include "../Util/Bitset.h"
#include "../GaZJUT/GaDeclaration.h"
#include "CCalcModeStruct.h"
#include "../Util/Point-Vector.h"

using GAZJUT::GeneVAR;
using util::Point3;
using util::Vector3;

namespace CATAZJUT{
  class CPeriodicFramework;
  class CPlane;
  class CAtom;
}

namespace CALCZJUT{

class CCalcMoleculeAdsorbent;
class CCalcSupportBase;
class CParameter;

class CCalc2DSupport:public CCalcModeStruct
{
    public:
        typedef enum LA_001{A_AXIS=0x30AB,
                            B_AXIS=0x30AC,
                            C_AXIS=0x30AD,
                          NONE_DIR=0x30AE}LATT_DIRECTION;

        CCalc2DSupport(CParameter*);
        virtual ~CCalc2DSupport();

        //virtual function from CCalcModeStruct
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        void getGeneValuefromStruct(std::vector<double>&);
        void GeneVARRange(std::vector<GeneVAR>&);

        void createSupport(Bitset &);
        void createMoleAdsorb(Bitset &);
        //overload the father class's function
        void setPeriodicFramekwork(CATAZJUT::CPeriodicFramework* mbf);

//        void backUpStructure();
//        void fromBackupToCurrent();

        void IdentifyvacuumLayerDirection();
        void perceiveSupportSurface();
        void eliminateCloseContacts(double distanceCutOff=1.0);

    protected:

    private:
        CCalcSupportBase        *m_pSupport;
        CCalcMoleculeAdsorbent  *m_pAdsorbMolecule;
        LATT_DIRECTION           m_latticeDirection;
        CATAZJUT::CPlane        *m_support_surface;

        Bitset        m_BitbackupSupport;
        Bitset        m_BitbackupAdsorbMolecule;

};

}
#endif // CCalc2DSupport_H
