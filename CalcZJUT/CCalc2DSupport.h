#ifndef CCALC2DSUPPORT_H
#define CCALC2DSUPPORT_H

#include <vector>
#include "Bitset.h"
#include "GaDeclaration.h"
#include "CCalcModeStruct.h"
#include "Point-Vector.h"
using GAZJUT::GENEVAR;
using CATAZJUT::Point3;
using CATAZJUT::Vector3;
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
        std::vector<double>*  getGeneValuefromStruct()const;
        std::vector<GENEVAR>* GeneVarRange();
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
