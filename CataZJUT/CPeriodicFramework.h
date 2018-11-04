#ifndef CPERIODICFRAMEWORK_H
#define CPERIODICFRAMEWORK_H
#include "CConfigurationBase.h"

namespace CALCZJUT{
  class CModel2DSupport;
  class CIOCar;
  class CExeVASP;
  class CParameter;
 }
namespace CATAZJUT{

class CUnitCell;
class CAtom;
class CBondTolerance;
class CCartesianCoordinates;
//


class CPeriodicFramework:public CConfigurationBase
{
    public:
        CPeriodicFramework(CALCZJUT::CParameter*);
        virtual ~CPeriodicFramework();
        CPeriodicFramework(CPeriodicFramework&);
        CPeriodicFramework(CConfigurationBase&);

        CPeriodicFramework* clone();
        virtual CFractionCoordinates* Fractioncoordinates();
        CUnitCell* unitcell();

        std::string SymmetrySymbol();
        void perceiveBonds();
    protected:

    private:
        friend class CFractionCoordinates;
        friend class CAtom;
        friend class CBondTolerance;
        friend class CUnitCell;
        friend class CCartesianCoordinates;
        friend class CConfigurationPrivateData;

        friend class CModel2DSupport;
        friend class CIOPoscar;
        friend class CIOCellFile;
        friend class CExeVASP;

    private:
        CUnitCell *m_pUnitCell;
        //mutable CFractionCoordinates        *m_pFraction;

};


}
#endif // CSUPPORT2D_H
