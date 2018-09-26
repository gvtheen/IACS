#ifndef CPERIODICFRAMEWORK_H
#define CPERIODICFRAMEWORK_H
#include "CConfigurationBase.h"

namespace CALCZJUT{
  class CCalc2DSupport;
  class CIOCar;
  class CCalcVASP;
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
        virtual CFractionCoordinates* Fractioncoordinates();
        CUnitCell* unitcell();
    protected:

    private:
        friend class CFractionCoordinates;
        friend class CAtom;
        friend class CBondTolerance;
        friend class CUnitCell;
        friend class CCartesianCoordinates;
        friend class CConfigurationPrivateData;

        friend class CCalc2DSupport;
        friend class CIOPoscar;
        friend class CIOCellFile;
        friend class CCalcVASP;

    private:
        CUnitCell *m_pUnitCell;

};


}
#endif // CSUPPORT2D_H
