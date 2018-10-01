#ifndef CSPHERE_H
#define CSPHERE_H
#include <Eigen/Dense>
#include "../Util/Point-Vector.h"
using util::Point3;

namespace CATAZJUT{

class CConfigurationBase;

class CSphere
{
    public:
        CSphere();
        CSphere(CConfigurationBase*);
        CSphere(CSphere&);

        CSphere* Clone();
        void     CreateSphere();

        void setConfiguration(CConfigurationBase*);
        CConfigurationBase* configuration()const;
        double  Radius() const;
        Point3  SphereCenter() const;
        Eigen::Vector4d Equation()const;

        Point3  toPolarCoordinate(const Point3&);

        Point3  geneCartesianCoordAt(const Point3& posPolarPoint);
        virtual ~CSphere();

    protected:
        size_t maxDist(Point3);
        bool checkPointIs(const Point3&,const Point3&);
    private:
        CConfigurationBase *m_pConfiguration;
        //sphere (x-x0)**2 + (y-y0)**2 + (z-z0)**2 =R0**2;
        // x0=m_Equation[0]; y0=m_Equation[1];z0=m_Equation[2];R=m_Equation[3]
        Eigen::Vector4d m_Equation;
};


}
#endif // CSPHERE_H
