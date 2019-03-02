#ifndef CSPHERE_H
#define CSPHERE_H
#include <vector>
#include <Eigen/StdVector>
#include "../Util/Point-Vector.h"
#include "../Util/log.hpp"

using util::Log;
using util::Point3;

namespace CATAZJUT{

class CSphere
{
    public:
        CSphere();
        CSphere(const std::vector<Point3,Eigen::aligned_allocator<Point3>>&);
        CSphere(CSphere&);

        CSphere* Clone();
        void     CreateSphere();

//        void setConfiguration(CConfigurationBase*);
//        CConfigurationBase* configuration()const;
        double  Radius() const;
        Point3  SphereCenter() const;
        util::Vector4 Equation()const;

        Point3  toPolarCoordinate(const Point3&);

        Point3  CartesianCoordAtGeneOf(const Point3& posPolarPoint);
        virtual ~CSphere();

        void output();
    protected:
        size_t maxDist(const Point3&);
        bool checkPointIs(const Point3&,const Point3&);
    private:

        util::Vector4 m_Equation;

        std::vector<Point3,Eigen::aligned_allocator<Point3>>  m_PointsMat;
        //sphere (x-x0)**2 + (y-y0)**2 + (z-z0)**2 =R0**2;
        // x0=m_Equation[0]; y0=m_Equation[1];z0=m_Equation[2];R=m_Equation[3]

};


}
#endif // CSPHERE_H
