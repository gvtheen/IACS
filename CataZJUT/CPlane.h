#ifndef CPLANE_H
#define CPLANE_H
#include <Eigen/Dense>
#include "Point-Vector.h"
namespace CATAZJUT{

class CPlane
{
    public:
        CPlane();
        CPlane(Eigen::Vector4d);
        CPlane(Eigen::MatrixXd* currentMat);
        CPlane(CPlane&);
        virtual ~CPlane();

        void            CreatePlane(Eigen::MatrixXd* currentMat);
        double          Distance(Point3);
        Eigen::Vector4d Equation();
        void SetEquation(Eigen::Vector4d temp);

        Eigen::Vector3d  NormalLine();
        void SetNormalLine(Eigen::Vector3d temp);

        Circle3D       CircleIncludingPoints();
        void SetCircleOfIncludingPoints(Circle3D);
        Point3          Projection(Point3); //
        Eigen::MatrixXd* Projection(Eigen::MatrixXd*);
        //if dis>0 new plan locate above old plane along normal line
        //Otherwise, it locate below old plane.
        CPlane*         PlaneWithDistance(double);
        CPlane*         Clone();

        Point3          PointInCircleFromGene(double Height,double Angle_degrees,double Radius_ratio);
                        //Angle_degrees =[0,360] deg, Radius_ratio=[0,1]

        //Varible
        Eigen::MatrixXd* m_pPointsToPlane;
    protected:
        void            CreateCircleIncludingPoints(Eigen::MatrixXd);
        double          FindDisMax(Eigen::MatrixXd,Point3,int&);
       // bool            IsOKCreatedPlane(std::vector<int> *defaultPoint);
    private:
        Eigen::Vector4d *m_pEquation;
        Eigen::Vector3d *m_pNormalLine;
        Circle3D        *m_pCircleIncludingPoints;
        // Plane equation: ax + by + cz + d=0;
        // a=m_Equation[0]; b=m_Equation[1];c=m_Equation[2];d=m_Equation[3]
        //it must be normalize
};

}
#endif // CPLANE_H
