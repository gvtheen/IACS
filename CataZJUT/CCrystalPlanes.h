#ifndef CCREATEPLANE_H
#define CCREATEPLANE_H
#include <Eigen/Dense>
#include <vector>
#include "../Util/Point-Vector.h"
#include "CPlane.h"
/*
   This class main function is to identify crystal plane of the known cluster.


*/
using util::Point3;

namespace CATAZJUT{

class CCrystalPlanes
{
    public:
        CCrystalPlanes();
        CCrystalPlanes(std::vector<Point3>&,double m_Value=2.8);
        CCrystalPlanes(CCrystalPlanes&);
        virtual ~CCrystalPlanes();

        CCrystalPlanes* Clone();
        void CreateCrystalPlane();
        bool CheckIsPlane(CPlane&,std::vector<size_t>&,std::vector<bool>&);

        void SetLatticPlane(std::vector<CPlane*>&);
        std::vector<CPlane*>& LatticePlane();

        Eigen::MatrixXd* PointsOfMat();
        void SetPointsMat(Eigen::MatrixXd*);

        size_t crystalPlaneNum();

        Point3 CartesianCoordinateAtGene(size_t crystal_Plane_No,double height,double R_radio,double thea);

        double DistanceCutoff();
          void SetDistanceCutoff(double);

        CPlane* operator[](size_t index);

        void outputCrystalPlane();

    protected:
        void RemoveRow(Eigen::MatrixXd&, size_t);
        void RemoveColumn(Eigen::MatrixXd&, size_t);
        bool IsOnLine(Eigen::MatrixXd&);

    private:
              Eigen::MatrixXd*  m_PointsMat;
         std::vector<CPlane*>   m_Plane;
std::vector<Eigen::MatrixXd*>   m_PointsInIndividualPlanes;
                        double  mDistance_Cutoff;
};

}
#endif // CCREATEPLANE_H
