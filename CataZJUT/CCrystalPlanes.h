#ifndef CCREATEPLANE_H
#define CCREATEPLANE_H
#include <Eigen/Dense>
#include <vector>
#include "../Util/Point-Vector.h"
#include "CPlane.h"
#include "../Util/Bitset.h"
/*
   This class main function is to identify crystal plane of the known cluster.


*/
using util::Point3;

namespace CATAZJUT{

class CConfigurationBase;

class CCrystalPlanes
{
    public:
        CCrystalPlanes();
        CCrystalPlanes(std::vector<size_t>&,CConfigurationBase*);
        CCrystalPlanes(CCrystalPlanes&);
        virtual ~CCrystalPlanes();

        CCrystalPlanes* Clone();
        void CreateCrystalPlane();

        void SetLatticPlane(std::vector<CPlane*>&);
        std::vector<CPlane*>& LatticePlane();

        Eigen::MatrixXd* PointsOfMat();
        void SetPointsMat(Eigen::MatrixXd*);

        CConfigurationBase* ConfigurationBase();
        void setConfigurationBase( CConfigurationBase*);

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
        bool IsInclude3PointsInExitPlane(size_t i,size_t j,size_t k);

        bool IsAdjecentPointer(size_t i,size_t j,size_t k);
        bool CheckIsCrystalPlane(CPlane&,std::vector<size_t>&,util::Bitset&);
    private:
              Eigen::MatrixXd*  m_PointsMat;
         std::vector<CPlane*>   m_Plane;
std::vector<Eigen::MatrixXd*>   m_PointsInIndividualPlanes;
                        double  mDistance_Cutoff=2.8;
            CConfigurationBase *m_pCConfigurationBase;
            std::vector<size_t> m_indexInConfig;
};

}
#endif // CCREATEPLANE_H
