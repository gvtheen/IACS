#ifndef CCREATEPLANE_H
#define CCREATEPLANE_H
#include <Eigen/Dense>
#include <vector>
#include "Point-Vector.h"
#include "CPlane.h"
namespace CATAZJUT{

class CCrystalPlanes
{
    public:
        CCrystalPlanes();
        CCrystalPlanes(Eigen::MatrixXd*,double m_Value=2.8);
        CCrystalPlanes(CCrystalPlanes&);
        virtual ~CCrystalPlanes();

        CCrystalPlanes* Clone();
        void CreateCrystalPlane();
        bool CheckIsPlane(CPlane,std::vector<size_t>&,std::vector<bool>&);

        void SetLatticPlane(std::vector<CPlane*>*);
        std::vector<CPlane*>* LatticePlane();

        Eigen::MatrixXd* PointsOfMat();
        void SetPointsMat(Eigen::MatrixXd*);

        double DistanceCutoff();
          void SetDistanceCutoff(double);

    protected:
        void RemoveRow(Eigen::MatrixXd&, size_t);
        void RemoveColumn(Eigen::MatrixXd&, size_t);
        bool IsOnLine(Eigen::MatrixXd&);
    private:
              Eigen::MatrixXd *m_pPointsMat;
         std::vector<CPlane*> *m_pPlane;
std::vector<Eigen::MatrixXd*> *m_pPointsInIndividualPlanes;
                        double mDistance_Cutoff;
};

}
#endif // CCREATEPLANE_H
