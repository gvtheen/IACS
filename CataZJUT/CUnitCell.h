#ifndef CUNITCELL_H
#define CUNITCELL_H
#include <vector>
#include "Point-Vector.h"
namespace CATAZJUT{

class CUnitCell
{
    public:
        enum TYPE{
           Reciprocal,
           Bravais
        };
        typedef enum crystal_001{
           Cubic=0x40A1,
           Hexagonal=0x40A2,
           Trigonal=0x40A3,
           Tetragonal=0x40A4,
           Orthorhombic=0x40A5,
           Monoclinic=0x40A6,
           Triclinic=0x40A17
        }CRYSTALSYSTEM;
        CUnitCell();
        CUnitCell(Vector3 &, Vector3 &, Vector3 &,double );
        CUnitCell(CUnitCell &othr);
        virtual ~CUnitCell();

        void setVec(const size_t index, const Vector3&);
        void setscalingFactor(const double& m);
        Vector3 avec()const;
        Vector3 bvec()const;
        Vector3 cvec()const;

        double a()const;
        double b()const;
        double c()const;
        double alfa()const;
        double alfaRadians()const;
        double betaRadians()const;
        double beta()const;
        double gamaRadians()const;
        double gama()const;
        double scalingFactor()const;
        double volumn()const;

        Eigen::Matrix<double, 3, 3> MatrixOfBravaisLattice();
        Eigen::Matrix<double, 3, 3> NormilizedBravaisMatrix();

        CUnitCell* toReciprocal();
        CUnitCell* toBravais();
        void fromCellPara(double a, double b, double c,double alfa, double beta, double gamma,char labellattic);
        CRYSTALSYSTEM crystalSystem(double a,double b,double c,double alfa,double beta,double gama);
        CRYSTALSYSTEM crystalSystem();

        bool orthoCoordinateSystem();

        bool operator == (CUnitCell& othr);

        TYPE   m_type;
    protected:
        bool isEqual(double m, double n,double error=1.0);
        void orthosimpleLattice(double a1, double b1, double c1);
        void orthobodyLattice(double a1, double b1, double c1);
        void orthofaceLattice(double a1, double b1, double c1);
    private:
        Vector3 m_a;            //basic vector a<<ax,ay,az
        Vector3 m_b;            //basic vector b<<bx,by,bz
        Vector3 m_c;            //basic vector c<<cx,cy,cz
         double m_ScalingFactor;
};

}
#endif // CUNITCELL_H
