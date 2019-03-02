/******************************************************************************
**
** Copyright (C) 2019-2031 Dr.Gui-lin Zhuang <glzhuang@zjut.edu.cn>
** All rights reserved.
**
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/
#include <cmath>
#include <Eigen/Dense>
#include "CUnitCell.h"
#include "Geometry.h"
#include "../IACS.h"
namespace CATAZJUT{

CUnitCell::CUnitCell()
:m_type(CATAZJUT::DEFINED::Bravais)
{
    this->m_a<<0,0,0;
    this->m_b<<0,0,0;
    this->m_c<<0,0,0;
    this->m_ScalingFactor=1.0;
}
CUnitCell::CUnitCell(Vector3 &a, Vector3 &b, Vector3 &c, double Scaling=1.0)
{
    this->m_a = a;
    this->m_b = b;
    this->m_c = c;
    this->m_type=CATAZJUT::DEFINED::Bravais;
    this->m_ScalingFactor=Scaling;
}
CUnitCell::CUnitCell(CUnitCell &othr)
{
    this->m_a = othr.avec();
    this->m_b = othr.bvec();
    this->m_c = othr.cvec();
    this->m_type = othr.unitCellType();
    this->m_ScalingFactor =othr.scalingFactor();
}
CUnitCell::~CUnitCell()
{
    //dtor
}
Vector3 CUnitCell::avec()const
{
    return this->m_a;
}
Vector3 CUnitCell::bvec()const
{
    return this->m_b;
}
Vector3 CUnitCell::cvec()const
{
    return this->m_c;
}
double CUnitCell::a()const
{
    return m_ScalingFactor * m_a.norm();
}
double CUnitCell::b()const
{
    return m_ScalingFactor * m_b.norm();
}
double CUnitCell::c()const
{
    return m_ScalingFactor * m_c.norm();
}
double CUnitCell::alfaRadians()const
{
    return CATAZJUT::Geometry::angleRadians(m_b,m_c);
}
double CUnitCell::alfa()const
{
    return CATAZJUT::Geometry::angle(m_b,m_c);
}
double CUnitCell::beta()const
{
   return CATAZJUT::Geometry::angle(m_a,m_c);
}
double CUnitCell::betaRadians()const
{
   return CATAZJUT::Geometry::angleRadians(m_a,m_c);
}
double CUnitCell::gama()const
{
   return CATAZJUT::Geometry::angle(m_a,m_b);
}
double CUnitCell::gamaRadians()const
{
   return CATAZJUT::Geometry::angleRadians(m_a,m_b);
}
double CUnitCell::volumn()const
{
    return std::pow(m_ScalingFactor,3) * m_a.dot(m_b.cross(m_c));
}
double CUnitCell::scalingFactor()const
{
    return this->m_ScalingFactor;
}
void CUnitCell::setscalingFactor(const double& m)
{
    this->m_ScalingFactor = m;
}
CATAZJUT::DEFINED::UnitcellType CUnitCell::unitCellType()
{
    return this->m_type;
}
void CUnitCell::setVec(const size_t index, const Vector3& m)
{
    switch(index)
    {
        case 0:
           this->m_a= m;
           break;
        case 1:
           this->m_b= m;
           break;
        case 2:
           this->m_c= m;
           break;
        default:
            break;
    }
}
void CUnitCell::setCellType(CATAZJUT::DEFINED::UnitcellType _type)
{
    this->m_type=_type;
}
/*
transtion matrix A
ax   bx   cx
ay   by   cy
az   bz   cz
*/
Eigen::Matrix<double, 3, 3> CUnitCell::MatrixOfBravaisLattice()
{
    Eigen::Matrix<double, 3, 3> tmpMat;
//    if(this->m_type==CATAZJUT::DEFINED::Bravais){
//        std::cout<<"OK,OK!"<<std::endl;
//    }
    if(this->m_type==CATAZJUT::DEFINED::Bravais){
        tmpMat.col(0)=m_a;
        tmpMat.col(1)=m_b;
        tmpMat.col(2)=m_c;
        tmpMat = tmpMat*m_ScalingFactor;

        //std::cout<<"Bravais:"<<tmpMat<<std::endl;

        return tmpMat;
    }else{
        CUnitCell* tempCell=this->toBravais();

        tmpMat.col(0)=(tempCell->scalingFactor())*tempCell->avec();
        tmpMat.col(1)=(tempCell->scalingFactor())*tempCell->bvec();
        tmpMat.col(2)=(tempCell->scalingFactor())*tempCell->cvec();
        delete tempCell;
        //std::cout<<"no Bravais:"<<avec()<<std::endl;
        return tmpMat;
    }

}
Eigen::Matrix<double, 3, 3> CUnitCell::NormilizedBravaisMatrix()
{
   Eigen::Matrix<double, 3, 3> tmpMat;
   tmpMat = MatrixOfBravaisLattice();
   tmpMat.col(0)=tmpMat.col(0).normalized();
   tmpMat.col(1)=tmpMat.col(1).normalized();
   tmpMat.col(2)=tmpMat.col(2).normalized();
   return tmpMat;
}
CUnitCell* CUnitCell::toReciprocal()
{
    if(this->m_type==CATAZJUT::DEFINED::Reciprocal)
        return this;
    else{
        CUnitCell *tempUnitCell = new CUnitCell();
        tempUnitCell->m_type = CATAZJUT::DEFINED::Reciprocal;
        tempUnitCell->m_ScalingFactor = 1/m_ScalingFactor;

        Eigen::MatrixXd *BravaisMat = new (Eigen::MatrixXd)(3,3);
        BravaisMat->col(0)=m_a;
        BravaisMat->col(1)=m_b;
        BravaisMat->col(2)=m_c;

        *BravaisMat = BravaisMat->inverse();

        tempUnitCell->m_a = 2*(CATAZJUT::constants::Pi)*BravaisMat->col(0);
        tempUnitCell->m_b = 2*(CATAZJUT::constants::Pi)*BravaisMat->col(1);
        tempUnitCell->m_c = 2*(CATAZJUT::constants::Pi)*BravaisMat->col(2);

        delete BravaisMat;
        return tempUnitCell;
    }
}
CUnitCell* CUnitCell::toBravais()
{
    if(this->m_type==CATAZJUT::DEFINED::Bravais)
        return this;
    else{
        CUnitCell *tempUnitCell = new CUnitCell();
        tempUnitCell->m_type = CATAZJUT::DEFINED::Bravais;
        tempUnitCell->m_ScalingFactor = 1/m_ScalingFactor;

        Eigen::MatrixXd *ReciprocalMat = new (Eigen::MatrixXd)(3,3);
        ReciprocalMat->col(0)=m_a;
        ReciprocalMat->col(1)=m_b;
        ReciprocalMat->col(2)=m_c;

        *ReciprocalMat = ReciprocalMat->inverse();

        tempUnitCell->m_a = 2*(CATAZJUT::constants::Pi)*ReciprocalMat->col(0);
        tempUnitCell->m_b = 2*(CATAZJUT::constants::Pi)*ReciprocalMat->col(1);
        tempUnitCell->m_c = 2*(CATAZJUT::constants::Pi)*ReciprocalMat->col(2);

        delete ReciprocalMat;
        return tempUnitCell;
    }
}
bool CUnitCell::operator == (CUnitCell& othr)
{
    if(this->m_type != othr.unitCellType())
        return false;
    return (m_a == othr.avec()) && (m_b == othr.bvec()) && \
           (m_c == othr.cvec()) && (m_ScalingFactor == othr.scalingFactor()) ;
}
CATAZJUT::DEFINED::CrystalSystemType CUnitCell::crystalSystem(double a,double b,double c,double alfa,double beta,double gama)
{
    if(isEqual(alfa,90.0) && isEqual(beta,90.0) && isEqual(gama,90.0)){
        if(isEqual(a,b,0.1) && isEqual(a,c,0.1))
            return CATAZJUT::DEFINED::Cubic;
        else if(isEqual(a,b,0.1) && !isEqual(b,c,0.1))
            return CATAZJUT::DEFINED::Tetragonal;
        else if(!isEqual(a,b,0.1)&& !isEqual(c,b,0.1)&& !isEqual(a,c,0.1))
            return CATAZJUT::DEFINED::Orthorhombic;
    }else if(isEqual(alfa,90.0) && isEqual(beta,90.0) && isEqual(gama,120.0)
             && isEqual(a,b,0.1) && !isEqual(c,b,0.1))
             return CATAZJUT::DEFINED::Hexagonal;
     else if(isEqual(a,b,0.1) && isEqual(a,c,0.1)&& isEqual(alfa,beta)
             && isEqual(gama,beta)&& !isEqual(gama,90.0))
             return CATAZJUT::DEFINED::Trigonal;
     else if(!isEqual(a,b,0.1)&& !isEqual(c,b,0.1)&& !isEqual(a,c,0.1)
             && isEqual(alfa,90.0) && isEqual(beta,90.0) && !isEqual(gama,90.0))
             return CATAZJUT::DEFINED::Monoclinic;
     else
             return CATAZJUT::DEFINED::Triclinic;

     return CATAZJUT::DEFINED::ERROR;
}
CATAZJUT::DEFINED::CrystalSystemType CUnitCell::crystalSystem()
{
    return crystalSystem(a(),b(),c(),alfa(),beta(),gama());
}
bool CUnitCell::orthoCoordinateSystem()
{
    return crystalSystem()== CATAZJUT::DEFINED::Cubic      ||
           crystalSystem()== CATAZJUT::DEFINED::Tetragonal ||
           crystalSystem()== CATAZJUT::DEFINED::Orthorhombic;
}
bool CUnitCell::isEqual(double m, double n,double error)
{
    return std::fabs(m-n) < error;
}
void CUnitCell::fromCellPara(double a1, double b1, double c1,double alfa1, double beta1, double gama1,char labellattice)
{
    CATAZJUT::DEFINED::CrystalSystemType   system_res = crystalSystem(a1,b1,c1,alfa1,beta1,gama1);
    double cosA,sinA,cosB,sinB,cosC;

    switch(system_res){
        case CATAZJUT::DEFINED::Triclinic:
              if(labellattice=='p' || labellattice =='P'){
                 cosA = std::cos(alfa1*CATAZJUT::constants::DegreesToRadians);
                 sinA = std::sin(alfa1*CATAZJUT::constants::DegreesToRadians);
                 cosB = std::cos(beta1*CATAZJUT::constants::DegreesToRadians);
                 sinB = std::sin(beta1*CATAZJUT::constants::DegreesToRadians);
                 cosC = std::cos(gama1*CATAZJUT::constants::DegreesToRadians);
                 //sinC = std::sin(gama1*CATAZJUT::constants::DegreesToRadians);
                 double tempB=b1*(cosC-cosA*cosB)/sinB;
                 double tempB2=std::sqrt(sinA*sinA-(cosC-cosA*cosB)*(cosC-cosA*cosB)/(sinB*sinB));

                 m_a<<a1*sinB, 0,      a1*cosB;
                 m_b<<  tempB, tempB2, b1*cosA;
                 m_c<<0,       0,      c1;
              }
              break;
        case CATAZJUT::DEFINED::Monoclinic:
            cosB = std::cos(beta1*CATAZJUT::constants::DegreesToRadians);
            sinB = std::sin(beta1*CATAZJUT::constants::DegreesToRadians);
            if(labellattice=='p' || labellattice =='P'){
                m_a<<a1*sinB, 0,   a1*cosB;
                m_b<<0,       b1,  0;
                m_c<<0,       0,   c1;
            }else if(labellattice=='C' || labellattice =='c'){
                m_a<<0.5*a1, -0.5*b1,   0;
                m_b<<0.5*a1,  0.5*b1,   0;
                m_c<<c1*cosB,      0,  c1*sinB;
            }
            break;
        case CATAZJUT::DEFINED::Orthorhombic:
            if(labellattice=='P' || labellattice =='p'){
                orthosimpleLattice(a1,b1,c1);
            }else if(labellattice=='I' || labellattice =='i'){
                orthobodyLattice(a1,b1,c1);
            }else if(labellattice=='F' || labellattice =='f'){
                orthofaceLattice(a1,b1,c1);
            }else if(labellattice=='C' || labellattice =='c'){
                m_a<<a1,-0.5*b1, 0;
                m_b<<a1, 0.5*b1, 0;
                m_c<<0,  0,      c1;
            }
            break;
        case CATAZJUT::DEFINED::Trigonal:
            if(labellattice=='R' || labellattice =='r'){
              cosA = std::cos(alfa1*CATAZJUT::constants::DegreesToRadians);
              double alta,kesai;
              alta = a1*std::sqrt((1+2*cosA)/3.0);
              kesai = a1*std::sqrt((1-cosA)*2.0/3.0);
                m_a<<0, alta, kesai;
                m_b<<-0.5*alta*std::sqrt(3.0), -0.5*alta, kesai;
                m_c<< 0.5*alta*std::sqrt(3.0),  -0.5*alta, kesai;
            }
            break;
        case CATAZJUT::DEFINED::Hexagonal:
            if(labellattice=='P' || labellattice =='p'){
                m_a<<-0.5*a1*std::sqrt(3.0), -0.5*a1, 0;
                m_b<< 0, a1, 0;
                m_c<< 0, 0,  c1;
            }
            break;
        case CATAZJUT::DEFINED::Tetragonal:
            if(labellattice=='P' || labellattice =='p'){
                orthosimpleLattice(a1,b1,c1);
            }else if(labellattice=='I' || labellattice =='i'){
                orthobodyLattice(a1,b1,c1);
            }
            break;
        case CATAZJUT::DEFINED::Cubic:
            if(labellattice=='P' || labellattice =='p'){
                orthosimpleLattice(a1,b1,c1);
            }else if(labellattice=='I' || labellattice =='i'){
                orthobodyLattice(a1,b1,c1);
            }else if(labellattice=='F' || labellattice =='f'){
                orthofaceLattice(a1,b1,c1);
            }
            break;
        default:
            break;
    }
    m_ScalingFactor=1.0;
}
void CUnitCell::orthosimpleLattice(double a1, double b1, double c1)
{
    m_a<<a1,0,0;
    m_b<<0,b1,0;
    m_c<<0,0, c1;
}
void CUnitCell::orthobodyLattice(double a1, double b1, double c1)
{
    m_a<<-0.5*a1, 0.5*b1,  0.5*c1;
    m_b<< 0.5*a1,-0.5*b1,  0.5*c1;
    m_c<< 0.5*a1, 0.5*b1, -0.5*c1;
}
void CUnitCell::orthofaceLattice(double a1, double b1, double c1)
{
   m_a<< 0,      0.5*b1,  0.5*c1;
   m_b<< 0.5*a1, 0,       0.5*c1;
   m_c<< 0.5*a1, 0.5*b1,  0;
}






}
