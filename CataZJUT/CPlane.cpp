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
#include<cmath>
#include <iostream>
#include "CPlane.h"
#include "Constant.h"
namespace CATAZJUT{

CPlane::CPlane()
{
    m_pEquation=nullptr;
    m_pNormalLine=nullptr;
    m_pCircleIncludingPoints = nullptr;
}
CPlane::CPlane(Eigen::Vector4d currentPlane)
{
    m_pEquation     = new (Eigen::Vector4d);
    *m_pEquation = currentPlane;
    m_pNormalLine = new (Eigen::Vector3d);
   (*m_pNormalLine)(0) = currentPlane(0);
   (*m_pNormalLine)(1) = currentPlane(1);
   (*m_pNormalLine)(2) = currentPlane(2);
}
CPlane::CPlane(CPlane& currentPlane)
{
      m_pEquation = new (Eigen::Vector4d);
    m_pNormalLine = new (Eigen::Vector3d);
     *m_pEquation = currentPlane.Equation();
    *m_pNormalLine= currentPlane.NormalLine();
}
CPlane::CPlane(Eigen::MatrixXd* currentMat)
{
    this->CreatePlane(currentMat);
}

void CPlane::CreatePlane(Eigen::MatrixXd* CurrentPointMat)
{
    assert(CurrentPointMat);
    assert(CurrentPointMat->cols()==3);

    //int rowNum_org =CurrentPointMat->rows();
     m_pPointsToPlane = new (Eigen::MatrixXd)(*CurrentPointMat);


    if(m_pEquation==nullptr)
      m_pEquation= new (Eigen::Vector4d);

    if(m_pNormalLine==nullptr)
      m_pNormalLine = new (Eigen::Vector3d);
    //

        int rowNum = m_pPointsToPlane->rows();
        Eigen::MatrixXd* xy1Mat = new (Eigen::MatrixXd)(rowNum,3);
          Eigen::MatrixXd* zMat = new (Eigen::MatrixXd)(rowNum,1);
        xy1Mat->setOnes();
        xy1Mat->col(0) = m_pPointsToPlane->col(0);
        xy1Mat->col(1) = m_pPointsToPlane->col(1);
          zMat->col(0) = m_pPointsToPlane->col(2);          //ax+by+c=z    res = [a,b,c]^T
        Eigen::Vector3d res = xy1Mat->colPivHouseholderQr().solve(*zMat);

          (*m_pEquation)<< res(0),res(1),-1.0,res(2);     //c
        (*m_pNormalLine)<< res(0),res(1),-1.0;
        this->m_pNormalLine->normalize();
    delete xy1Mat;
    delete zMat;
}
//bool CPlane::IsOKCreatedPlane(std::vector<int> *defaultPoint)
//{
//
//}
double CPlane::Distance(const Point3& P1) // projection of P0P1 vector on normal line of plane.
{
    // if distance<0, projection vector of P0P1 vector on normal line lie contrary direction with
    // normal line. Thus, P1 point lie below the plane
    // Otherwise  P1 point lie above the plane.
    Point3  P0;
    P0<<0,0,(*m_pEquation)(3);   // any point on the plane.
    Vector3 P0P1 = P1 - P0;
    return  P0P1.dot(*m_pNormalLine)/ m_pNormalLine->norm();
}
Point3 CPlane::Projection(Point3 P1)
{
   Point3  P0;
   P0<<0,0,(*m_pEquation)(3);   // Any Point on the plane
   Vector3 V = P1 - P0;
   Vector3 V_Para = V.dot(*m_pNormalLine) * (*m_pNormalLine);  //V
   Vector3 V_Vert = V - V_Para;
   return  P0 + V_Vert;
}
Eigen::MatrixXd* CPlane::Projection(Eigen::MatrixXd* currentMat)
{
   Eigen::MatrixXd* newMat = new (Eigen::MatrixXd)(currentMat->rows(),currentMat->cols());
   int row_Num = currentMat->rows();
   for(int i=0;i<row_Num;i++)
      newMat->row(i)= Projection(currentMat->row(i));
   return newMat;
}
CPlane* CPlane::Clone()
{
    return new CPlane(*this);
}
CPlane* CPlane::PlaneWithDistance(double dis)
{
    CPlane* newPlane = this->Clone();
//     Plane equation: ax + by + cz + d=0;
//     a=m_Equation[0]; b=m_Equation[1];c=m_Equation[2];d=m_Equation[3]
    Vector3 temp;
    temp<<(*m_pEquation)(0),(*m_pEquation)(1),(*m_pEquation)(2);
    Eigen::Vector4d tempVect;
    tempVect = newPlane->Equation();
    tempVect(3)= (*m_pEquation)(3) - dis*temp.norm();
    newPlane->SetEquation(tempVect);
    //if dis>0 new plan locate above old plane along normal line
    //Otherwise, it locate below old plane.
    return newPlane;
}
void CPlane::CreateCircleIncludingPoints(Eigen::MatrixXd ScatPointCoord)
{
    if(ScatPointCoord.cols()!=3)
        ;
    double ans  =1e12,eps=1e-12,R;
    Point3 CirclePoint= ScatPointCoord.row(0);
    Point3 temPoint;
    int pos;
    double step =FindDisMax(ScatPointCoord,CirclePoint,pos);
    while(step>eps)
    {
        R=FindDisMax(ScatPointCoord,CirclePoint,pos);
        if(R < ans)
            ans=R;
        temPoint<<ScatPointCoord(pos,0),
                  ScatPointCoord(pos,1),
                  ScatPointCoord(pos,2);
        temPoint=(CirclePoint - temPoint);
        temPoint.normalize();
        CirclePoint = CirclePoint + step*temPoint;
        step=step*0.9;
    }
    if( m_pCircleIncludingPoints==nullptr )
        m_pCircleIncludingPoints = new Circle3D;

    Vector3 U;
    if((*m_pNormalLine)(0)!=0 || (*m_pNormalLine)(1)!=0)
        U<<(*m_pNormalLine)(1),-1*(*m_pNormalLine)(0),0;
    else
        U<<0,-1*(*m_pNormalLine)(2),(*m_pNormalLine)(1);
    U.normalize();
    m_pCircleIncludingPoints->col(0) = CirclePoint;
    m_pCircleIncludingPoints->col(1) = R*U;
    m_pCircleIncludingPoints->col(2) = R*(m_pNormalLine->cross(U));;

}
double CPlane::FindDisMax(Eigen::MatrixXd ScatPointCoord,Point3 cirPoint,int& index)
{
    double res=0,temp;
    index=0;
    int row_num=ScatPointCoord.rows();
    Point3 tem_P;
    for(int i=0;i<row_num;i++)
    {
        tem_P<<ScatPointCoord(i,0),ScatPointCoord(i,1),ScatPointCoord(i,2);
        temp = (tem_P - cirPoint).norm();
        if(temp>res){
            res=temp;
            index = i;
        }
    }
    return res;
}
Point3 CPlane::PointInCircleFromGene(double Height,double Radius_ratio,double Angle_Radian) // Gene = Angle_degrees
{
    CPlane *newPlane = std::move(PlaneWithDistance(Height));

    Eigen::MatrixXd* ProjPointsInNewPlane = std::move(newPlane->Projection(this->m_pPointsToPlane));
    newPlane->CreateCircleIncludingPoints(*ProjPointsInNewPlane);

    Vector3 Angle_Vect;
//    double Angle_Radian = Angle_degrees*CATAZJUT::constants::DegreesToRadians;

    Angle_Vect<<1.0,std::cos(Angle_Radian),std::sin(Angle_Radian);

    Point3 PointOnCircle = (*(newPlane->m_pCircleIncludingPoints))*Angle_Vect;
    Point3 CircleCenter = newPlane->m_pCircleIncludingPoints->col(0);
    return CircleCenter + Radius_ratio*(PointOnCircle - CircleCenter);
}
Eigen::Vector4d CPlane::Equation()
{
    //if(this->m_pEquation!=nullptr)
       return *(m_pEquation);
}
void CPlane::SetEquation(Eigen::Vector4d  temp)
{
   if(this->m_pEquation==nullptr)
       m_pEquation = new (Eigen::Vector4d);
   *(m_pEquation)=temp;
}
Eigen::Vector3d CPlane::NormalLine()
{
    //if(this->m_pNormalLine!=nullptr)
       return *(m_pNormalLine);
}
void CPlane::SetNormalLine(Eigen::Vector3d  temp)
{
   if(this->m_pNormalLine==nullptr)
       m_pNormalLine = new (Eigen::Vector3d);
   *(m_pNormalLine)=temp;
}
Circle3D CPlane::CircleIncludingPoints()
{
    //if(this->m_pCircleIncludingPoints!=nullptr)
       return *m_pCircleIncludingPoints;
}
void CPlane::SetCircleOfIncludingPoints(Circle3D currCircle)
{
    if( m_pCircleIncludingPoints==nullptr )
        m_pCircleIncludingPoints = new Circle3D;
    *m_pCircleIncludingPoints = currCircle;

}
CPlane::~CPlane()
{
    if(m_pEquation!=nullptr)
        delete m_pEquation;
    if(m_pNormalLine!=nullptr)
        delete m_pNormalLine;
    if(m_pCircleIncludingPoints!=nullptr)
        delete m_pCircleIncludingPoints;
    if(m_pPointsToPlane!=nullptr)
        delete m_pPointsToPlane;
}

//Circle3D CPlane::CircleEquation(Point3 centerP,double R)
//{
//    Vector3 U;
//    if((*m_pNormalLine)(0)!=0 || (*m_pNormalLine)(1)!=0)
//        U<<(*m_pNormalLine)(1),-1*(*m_pNormalLine)(0),0;
//    else
//        U<<0,-1*(*m_pNormalLine)(2),(*m_pNormalLine)(1);
//    U.normalize();
//    Circle3D MyCircle;
//    MyCircle.col(0) = centerP;
//    MyCircle.col(1) = R*U;
//    MyCircle.col(2) = R*(m_pNormalLine->cross(U));
//    (x,y,z)^T = MyCircle.col(0) + cost*MyCircle.col(1) + sint*MyCircle.col(2);
//    Or
//    (x,y,z)^T = MyCircle*(1,cost,sint)^T
//*/
//    return MyCircle;
//}

}
