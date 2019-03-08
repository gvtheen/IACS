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
#include "CSphere.h"
#include "CPlane.h"
#include "Constant.h"
#include "Geometry.h"
#include "CCartesianCoordinates.h"
#include "CConfigurationBase.h"
namespace CATAZJUT{

CSphere::CSphere()
{
    //ctor
}
CSphere::CSphere(Eigen::MatrixXd* currentPoints)
{

    assert(currentPoints->rows()>0);

    this->m_PointsMat=new (Eigen::MatrixXd)(*currentPoints);

}
CSphere::CSphere(CSphere& mht)
{
    this->m_Equation = new util::Vector4(mht.Equation());
}
CSphere::~CSphere()
{
    delete m_PointsMat;
    delete this->m_Equation;
}
/*
    identify the sphere equation in the given discrete points ( cluster )
    by annealing method
*/
void CSphere::CreateSphere()
{
   double ans  =1e12,eps=1e-12,R,step;
   //CCartesianCoordinates* currentCart=this->m_pConfiguration->coordinates();

   Point3 P_pos,tmpP;
   P_pos<<0,0,0;
   size_t pos;
   step=100;
   #ifdef DEBUG
            Log::Debug<<"CSphere::CreateSphere():" <<m_PointsMat->rows()<< std::endl;
   #endif
   while(step>eps)
   {
       pos = maxDist(P_pos);
       tmpP<<(*m_PointsMat)(pos,0),(*m_PointsMat)(pos,1),(*m_PointsMat)(pos,2);
       R = CATAZJUT::Geometry::distance(P_pos,tmpP);
       if(R<ans)
          ans=R;
       if(R==0)
          R=R+eps;
       P_pos = P_pos + (step/R)*(tmpP-P_pos);
       step=step*0.98;
   }
   #ifdef DEBUG
                     Log::Debug<<"2-CSphere::CreateSphere():" <<R<< std::endl;
   #endif
   m_Equation = new (util::Vector4);
   (*m_Equation)<<P_pos[0],P_pos[1],P_pos[2],R;
}
/*
    sub method in annealing method
*/
size_t CSphere::maxDist(const Point3& p1)
{
    double res=0,temp;
    size_t res_index=0;
    Point3 tmpPonter;

    for(size_t i=0;i<m_PointsMat->rows();i++)
    {
        tmpPonter<<(*m_PointsMat)(i,0),(*m_PointsMat)(i,1),(*m_PointsMat)(i,2);
        temp=CATAZJUT::Geometry::distance(p1,tmpPonter);
        if(temp>res)
        {
            res_index=i;
            res= temp;
        }
    }
    return res_index;
}
util::Vector4 CSphere::Equation()const
{
    return *m_Equation;
}
//void CSphere::setConfiguration(CConfigurationBase* configure)
//{
//    this->m_pConfiguration=configure;
//}
//CConfigurationBase* CSphere::configuration()const
//{
//    return this->m_pConfiguration;
//}
double CSphere::Radius() const
{
    return (*m_Equation)[3];
}
Point3 CSphere::SphereCenter() const
{
    Point3 temp;
    temp<<(*m_Equation)[0],(*m_Equation)[1],(*m_Equation)[2];
    return temp;
}
/*
 convert cartesian coordinate to polar coordinate.
 return one point ( r, phi, thea )  angle unit:  Radians

 original point in polar coordinate system is center of sphere.
    >>Point3 newCartCoord = CartCoord - SphereCenter();

*/
Point3 CSphere::toPolarCoordinate( const Point3& CartCoord)  // res= r,phi,thea;  unit:  Radians
{
    Point3 newCartCoord = CartCoord - SphereCenter();

    double r = newCartCoord.norm();
    Point3 res;

    double x=newCartCoord[0], y=newCartCoord[1], z=newCartCoord[2];

    if (r==0)
    {
        res<<0,0,0;
        return res;
    }else if( x==0 && y==0 ){
        if(z>0)
           res<<z,0,0;
        else
           res<<z,0,180;
        return res;
    }
    res[0]=r;
    if(y>=0)
      res[1]= std::acos(x/std::sqrt(x*x+y*y));
    else
      res[1]= (360 - CATAZJUT::constants::RadiansToDegrees*std::acos(x/std::sqrt(x*x+y*y)))*\
              CATAZJUT::constants::DegreesToRadians;
    res[2]=std::acos(z/r);
    return res;
}
/*
  Gene polar coordinate on the point ( r,phi,thea  ) with the range of deta_radian;
   r : the distance between adsorbent molecule and surface of sphere.
  phi: [0,2*PI]
 thea: [0,PI]

 The function is to convert the polar coordinate gene to cartesian coordinate.
*/
Point3 CSphere::CartesianCoordAtGeneOf(const Point3& posPolarPoint)
{
    double maxDist=0;
    std::vector<size_t> res_Point;
    Point3 CurrentpolarPoint,tmpPonter,tmpPonter_1;
    // get all pointer from cluster configuration

    for(size_t i=0;i<m_PointsMat->rows();i++)
    {
        tmpPonter<<(*m_PointsMat)(i,0),(*m_PointsMat)(i,1),(*m_PointsMat)(i,2);
        CurrentpolarPoint = toPolarCoordinate( tmpPonter );
        if(checkPointIs(CurrentpolarPoint,posPolarPoint)==true)
        {
            res_Point.push_back(i);
            if(CurrentpolarPoint[0]>maxDist)
                maxDist = CurrentpolarPoint[0];
        }
    }
    std::vector<size_t>::iterator it;
    double deta_dist = 0.5;
    for(it=res_Point.begin();it!=res_Point.end();)
    {
        tmpPonter<<(*m_PointsMat)(*it,0),(*m_PointsMat)(*it,1),(*m_PointsMat)(*it,2);
        CurrentpolarPoint = toPolarCoordinate(tmpPonter);
        if(CurrentpolarPoint[0]< maxDist -deta_dist )
          it =  res_Point.erase(it);
        else
          it++;
    }
    if(res_Point.size()>=3){
       Eigen::MatrixXd matPoint(res_Point.size(),3);
       for(size_t i=0;i<res_Point.size();i++)
           matPoint.row(i)=m_PointsMat->row(res_Point[i]);

       CPlane *newPlane = new CPlane(matPoint);
       Point3 SphereCenterP=SphereCenter();
       Point3 projectPointOnPlane = newPlane->Projection(SphereCenterP);
       Point3 temp=projectPointOnPlane-SphereCenterP;
       delete newPlane;
       return SphereCenterP + temp*( 1.0 + temp.norm()/posPolarPoint[0] );
    }else if (res_Point.size()==2){
       Point3 SphereCenterP=SphereCenter();
       tmpPonter<<(*m_PointsMat)(res_Point[0],0),(*m_PointsMat)(res_Point[0],1),(*m_PointsMat)(res_Point[0],2);
       tmpPonter_1<<(*m_PointsMat)(res_Point[1],0),(*m_PointsMat)(res_Point[1],1),(*m_PointsMat)(res_Point[1],2);
       Point3 temp = CATAZJUT::Geometry::midpoint(tmpPonter,tmpPonter_1);
       temp=temp-SphereCenterP;
       return SphereCenterP + temp*( 1.0 + temp.norm()/posPolarPoint[0] );
    }else{
       Point3 SphereCenterP=SphereCenter();

       tmpPonter<<(*m_PointsMat)(res_Point[0],0),(*m_PointsMat)(res_Point[0],1),(*m_PointsMat)(res_Point[0],2);
       Point3 temp = SphereCenterP - tmpPonter;
       temp = temp-SphereCenterP;
       return SphereCenterP + temp*( 1.0 + temp.norm()/posPolarPoint[0] );
    }
}
bool CSphere::checkPointIs(const Point3& PolartempP,const Point3& posPolarPoint)
{
    double phi,thea;

    double PI = CATAZJUT::constants::Pi;

    double deta_radian = PI*30.0/360.0;
//    r=posPolarPoint[0];
     phi = posPolarPoint[1];
    thea = posPolarPoint[2];

    bool res1=false,res2=false;

    if( (phi+deta_radian)<=2*PI && (phi-deta_radian) >=0 )
    {
        if(std::fabs(phi-PolartempP[1])<deta_radian)
            res1=true;
    }else if((phi+deta_radian)>2*PI){
        if(PolartempP[1]>phi || PolartempP[1] < phi+deta_radian-2*PI )
            res1=true;
    }else{//(phi-deta_radian)<0
        if(PolartempP[1]>(phi-deta_radian + 2*PI) || PolartempP[1] < phi+deta_radian )
            res1=true;
    }

    if((thea+deta_radian)<=PI && (thea-deta_radian)>= 0 )
    {
        if(std::fabs(thea-PolartempP[2])<deta_radian)
            res2=true;
    }else if((thea+deta_radian)>PI){
        double temp = thea-deta_radian < PI-(deta_radian- (PI-thea))? (thea-deta_radian) : PI-(deta_radian- (PI-thea));
        if(PolartempP[2]>temp && PolartempP[2]<PI)
            res2=true;
    }else{//(thea-deta_radian)<0
        double temp = (thea+deta_radian)>deta_radian-thea? thea+deta_radian:deta_radian-thea;
        if(PolartempP[2]>0 && PolartempP[2]<temp )
            res2=true;
    }
    return res1 && res2;
}
void CSphere::output()
{
    Log::Output<<"****Sphere model equation of current cluster****\n";
    Log::Output<<">>( x - "<<(*m_Equation)[0]<<" )**2 + ( y - " <<(*m_Equation)[1]<<" )**2 + ( z - " ;
    Log::Output<<(*m_Equation)[2]<< " )**2 = " << (*m_Equation)[3]*(*m_Equation)[3] << std::endl;
    Log::Output<<"****End model equation****\n";
}




}
