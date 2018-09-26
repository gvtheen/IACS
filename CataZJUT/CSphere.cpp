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
CSphere::CSphere(CConfigurationBase* configure)
{
    this->m_pConfiguration=configure;
}
CSphere::CSphere(CSphere& mht)
{
    this->m_pConfiguration=mht.configuration();
    this->m_Equation = mht.Equation();
}
CSphere::~CSphere()
{
    //dtor
}
void CSphere::CreateSphere()
{
   double ans  =1e12,eps=1e-12,R,step;
   CCartesianCoordinates* currentCart=this->m_pConfiguration->coordinates();
   Point3 P_pos;
   P_pos<<0,0,0;
   size_t pos;
   step=100;
   while(step>eps)
   {
       pos = maxDist(P_pos);
       R = CATAZJUT::Geometry::distance(P_pos,currentCart->position(pos));
       if(step<ans)
          ans=step;
       P_pos = P_pos + step*(P_pos - currentCart->position(pos))/R;
       step=step*0.98;
   }
   m_Equation<<P_pos[0],P_pos[1],P_pos[2],R;
}
size_t CSphere::maxDist(Point3 p1)
{
    double res=0,temp;
    size_t res_index;
    CCartesianCoordinates* currentCart=this->m_pConfiguration->coordinates();
    for(size_t i=0;i<currentCart->size();i++)
    {
        temp=CATAZJUT::Geometry::distance(p1,currentCart->position(i));
        if(res>temp)
        {
            res_index=i;
            res= temp;
        }
    }
    return res_index;
}
Eigen::Vector4d CSphere::Equation()const
{
    return this->m_Equation;
}
void CSphere::setConfiguration(CConfigurationBase* configure)
{
    this->m_pConfiguration=configure;
}
CConfigurationBase* CSphere::configuration()const
{
    return this->m_pConfiguration;
}
double CSphere::Radius() const
{
    return m_Equation[3];
}
Point3 CSphere::SphereCenter() const
{
    Point3 temp;
    temp<<m_Equation[0],m_Equation[1],m_Equation[2];
    return temp;
}
Point3 CSphere::toPolarCoordinate( const Point3& CartCoord)  // res= r,phi,thea;  unit:  Radians
{
    Point3 newCartCoord = CartCoord - SphereCenter();
    double r = newCartCoord.norm();
    Point3 res;
    double x=newCartCoord[0],y=newCartCoord[1],z=newCartCoord[2];
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
      res[1]= (360-CATAZJUT::constants::RadiansToDegrees*std::acos(x/std::sqrt(x*x+y*y)))*\
              CATAZJUT::constants::DegreesToRadians;
    res[2]=CATAZJUT::constants::RadiansToDegrees*std::acos(z/r);
    return res;
}
/*
  gene polarcoordinate on the point ( r,phi,thea  ) with the range of deta_radian;
*/
Point3 CSphere::geneCartesianCoordAt(const Point3& posPolarPoint)
{
    double maxDist=0;
    std::vector<size_t> res_Point;
    Point3 CurrentpolarPoint;
    CCartesianCoordinates* currentCart=this->m_pConfiguration->coordinates();
    for(size_t i=0;i<currentCart->size();i++)
    {
        CurrentpolarPoint = toPolarCoordinate(currentCart->position(i));
        if(checkPointIs(CurrentpolarPoint,posPolarPoint)==true)
        {
            res_Point.push_back(i);
            if(CurrentpolarPoint[0]>maxDist)
                maxDist=CurrentpolarPoint[0];
        }
    }
    std::vector<size_t>::iterator it;
    double deta_dist = 0.5;
    for(it=res_Point.begin();it!=res_Point.end();)
    {
        CurrentpolarPoint = toPolarCoordinate(currentCart->position(*it));
        if(CurrentpolarPoint[0]< maxDist -deta_dist )
          it =  res_Point.erase(it);
        else
          it++;
    }
    if(res_Point.size()>=3){
       Eigen::MatrixXd matPoint(res_Point.size(),3);
       for(size_t i=0;i<res_Point.size();i++)
           matPoint.row(i)=currentCart->position(res_Point[i]);
       CPlane *newPlane = new CPlane(matPoint);
       Point3 SphereCenterP=SphereCenter();
       Point3 projectPointOnPlane = newPlane->Projection(SphereCenterP);
       Point3 temp=projectPointOnPlane-SphereCenterP;
       delete newPlane;
       return SphereCenterP + temp*temp.norm()/posPolarPoint[0];
    }else if (res_Point.size()==2){
       Point3 SphereCenterP=SphereCenter();
       Point3 temp=CATAZJUT::Geometry::midpoint(currentCart->position(res_Point[0]),currentCart->position(res_Point[1]));
       temp=temp-SphereCenterP;
       return SphereCenterP + temp*temp.norm()/posPolarPoint[0];
    }else{
       Point3 SphereCenterP=SphereCenter();
       Point3 temp = SphereCenterP - currentCart->position(res_Point[0]);
       temp = temp-SphereCenterP;
       return SphereCenterP + temp*temp.norm()/posPolarPoint[0];
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





}
