#include <cstddef>
#include <boost/config.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/algorithm/string.hpp>
#include "Point-Vector.h"
#include "utilFunction.h"
#include "../CataZJUT/Geometry.h"
#include "../GaZJUT/GaDeclaration.h"

using GAZJUT::GeneVAR;
using util::Bitset;

namespace util{
//
Vector4 SphereEquationFromPoints(const std::vector<Point3>& coordinate)
{
   double ans  =1e12,eps=1e-12,R,step,res,temp;

   Point3 P_pos;
   P_pos<<0,0,0;
   size_t pos;
   step=100;
   while(step>eps)
   {
       //pos = maxDist(P_pos
       res=0.0;
       for(size_t i=0;i<coordinate.size();i++)
       {
           temp=CATAZJUT::Geometry::distance(P_pos,coordinate[i]);
           if(res>temp)
           {
              pos=i;
              res= temp;
           }
       }
       R = CATAZJUT::Geometry::distance(P_pos,coordinate[pos]);
       if(step<ans)
          ans=step;
       P_pos = P_pos + step*(P_pos - coordinate[pos])/R;
       step=step*0.98;
   }
   Vector4 equation;
   equation<<P_pos[0],P_pos[1],P_pos[2],R;
   return equation;
}
double binaryDecode(const Bitset & myCode,GeneVAR myGenVar)
{
     int i,m;
     double low,high,sum;

     low= myGenVar.min;
     high=myGenVar.max;
     sum=0.0;
     m=myCode.size();
     for(i=0;i<m;i++)
          sum=sum+(double)(myCode[i])*std::pow(2.0,m-i-1);
     return  low+(high-low)*sum/(std::pow(2.0,m)-1);
}
int calcBitNum(GeneVAR myGeneVAR)
{
    double c,m;
    m=(myGeneVAR.max - myGeneVAR.min)/myGeneVAR.accuracy;
    c=std::log10(m)/std::log10(2.0);
    return (int)c+1;
}
void grayTobit(Bitset& data)
{
     int num=data.size();
     Bitset temp= data;
     data[0]=temp[0];
     for(int i=1;i<num;i++)
        data[i]=(data[i-1])^(temp[i]);
     temp.clear();
}
void bitTogray(Bitset& data)
{
     int num=data.size();
     Bitset temp= data;

     data[0]=temp[0];
     for(int i=1;i<num;i++)
        data[i]=(temp[i])^(temp[i-1]);
     temp.clear();
}
bool strcasecmp(const std::string& s1, const std::string& s2)
{
    std::string t1=s1;
    std::string t2=s2;
    boost::to_upper(t1);
    boost::to_upper(t2);
    return t1==t2;
}

}
