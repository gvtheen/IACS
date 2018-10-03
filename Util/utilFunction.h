#ifndef GACATA_H
#define GACATA_H

#include <cstddef>
#include <boost/config.hpp>
#include <boost/preprocessor/stringize.hpp>
#include "Point-Vector.h"
#include "../CataZJUT/Geometry.h"
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



}



#endif
