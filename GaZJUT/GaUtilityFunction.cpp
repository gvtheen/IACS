#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<ctime>
#include "GaUtilityFunction.h"


namespace GAZJUT{

//
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

}  //namespace

