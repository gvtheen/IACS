#ifndef CBONDPRIVATE_H
#define CBONDPRIVATE_H
#include <string>
#include<iostream>
namespace CATAZJUT{

class CBondPrivate
{
     public:
        CBondPrivate();
        CBondPrivate(std::string,std::string,double,double);
        virtual ~CBondPrivate();
     public:
       std::string m_1stAtom;
       std::string m_2ndAtom;
       double m_MinBondLength;
       double m_MaxBondLength;

};


}
#endif
