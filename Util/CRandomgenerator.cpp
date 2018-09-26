#include <boost/random.hpp>
#include<time.h>
#include <limits>
#include "CRandomgenerator.h"

namespace util{

unsigned int CRandomgenerator::classcalled=0;

CRandomgenerator::CRandomgenerator()
{
    classcalled++;
    //cout<<"new Ranm: " <<classcalled<<endl;
    m_seedInt=classcalled;
    if (classcalled==std::numeric_limits<unsigned int>::max())
        classcalled=0;
}
CRandomgenerator::~CRandomgenerator()
{
    //dtor
}
double CRandomgenerator::uniformRandom01(int seedNum)
{
    boost::mt19937 rng(time(0)+seedNum*seedNum+m_seedInt);
    boost::uniform_01<boost::mt19937&> u01(rng);
    boost::uniform_int<> ui(0,MNNN);
    double *myRandom=new double[MNNN];
    for(int i=0;i<MNNN;i++)
        myRandom[i]=u01();
    double res= myRandom[ui(rng)];
    delete myRandom;
    return res;
}
int CRandomgenerator::uniformRandomRandge(int start, int stop,int seedNum)
{
     boost::mt19937 rng(time(0)+seedNum*seedNum*10+m_seedInt*m_seedInt*10);
     boost::uniform_int<> ui(start,stop);
     return ui(rng);
}


}
