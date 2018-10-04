#include "CRealgene.h"
namespace GAZJUT{

CRealgene::CRealgene():CGenebase()
{
    this->codeType=REAL;
}
CRealgene::CRealgene(GeneVAR myVal):CGenebase(myVal)
{
    this->init(myVal);
    this->codeType=REAL;
}

CRealgene::~CRealgene()
{
    //dtor
}
double CRealgene::decode()
{
    return this->m_value;
}
void CRealgene::init(GeneVAR myVal)
{
    double low,high;
    assert(this->m_GeneVAR);

    if(this->m_GeneVAR==nullptr)
    {
        this->m_GeneVAR=new GeneVAR();
        *m_GeneVAR=myVal;
    }
    CRandomgenerator *rndgenerator=new CRandomgenerator();
    low =this->m_GeneVAR->min;
    high=this->m_GeneVAR->max;
    this->m_value=low + (high-low)*rndgenerator->uniformRandom01(100);
}
void CRealgene::updatecode(double value)
{
    this->m_value=value;
}
double& CRealgene::realGene()
{
    return this->m_value;
}
int CRealgene::bitNum()
{
    return 1;
}

}
