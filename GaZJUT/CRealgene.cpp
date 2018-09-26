#include "CRealgene.h"
namespace GAZJUT{

CRealgene::CRealgene():CGenebase()
{

}
CRealgene::CRealgene(GENEVAR myVal):CGenebase(myVal)
{
    this->init(myVal);
}

CRealgene::~CRealgene()
{
    //dtor
}
double CRealgene::decode()
{
    return this->m_value;
}
void CRealgene::init(GENEVAR myVal)
{
    double low,high;
    assert(this->m_geneVar);

    if(this->m_geneVar==nullptr)
    {
        this->m_geneVar=new GENEVAR();
        *m_geneVar=myVal;
    }
    CRandomgenerator *rndgenerator=new CRandomgenerator();
    low =this->m_geneVar->min;
    high=this->m_geneVar->max;
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
