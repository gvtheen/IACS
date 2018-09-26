#include "CGenebase.h"
#include "GaUtilityFunction.h"

namespace GAZJUT{
CGenebase::CGenebase()
{
    this->m_geneVar=nullptr;
}
CGenebase::CGenebase(GENEVAR var)
{
    this->init(var);
}
double CGenebase::decode()
{
    return this->m_value;
}
void CGenebase::init(GENEVAR var)
{
    this->m_geneVar=new GENEVAR();
    *m_geneVar=var;
}
void CGenebase::updatecode(double m)
{

}
std::vector<unsigned int>* CGenebase::bitGene()
{
    return nullptr;
}
double CGenebase::realGene()
{
    return this->m_value;
}
double CGenebase::value()
{
    return this->m_value;
}
int CGenebase::bitNum()
{
    return 0;
}
CGenebase::~CGenebase()
{
    delete this->m_geneVar;
}
}
