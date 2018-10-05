#include "CGenebase.h"
#include "GaUtilityFunction.h"

namespace GAZJUT{
CGenebase::CGenebase()
{
    this->m_GeneVAR=nullptr;
}
CGenebase::CGenebase(GeneVAR* var)
{
    this->init(var);
}
double CGenebase::decode()
{
    return this->m_value;
}
void CGenebase::init(GeneVAR* var)
{
    m_GeneVAR=var;
}
void CGenebase::updatecode(double m)
{

}
Bitset& CGenebase::bitGene()
{
    Bitset a1;        // nothing
    return a1;
}
double CGenebase::realGene()
{
    return this->m_value;
}
double CGenebase::value()
{
    return this->m_value;
}
size_t CGenebase::bitNum()
{
    return 0;
}
CGenebase::~CGenebase()
{
    //delete this->m_GeneVAR;
}


}
