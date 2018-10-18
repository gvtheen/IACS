#include "CCalcStructBasePool.h"



namespace CALCZJUT{

CStructPoolBase::CStructPoolBase(CParameter* othr)
:m_Parameter(othr)
{
    //ctor
}

CStructPoolBase::~CStructPoolBase()
{
    //dtor
}
CModelBase* CStructPoolBase::operator[](size index)
{
    if(index>=this->m_CalcStructPool.size())
    {
         Log::Error<<"Index is out of range! CStructPoolBase::operator[]!\n";
         boost::throw_exception(std::runtime_error("Index is out of range! CStructPoolBase::operator[]!"));
    }
    return this->m_CalcStructPool[index];
}
void CStructPoolBase::GeneVARRange(std::vector<GeneVAR>& mth)
{

}
void CStructPoolBase::init()
{

}




}
