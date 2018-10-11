#include "CCalcStructBasePool.h"



namespace CALCZJUT{

CCalcStructureBasePool::CCalcStructureBasePool(CParameter* othr)
:m_Parameter(othr)
{
    //ctor
}

CCalcStructureBasePool::~CCalcStructureBasePool()
{
    //dtor
}
CCalcModeStruct* CCalcStructureBasePool::operator[](size index)
{
    if(index>=this->m_CalcStructPool.size())
    {
         Log::Error<<"Index is out of range! CCalcStructureBasePool::operator[]!\n";
         boost::throw_exception(std::runtime_error("Index is out of range! CCalcStructureBasePool::operator[]!"));
    }
    return this->m_CalcStructPool[index];
}


}
