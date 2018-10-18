#include "CIOBase.h"

namespace CALCZJUT{

CIOBase::CIOBase(CATAZJUT::CPeriodicFramework* mbf)
:m_pPeriodicFramework(mbf)
{
    //ctor
}

CIOBase::~CIOBase()
{
    //dtor
}
Bitset CIOBase::input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE)
{
    Bitset a;
    return a;
}
void CIOBase::setConfiguration(CATAZJUT::CPeriodicFramework* mbf)
{
    this->m_pPeriodicFramework = mbf;
}

}


