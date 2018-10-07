#include "Cios.h"

namespace CALCZJUT{

Cios::Cios(CATAZJUT::CPeriodicFramework* mbf)
:m_pPeriodicFramework(mbf)
{
    //ctor
}

Cios::~Cios()
{
    //dtor
}
Bitset Cios::input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE)
{
    Bitset a;
    return a;
}
void Cios::setConfiguration(CATAZJUT::CPeriodicFramework* mbf)
{
    this->m_pPeriodicFramework = mbf;
}

}


