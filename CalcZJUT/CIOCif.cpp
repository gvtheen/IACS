#include <iostream>
#include "unistd.h"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <string>
#include "../Util/Point-Vector.h"
#include "CIOCif.h"

namespace CALCZJUT{

CIOCif::CIOCif(CATAZJUT::CPeriodicFramework* mpa)
:CIOBase(mpa)
{
    //ctor
}

CIOCif::~CIOCif()
{
    //dtor
}
void CIOCif::output(const std::string& file)
{

}
void CIOCif::input(std::string file)
{

}
Bitset CIOCif::input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE mode)
{

}

}
