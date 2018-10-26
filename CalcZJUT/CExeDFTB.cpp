#include "CExeDFTB.h"

namespace CALCZJUT{

CExeDFTB::CExeDFTB(CParameter* mpara)
:CExeFitnessInterface(mpara)
{
    //ctor
}

CExeDFTB::~CExeDFTB()
{
    //dtor
}
CExeFitnessInterface* CExeDFTB::clone()
{

}
void CExeDFTB::init()
{

}
double CExeDFTB::CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist)
{

}
void CExeDFTB::ConvOrigToRawScore(std::vector<double>&)
{

}
void CExeDFTB::CheckInputFile()
{

}
bool CExeDFTB::IsNormalComplete()
{

}
double CExeDFTB::readFinalEnergy()
{

}
void CExeDFTB::getRelaxedGeometryCoord()
{

}


char* CExeDFTB::ExeName()
{
    return "DFTB";
}
}
