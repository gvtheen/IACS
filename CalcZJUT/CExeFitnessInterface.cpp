#include <boost/algorithm/string.hpp>
#include "CExeFitnessInterface.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "CParameter.h"
#include "../CataZJUT/CFragment.h"
#include "CModelCluster.h"
#include "CModelClusterSupport.h"
#include "CModel2DSupport.h"
#include "CIOBase.h"
#include "CIOCar.h"
#include "CIOMol.h"
#include "CIOPoscar.h"
#include "CIOGjf.h"
#include "CIOCif.h"
#include "CIOCellFile.h"
#include "../Util/Point-Vector.h"
#include "../GaZJUT/CGaparameter.h"
#include "../Util/log.hpp"
#include "../CataZJUT/CBondTolerance.h"

using util::Log;
using util::Bitset;

namespace CALCZJUT{

CExeFitnessInterface::CExeFitnessInterface(CParameter* mpara)
:m_Parameter(mpara)
{

}

CExeFitnessInterface::~CExeFitnessInterface()
{
}

void CExeFitnessInterface::init()
{

}

double CExeFitnessInterface::CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist)
{
    return 0;
}

void CExeFitnessInterface::ConvOrigToRawScore(std::vector<double>& othr)
{

}
void CExeFitnessInterface::setCalcModeStruct(CModelBase* Temp_calcModeStruct)
{
    this->m_pCalcModeStruct = Temp_calcModeStruct;
}
CModelBase* CExeFitnessInterface::calcModeStruct()
{
    return this->m_pCalcModeStruct;
}
void CExeFitnessInterface::setIO(CIOBase* m_IO)
{
    this->m_pIO=m_IO;
}
CIOBase* CExeFitnessInterface::IO()const
{
    return this->m_pIO;
}
CIOBase* getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework)
{
    std::vector<std::string> vectstr;
    boost::algorithm::split(vectstr,file_name,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
    boost::algorithm::trim(vectstr[1]);
    if(vectstr[1]=="mol")
        return new CIOMol(currentPeriodicFramework);
    else if(vectstr[1]=="car")
        return new CIOCar(currentPeriodicFramework);
    else if(vectstr[1]=="poscar")
        return new CIOPoscar(currentPeriodicFramework);
    else if(vectstr[1]=="gjf")
        return new CIOGjf(currentPeriodicFramework);
    else if(vectstr[1]=="cif")
        return new CIOCif(currentPeriodicFramework);
    else if(vectstr[1]=="cell")
        return new CIOCellFile(currentPeriodicFramework);
    else
        ;

    return nullptr;
}

}
