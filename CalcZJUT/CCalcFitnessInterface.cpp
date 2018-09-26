#include "CCalcFitnessInterface.h"
#include "CPeriodicFramework.h"
#include "CParameter.h"
#include "CGaparameter.h"
#include "CCalcCluster.h"
#include "CCalcClusterSupport.h"
#include "CCalc2DSupport.h"
#include <boost/algorithm/string.hpp>

namespace CALCZJUT{

CCalcFitnessInterface::CCalcFitnessInterface(CParameter* mpara)
:m_Parameter(mpara)
{
    pop_run_state.resize(m_Parameter->GaParameter()->PopNum(),false);
}

CCalcFitnessInterface::~CCalcFitnessInterface()
{
}
void CCalcFitnessInterface::setCalcModeStruct(CCalcModeStruct* Temp_calcModeStruct)
{
    this->m_pCalcModeStruct = Temp_calcModeStruct;
}
CCalcModeStruct* CCalcFitnessInterface::calcModeStruct()
{
    return this->m_pCalcModeStruct;
}
void CCalcFitnessInterface::setIO(Cios* m_IO)
{
    this->m_pIO=m_IO;
}
Cios* CCalcFitnessInterface::IO()const
{
    return this->m_pIO;
}
Cios* getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework)
{
    std::vector<std::string> vectstr;
    boost::algorithm::split(vectstr,file_name,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
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

    return nullptr;
}

}
