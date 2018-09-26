#include "CPeriodicFramework.h"
#include "CCalcCluster.h"
#include "CConfigurationBase.h"


namespace CALCZJUT{

CCalcCluster::CCalcCluster(CParameter* mPara)
:CCalcModeStruct(mPara)
{
    CConfigurationBase* temp = new CConfigurationBase(mPara);

    this->m_pPeriodicFramework = new CATAZJUT::CPeriodicFramework(*temp);
     // In this situation, this->m_pPeriodicFramework->m_pUnitCell is empty pointer NULL;
    delete temp;

}

CCalcCluster::~CCalcCluster()
{
    //dtor
}
void CCalcCluster::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{

}
std::vector<double>*  CCalcCluster::getGeneValuefromStruct()const
{

}
std::vector<GENEVAR>* CCalcCluster::GeneVarRange()
{

}
void CCalcCluster::Initialization(std::string& mth)
{
       // initialize from chemical
       // chemical formula = C10 H10 O12

    std::vector<std::pair<std::string,size_t>>* chemicalFormula = new (std::vector<std::pair<std::string,size_t>>);
}
void CCalcCluster::Initialization(char* mth)
{      // initialize from chemical formula
     std::string tmp(mth);
     this->Initialization(tmp);
}
void CCalcCluster::Initialization(CConfigurationBase& mth)  // initialize from exit
{

}




}
