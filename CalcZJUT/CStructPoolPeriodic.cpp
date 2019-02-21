#include <boost/algorithm/string.hpp>
#include "CStructPoolPeriodic.h"
#include "CModelPeriodicStruct.h"
#include "CParameter.h"
#include "../GaZJUT/CGaparameter.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../util/log.hpp"
#include "CIOBase.h"
#include "CIOMol.h"
#include "CIOCar.h"
#include "CIOGjf.h"
#include "CIOPoscar.h"

using util::Log;

namespace CALCZJUT{

CStructPoolPeriodic::CStructPoolPeriodic(CParameter* othr)
:CStructPoolBase(othr)
{
    for(size_t i=0;i<this->m_pParameter->GaParameter()->PopNum();i++){
        this->m_CalcStructPool.push_back(new CModelPeriodicStruct(this->m_pParameter,i));
        m_CalcStructPool[m_CalcStructPool.size()-1]->periodicFramework()->setExcludeBond(m_pParameter->excludeBond);
        m_CalcStructPool[m_CalcStructPool.size()-1]->periodicFramework()->setTolerancefactor(m_pParameter->bondToleranceFactor);
    }
}

CStructPoolPeriodic::~CStructPoolPeriodic()
{

}
void CStructPoolPeriodic::init()
{
    if(this->m_pParameter->periodic_structure_Input_File.size()!=0){
        this->Initialization(m_pParameter->cluster_Input_File);
    }else if(this->m_pParameter->periodic_structure_Formula!=""){
        this->Initialization(m_pParameter->cluster_Formula);
    }else{
       Log::Error<< " Chemical formula or initially structural files is required. CStructPoolPeriodic::init()!\n";
       boost::throw_exception(std::runtime_error("Chemical formula and structural files is required. CStructPoolPeriodic::init()!!\n"));//ERROR TREATMENT;
    }
}
void CStructPoolPeriodic::GeneVARRange(std::vector<GeneVAR>&  mht)
{
   double maxRadius=0.0;
   std::vector<GeneVAR> resTemp;
   for(size_t i=0;i<this->m_CalcStructPool.size();i++){
        this->m_CalcStructPool[i]->GeneVARRange(resTemp);
        if(maxRadius > resTemp[0].max)
            maxRadius = resTemp[0].max;
   }
   maxRadius = maxRadius + 1.0;
   size_t num=this->m_CalcStructPool[0]->m_pPeriodicFramework->size();

   for(size_t i=0;i<3*num;i++)
      mht.push_back({-1*maxRadius,maxRadius,0.001});
}
void CStructPoolPeriodic::Initialization(const std::string& mth)         // initialize from chemical formula
{
          // initialize from chemical
       // chemical formula = C10 H10 O12
/*
    e.g.:   mth = C10 N S10        // delimiter = blank " "
*/
    //dealwith chemical formula by blank
    std::string str=mth;
    std::vector<std::string> vecStr;
    boost::algorithm::trim(str);
    boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);

    std::pair<std::string,size_t> ChemicalStr;
    std::vector<std::pair<std::string,size_t>> tempChemFormula;

    for(size_t i=0;i<vecStr.size();i++){
        ChemicalStr.first  = boost::algorithm::trim_copy_if(vecStr[i], boost::algorithm::is_digit());
        //check whether the label is OK.
        if( CATAZJUT::CElement::isValidSymbol(ChemicalStr.first) == false ){
           Log::Error<<vecStr[i]<< " is error. Initialization_CModelCluster!\n";
           boost::throw_exception(std::runtime_error(vecStr[i] + " is error. Initialization_CModelCluster!\n"));//ERROR TREATMENT;
        }
        str = boost::algorithm::trim_copy_if(vecStr[i], boost::algorithm::is_alpha());
        if ( str=="" )
            ChemicalStr.second = 1;
        else
            ChemicalStr.second = std::stoi(str);
        tempChemFormula.push_back(ChemicalStr);
    }
    for(size_t i=0;i<this->m_pParameter->GaParameter()->PopNum();i++)
        RandomBuildFromChemicalFormula(this->m_CalcStructPool[i]->periodicFramework(),tempChemFormula);
    tempChemFormula.clear();
}
void CStructPoolPeriodic::Initialization(const char* mth)               // initialize from chemical formula
{
     std::string tmp(mth);
     this->Initialization(tmp);
}
void CStructPoolPeriodic::Initialization(const std::vector<std::string*>& inputfiles)     // initialize from exit structure
{
    std::string str;
    std::vector<std::string> vecStr;
    CIOBase* inputIO;
    size_t pos=0;
    for(size_t i=0;i<inputfiles.size();i++){
       str=*(inputfiles[i]);
       boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
       if(boost::iequals(vecStr[1],"gjf") || boost::iequals(vecStr[1],"mol") ){
          Log::Error<<str << " file is not supported by periodic model! CStructPoolPeriodic::Initialization\n";
          boost::throw_exception(std::runtime_error(str+ " file is not supported by periodic model! CStructPoolPeriodic::Initialization!"));
       }else if(boost::iequals(vecStr[1],"car"))
          inputIO = new CIOCar(this->m_CalcStructPool[pos]->periodicFramework());
       else if(boost::iequals(vecStr[1],"poscar")==0)
          inputIO = new CIOPoscar(this->m_CalcStructPool[pos]->periodicFramework());
       else{
          Log::Error<< str << " file is not supported by cluster model! Initialization_CModelCluster\n";
          boost::throw_exception(std::runtime_error( str + " file is not supported by cluster model! Initialization_CModelCluster!"));
       }
       // read structure and save m_PopuPeriodicFramework[pos++]
       inputIO->input(str);

       m_CalcStructPool[pos]->setChemicalFormula( m_CalcStructPool[pos]->periodicFramework()->composition() );
       pos++;
    }

    for(size_t i=pos;i<this->m_pParameter->GaParameter()->PopNum();i++){
      this->m_CalcStructPool[i]->setChemicalFormula(this->m_CalcStructPool[0]->chemicalFormula());
      RandomBuildFromChemicalFormula(this->m_CalcStructPool[i]->periodicFramework(),
                                     this->m_CalcStructPool[i]->chemicalFormula());
    }
}
void CStructPoolPeriodic::RandomBuildFromChemicalFormula(CATAZJUT::CPeriodicFramework* strut,
                                            std::vector<std::pair<std::string,size_t>>& chemFormula)
{

}



}
