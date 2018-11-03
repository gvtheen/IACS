#include "unistd.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "../GaZJUT/GaUtilityFunction.h"
#include "CExeCastep.h"
#include "CParameter.h"
#include "../Util/log.hpp"
#include "CIOCar.h"
#include "CModel2DSupport.h"
#include "CModelCluster.h"
#include "CModelClusterSupport.h"
#include "../CataZJUT/CFragment.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../Util/Point-Vector.h"
#include "../Util/log.hpp"
#include "../GaZJUT/CGaparameter.h"
#include "../CataZJUT/CBondTolerance.h"
using util::Log;

namespace CALCZJUT{

CExeCastep::CExeCastep(CParameter* mpara)
:CExeFitnessInterface(mpara)
{
    std::string sys_name=this->m_Parameter->sysName;
    m_pInputFile.push_back(sys_name+".cell");
    m_pInputFile.push_back(sys_name+".param");
}

CExeCastep::~CExeCastep()
{
    //dtor
}
void CExeCastep::init()
{
    if(m_Parameter->output_struct_format=="")
        m_Parameter->output_struct_format = "cell";
}
double CExeCastep::CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalExist)
{
    pid_t pid;
     double res;
     size_t currGeneration =m_Parameter->GaParameter()->Curr_Generation;

     Log::Info<<" Run CASTEP calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";

     //construct new object of CPeriodicFramework class
     //transfer gene value to structure file
     m_pCalcModeStruct->setGeneValueToStruct(RealValueOfGenome);
     //
     m_pIO->setConfiguration(m_pCalcModeStruct->m_pPeriodicFramework);
     //write structure in the car file
     m_pIO->output(m_pInputFile[0]); //.cell file
     //
     //Check whether input files is OK?
     CheckInputFile();

     //pid=fork();

     if (pid<0){
         Log::Error<<"Building new process is error! CExeCastep::CalcuRawFit"<<std::endl;
         boost::throw_exception(std::runtime_error("Building new process is error! Check the file: Error_information.txt."));//ERROR TREATMENT;
     }else if(pid==0){
         //run DMOL program
     }else
         ;//wait(NULL);
     //
     //read CONTCAR file;
     getRelaxedGeometryCoord();
     std::string out_filename("relaxed_");

     if(IsNormalComplete()==true){
        isNormalExist=true;
        res=readFinalEnergy();
        out_filename = out_filename + std::to_string(currGeneration) + \
                   "_" + std::to_string(pop_index);
        Log::Info<<" Normally Finish CASTEP calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";
     }else{
        isNormalExist=false;
        res=9999999;
        out_filename = out_filename + std::to_string(currGeneration) + \
                   "_" + std::to_string(pop_index) + "_ERROR";
        Log::Info<<"InNormally Finish CASTEP calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";
     }

     CIOBase* tempIO = this->getIO(m_Parameter->output_struct_format,m_pCalcModeStruct->periodicFramework());
     out_filename = out_filename + m_Parameter->output_struct_format;
     tempIO->output(out_filename);
     delete tempIO;

     return res;
}
void CExeCastep::ConvOrigToRawScore(std::vector<double>& temporgValue)
{
    std::vector<double> tmpValue;
    tmpValue.assign(temporgValue.begin(),temporgValue.end());
    std::vector<double> ::iterator maxEnergy = std::max_element(tmpValue.begin(),tmpValue.end(),[](double a,double b){return a < b;});
    for(size_t i=0;i<tmpValue.size();i++)
         temporgValue[i] = *maxEnergy - tmpValue[i];
}
char* CExeCastep::ExeName()
{
    return "CASTEP";
}
void CExeCastep::CheckInputFile()
{
    bool res=true;
    for(size_t i=0;i<this->m_pInputFile.size();i++)
        if(access(m_pInputFile[i].c_str(),F_OK) != 0)
        {
           std::string error_info = m_pInputFile[i] + std::string(" file is not exist!");
           Log::Error<<error_info<<" CExeCastep::CheckInputFile()\n";
           res=false;
        }
    if(res==false)
    {
         Log::Error<<"Input files in DMOL is error! Check CExeCastep::CheckInputFile()"<<std::endl;
         boost::throw_exception(std::runtime_error("Input files in DMOL is error! CExeCastep::CheckInputFile()"));
    }
}

}
