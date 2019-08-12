/******************************************************************************
**
** Copyright (C) 2019-2031 Dr.Gui-lin Zhuang <glzhuang@zjut.edu.cn>
** All rights reserved.
**
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/
#include "unistd.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include "../GaZJUT/GaUtilityFunction.h"
#include "CExeVASP.h"
#include "CModel2DSupport.h"
#include "CIOPoscar.h"
#include "CModelCluster.h"
#include "CModelClusterSupport.h"
#include "../CataZJUT/CFragment.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../Util/Point-Vector.h"
#include "../GaZJUT/CGaparameter.h"
#include "../Util/log.hpp"
#include "../CataZJUT/CBondTolerance.h"

namespace CATAZJUT{
  class CFragment;
  class CConfigurationBase;
}

using util::Point3;
using util::Log;
namespace CALCZJUT{

CExeVASP::CExeVASP(CParameter* mpara):CExeFitnessInterface(mpara)
{
    m_pInputFile.push_back(new std::string("INCAR"));
    m_pInputFile.push_back(new std::string("KPOINTS"));
    m_pInputFile.push_back(new std::string("POTCAR"));
}

CExeVASP::~CExeVASP()
{
    for(size_t i=0;i<m_pInputFile.size();i++)
        delete m_pInputFile[i];
    m_pInputFile.clear();

}

CExeFitnessInterface* CExeVASP::clone()
{
     CExeVASP* res =new CExeVASP(this->m_Parameter);
     return res;
}
void CExeVASP::init()
{
    if(m_Parameter->output_struct_format=="")
       m_Parameter->output_struct_format="poscar";   //default value;
    // check whether other parameter files is exist.
    m_Parameter->checkExeNecessaryFiles(m_pInputFile);
    m_pParaFileAbsPath=new std::string(m_Parameter->parameterPath());
    this->m_pIO=new CIOPoscar(nullptr);

}
std::string CExeVASP::ExeName()
{
    return "VASP";
}
void CExeVASP::ConvOrigToRawScore(std::vector<double>& temporgValue)
{
    std::vector<double> tmpValue;
    tmpValue.assign(temporgValue.begin(),temporgValue.end());
    std::vector<double> ::iterator maxEnergy = std::max_element(tmpValue.begin(),tmpValue.end(),[](double a,double b){return a < b;});
    for(size_t i=0;i<tmpValue.size();i++)
         temporgValue[i] = *maxEnergy - tmpValue[i];
}
double CExeVASP::CalcuRawFit(std::vector<double>&RealValueOfGenome,size_t& pop_index, bool& isNormalExist)
{
      //pid_t pid;
     size_t pid;   // replaced by last definition
     double res;

     size_t currGeneration = m_Parameter->GaParameter()->Curr_Generation;

     //construct new object of CConfigurationBase class
     //transfer gene value to POSCAR file
     // if( currGeneration != 0 )
     m_pCalcModeStruct->setGeneValueToStruct(RealValueOfGenome);
     #ifdef DEBUG
       Log::Debug<<"*********** 1-CExeVASP::CalcuRawFit***********"<< std::endl;
     #endif // DEBU
     m_pIO->setConfiguration(m_pCalcModeStruct->periodicFramework());
     #ifdef DEBUG
       Log::Debug<<"*********** 2-CExeVASP::CalcuRawFit***********"<< std::endl;
     #endif // DEBU
     m_Parameter->setCurrentWorkPathAt(this->m_Parameter->currentGenerationNum()+1,pop_index+1);
     #ifdef DEBUG
       Log::Debug<<"*********** 3-CExeVASP::CalcuRawFit***********"<< std::endl;
     #endif // DEBU
     //copy parameter files into current path
     std::string new_path_str=m_Parameter->currentWorkPath();
     std::string file_path_str;
     for(size_t i=0;i<this->m_pInputFile.size();i++)
     {
         file_path_str=*m_pParaFileAbsPath + "/"+ *(m_pInputFile[i]);
         this->m_Parameter->copyFileToPath(file_path_str,new_path_str);
     }

     //output initial structure
     m_pIO->output("POSCAR");
     //
     //Check whether input files is OK?
     CheckInputFile();

    // pid=fork();
     if (pid<0){
         Log::Error<<"Building new process is error! CalcuRawFit_CExeVASP"<<std::endl;
         boost::throw_exception(std::runtime_error("Building new process is error! Check the file: Error_information.txt."));//ERROR TREATMENT;
     }else if(pid==0){
         //run VASP program
     }else
         ;//wait(NULL);
     //
     //read CONTCAR file;
     m_Parameter->setCurrentWorkPathAt(this->m_Parameter->currentGenerationNum()+1,pop_index+1);
     getRelaxedGeometryCoord();
     std::string out_filename("VASP_");

     if(this->IsNormalComplete()==true){
        isNormalExist=true;
        res=readFinalEnergy();
        out_filename = out_filename + std::to_string(currGeneration) + \
                   "_" + std::to_string(pop_index);
        Log::Info<<"Normally finish VASP of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";
     }else{
        isNormalExist=false;
        res=9999999;
        out_filename = out_filename + std::to_string(currGeneration) + \
                   "_" + std::to_string(pop_index);
        Log::Info<<"InNormally finish VASP calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";
     }

     // Use the output file format to output structure;
     CIOBase* tempIO=nullptr;
     this->getIO(m_Parameter->output_struct_format,m_pCalcModeStruct->periodicFramework(),tempIO);
     out_filename = out_filename + m_Parameter->output_struct_format;

     //change current path
     this->m_Parameter->setCurrentWorkPathAt(CParameter::SCRATCH);
     tempIO->output(out_filename);

     delete tempIO;
     //release IO memory!

     return res;
}
void CExeVASP::CheckInputFile()
{
    bool res=true;
    for(unsigned int i=0;i<this->m_pInputFile.size();i++)
        if(access(m_pInputFile[i]->c_str(),F_OK) != 0)
        {
           std::string error_info = *(m_pInputFile[i]) + std::string(" file is not exist!");
           Log::Error<<error_info <<" CheckInputFile_CExeVASP\n";
           res=false;
        }
    if(res==false)
    {
         Log::Error<<"Input files in VASP is error!CheckInputFile_CExeVASP"<<std::endl;
         boost::throw_exception(std::runtime_error("Input files in VASP is error! Check the file: Error_information.txt."));
    }
}
bool CExeVASP::IsNormalComplete()
{
      std::string NormalLabel("General timing and accounting informations");

      bool res=false;
      if(access("OUTCAR",F_OK) != 0 )
      {
           Log::Error<<"OUTCAR file is no exist! IsNormalComplete_CExeVASP\n";
           return res;
      }
      std::ifstream *in;
      try{
          in= new std::ifstream("OUTCAR",std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );
          std::string str;
          while(!in->eof())
          {
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str, NormalLabel)>0)
             {
                 res = true;
                 break;
             }
          }
      }catch(const std::ifstream::failure& e){
          Log::Error<<e.what()<<" IsNormalComplete_CExeVASP\n";
          exit(-1);
      }
      in->close();
      return res;
}
double CExeVASP::readFinalEnergy()
{
      std::string NormalLabel("General timing and accounting informations");
      std::vector<double> *energy_Vect = new (std::vector<double>);
      std::ifstream *in;
      try{
          in= new std::ifstream("OUTCAR",std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );
          std::string str;
          std::vector<std::string> strVec;
          while(!in->eof())
          {
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str, "free  energy   TOTEN")>0)
             {
                  boost::algorithm::split(strVec,str,boost::algorithm::is_any_of("="),boost::algorithm::token_compress_on);
                  boost::algorithm::trim(strVec[1]);
                  boost::algorithm::split(strVec,strVec[1],boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
                  energy_Vect->push_back(std::stod(strVec[0]));
             }
             if(boost::algorithm::contains(str, NormalLabel)>0)
                  break;
          }
      }catch(const std::ifstream::failure& e){
          Log::Error<<e.what()<<"readFinalEnergy_CExeVASP\n";
          exit(-1);
      }
      in->close();
      double FinalEnergy = energy_Vect->at(energy_Vect->size()-1);

      delete energy_Vect;

      return FinalEnergy;
}
void CExeVASP::getRelaxedGeometryCoord()
{
      m_pIO->m_pPeriodicFramework->clear();
      m_pIO->input("CONTCAR");
}



}
