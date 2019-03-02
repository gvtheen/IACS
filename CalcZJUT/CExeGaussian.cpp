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
#include <boost/filesystem/fstream.hpp>
#include "../GaZJUT/GaUtilityFunction.h"
#include "CExeVASP.h"
#include "CModel2DSupport.h"
#include "CIOPoscar.h"
#include "CModelCluster.h"
#include "CModelClusterSupport.h"
#include "CExeGaussian.h"
#include "../CataZJUT/CFragment.h"
#include "../CataZJUT/CElement.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../Util/Point-Vector.h"
#include "../GaZJUT/CGaparameter.h"
#include "../Util/log.hpp"
#include "../CataZJUT/CBondTolerance.h"

using util::Log;
using util::Vector3;

namespace CALCZJUT{

//
CExeGaussian::CExeGaussian(CParameter* mpara)
:CExeFitnessInterface(mpara)
{
    m_pInputFile =new std::string("");
    *m_pInputFile=m_Parameter->sysName+".com";
    m_pParaFileAbsPath= new std::string("");
}
//
CExeGaussian::~CExeGaussian()
{
    if(m_pInputFile!=nullptr)
      delete m_pInputFile;

    if(m_pParaFileAbsPath!=nullptr)
      delete m_pParaFileAbsPath;
}
CExeFitnessInterface* CExeGaussian::clone()
{
    CExeGaussian* res =new CExeGaussian(this->m_Parameter);
    res->setInputFile(this->inputFile());

    return res;
}
void CExeGaussian::init()
{
    if(m_Parameter->output_struct_format=="")
       m_Parameter->output_struct_format="gjf";   //default value;

    std::string old_path_str=m_Parameter->parameterPath();
    old_path_str = old_path_str+"/"+ m_Parameter->sysName;
    old_path_str=old_path_str+".com";

    boost::filesystem::path tempPath(old_path_str);
    if(boost::filesystem::is_regular_file(tempPath))
        m_pParaFileAbsPath=new std::string(tempPath.string());
    else{
        tempPath.replace_extension(".gjf");
        if(boost::filesystem::is_regular_file(tempPath))
           m_pParaFileAbsPath= new std::string(tempPath.string());
        else{
           Log::Error<<tempPath.string()<<" file isnot exist! CExeGaussian::init()\n";
           boost::throw_exception(std::runtime_error("Parameter file isnot exist! CExeGaussian::init()."));//ERROR TREATMENT;
        }
    }
}
double CExeGaussian::CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalExist)
{
      //pid_t pid;
     size_t pid;   // replaced by last definition
     double res;
     //
     size_t currGeneration =m_Parameter->GaParameter()->Curr_Generation;

     Log::Info<<" Run Gaussian calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";
     //construct new object of CConfigurationBase class
     //transfer gene value to POSCAR file
     m_pCalcModeStruct->setGeneValueToStruct(RealValueOfGenome);
     m_pIO->setConfiguration(m_pCalcModeStruct->m_pPeriodicFramework);

     // change working path
     this->m_Parameter->setCurrentWorkPathAt(this->m_Parameter->currentGenerationNum(),pop_index);
     //copy parameter file to current working path

     std::string new_path_str=m_Parameter->currentWorkPath();

     this->m_Parameter->copyFileToPath(*m_pParaFileAbsPath,new_path_str);

     new_path_str=new_path_str+"/"+*m_pInputFile+".com";
     m_pIO->output( new_path_str ); //.mol file
     //
     //Check whether input files is OK?
     CheckInputFile();

     //pid=fork();
     if (pid<0){
         Log::Error<<"Building new process is error! CalcuRawFit_CExeDMol\n";
         boost::throw_exception(std::runtime_error("Building new process is error! Check the file: Error_information.txt."));//ERROR TREATMENT;
     }else if(pid==0){
         //run VASP program
     }else
        ;// wait(NULL);
     //
     //read CONTCAR file;
     this->m_Parameter->setCurrentWorkPathAt(this->m_Parameter->currentGenerationNum(),pop_index);
     getRelaxedGeometryCoord();
     std::string out_filename("gaussian_");

     if(IsNormalComplete()==true){
        isNormalExist=true;
        res=readFinalEnergy();
        out_filename = out_filename + std::to_string(currGeneration) + \
                   "_" + std::to_string(pop_index);
        Log::Info<<" Normally Finish Gaussian calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";
     }else{
        isNormalExist=false;
        res=9999999;
        out_filename = out_filename + std::to_string(currGeneration) + \
                   "_" + std::to_string(pop_index) + "_ERROR";
        Log::Info<<"InNormally Finish Gaussian calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";
     }

     CIOBase* tempIO=nullptr;
     this->getIO(m_Parameter->output_struct_format,m_pCalcModeStruct->periodicFramework(),tempIO);
     out_filename = out_filename + m_Parameter->output_struct_format;
     tempIO->output(out_filename);
     delete tempIO;

     return res;
}
void CExeGaussian::ConvOrigToRawScore(std::vector<double>& temporgValue)
{
    std::vector<double> tmpValue;
    tmpValue.assign(temporgValue.begin(),temporgValue.end());
    std::vector<double> ::iterator maxEnergy = std::max_element(tmpValue.begin(),tmpValue.end(),[](double a,double b){return a < b;});
    for(size_t i=0;i<tmpValue.size();i++)
         temporgValue[i] = *maxEnergy - tmpValue[i];
}
std::string CExeGaussian::ExeName()
{
    return "Gaussian";
}
void CExeGaussian::CheckInputFile()
{
   if(access(m_pInputFile->c_str(),F_OK) != 0)
    {
         std::string error_info = *(m_pInputFile) + std::string(" file is not exist!");
         Log::Error<<error_info<<" CheckInputFile_CExeGaussian\n";
         boost::throw_exception(std::runtime_error("Input files in Gaussian isnot exist! Check the file."));
    }
}
bool CExeGaussian::IsNormalComplete()
{
      std::string NormalLabel("Normal termination of Gaussian");

      bool res=false;
      std::vector<std::string> strVec;
      boost::algorithm::split(strVec,*m_pInputFile,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
      std::string file = strVec[0] + ".log";

      if(access(file.c_str(),F_OK) != 0 )
      {
           Log::Error<<file<<" file is no exist! CExeGaussian::IsNormalComplete\n";
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
          Log::Error<<e.what()<<" CExeGaussian::IsNormalComplete\n";
          return res;
      }
      in->close();
      return res;
}
double CExeGaussian::readFinalEnergy()
{
      std::string NormalLabel("Normal termination of Gaussian");
      std::vector<std::string> strVec;
      boost::algorithm::split(strVec,*m_pInputFile,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
      std::string file = strVec[0] + ".log";

      std::vector<double> energy_Vect;
      std::ifstream *in;
      try{
          in= new std::ifstream(file.c_str(),std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );
          std::string str;
          std::vector<std::string> strVec;
          while(!in->eof())
          {
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str, "SCF Done")>0)
             {
                  boost::algorithm::split(strVec,str,boost::algorithm::is_any_of("="),boost::algorithm::token_compress_on);
                  boost::algorithm::trim(strVec[1]);
                  boost::algorithm::split(strVec,strVec[1],boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
                  energy_Vect.push_back(std::stod(strVec[0]));
             }
             if(boost::algorithm::contains(str, NormalLabel)>0)
                  break;
          }
      }catch(const std::ifstream::failure& e){
          Log::Error<<e.what()<<"  see CExeGaussian::readFinalEnergy\n";
      }
      in->close();
      double temEnergy = energy_Vect[energy_Vect.size()-1];
      energy_Vect.clear();
      return temEnergy;
}
void CExeGaussian::getRelaxedGeometryCoord()
{
      std::vector<std::string> strVec;
      boost::algorithm::split(strVec,*m_pInputFile,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
      std::string file = strVec[0] + ".log";

      std::ifstream *in;
      Vector3 poptionVect;
      CATAZJUT::CElement* tempEle=nullptr;
      try{
          in= new std::ifstream(file.c_str(),std::ios::ate);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );
          std::string str;
          std::vector<std::string> strVec;
          while(!in->eof()){
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str, "Optimization completed")>0)
                  break;
          }
          while(!in->eof()){
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str, "Standard orientation")>0){
                 std::getline(*in,str,'\n');
                 std::getline(*in,str,'\n');
                 std::getline(*in,str,'\n');
                 std::getline(*in,str,'\n');
                 m_pCalcModeStruct->m_pPeriodicFramework->clear();
                 while(!in->eof()){
                     std::getline(*in,str,'\n');
                     boost::algorithm::trim(str);
                     boost::algorithm::split(strVec,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
                     poptionVect<<std::stod(strVec[3]),std::stod(strVec[4]),std::stod(strVec[5]);
                     tempEle = new CATAZJUT::CElement(std::stoi(strVec[1]));
                     m_pCalcModeStruct->m_pPeriodicFramework->addAtom(tempEle->symbol(),poptionVect);
                     delete tempEle;
                 }
             }
          }
      }catch(const std::ifstream::failure& e){
          Log::Error<<e.what()<<"  see CExeGaussian::readFinalEnergy\n";
          boost::throw_exception(std::runtime_error("ERROR in reading the file! Check the file:CExeGaussian::readFinalEnergy"));
      }
      in->close();
}
std::string& CExeGaussian::inputFile()const
{
   return *(this->m_pInputFile);
}
void CExeGaussian::setInputFile(const std::string& filename)
{
    if(this->m_pInputFile!=nullptr)
        delete this->m_pInputFile;

    this->m_pInputFile = new std::string(filename);
}



}
