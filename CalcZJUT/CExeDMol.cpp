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
#include <boost/algorithm/string.hpp>
#include "../GaZJUT/GaUtilityFunction.h"
#include "CExeDMol.h"
#include "CParameter.h"
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

CExeDMol::CExeDMol(CParameter* mpara)
:CExeFitnessInterface(mpara)
{
    dmol_inputfile =new std::string(m_Parameter->runCmd[m_Parameter->runCmd.size()-1]);

    m_pInputFile.push_back(new std::string(*dmol_inputfile+".car"));
    m_pInputFile.push_back(new std::string(*dmol_inputfile+".input"));
    m_pInputFile.push_back(new std::string(*dmol_inputfile+".mdf"));
}
CExeDMol::~CExeDMol()
{
    if(dmol_inputfile!=nullptr)
        delete dmol_inputfile;

    for(size_t i=0;i<m_pInputFile.size();i++)
        delete m_pInputFile[i];
    m_pInputFile.clear();
}
CExeFitnessInterface* CExeDMol::clone()
{
    CExeDMol* res =new CExeDMol(this->m_Parameter);
    res->setInputFile(this->inputFile());

    return res;
}
void CExeDMol::init()
{
     #ifdef DEBUG
       Log::Debug<<" Initialize CExeDMol object!" <<std::endl;
     #endif // DEBUG

    if(m_Parameter->output_struct_format=="")
       m_Parameter->output_struct_format="car";   //default value;

    m_Parameter->checkExeNecessaryFiles(m_pInputFile);
}
std::string& CExeDMol::inputFile()const
{
    return *(this->dmol_inputfile);
}
char* CExeDMol::ExeName()
{
    return "DMol3";
}
void CExeDMol::setInputFile(const std::string& filename)
{
     if(dmol_inputfile!=nullptr)
        delete dmol_inputfile;
     dmol_inputfile =new std::string(filename);
     m_pInputFile.push_back(new std::string(*dmol_inputfile+".car"));
     m_pInputFile.push_back(new std::string(*dmol_inputfile+".input"));
     m_pInputFile.push_back(new std::string(*dmol_inputfile+".mdf"));
}
double CExeDMol::CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalExist)
{
     pid_t pid;
     double res;
     size_t currGeneration =m_Parameter->GaParameter()->Curr_Generation;

     Log::Info<<" Run DMol calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";

     //construct new object of CPeriodicFramework class
     //transfer gene value to structure file
     m_pCalcModeStruct->setGeneValueToStruct(RealValueOfGenome);
     //change current path
     m_Parameter->setCurrentWorkPathAt(this->m_Parameter->currentGenerationNum(),pop_index);
     //set outputing object
     m_pIO->setConfiguration(m_pCalcModeStruct->m_pPeriodicFramework);
     //write structure in the car file
     m_pIO->output(*(m_pInputFile[0])+".car"); //.car file
     //
     //Check whether input files is OK?
     CheckInputFile();

     //pid=fork();

     if (pid<0){
         Log::Error<<"Building new process is error! CalcuRawFit_CExeDMol"<<std::endl;
         boost::throw_exception(std::runtime_error("Building new process is error! Check the file: Error_information.txt."));//ERROR TREATMENT;
     }else if(pid==0){
         //run DMOL program
     }else
         ;//wait(NULL);
     //
     //read CONTCAR file;
     m_Parameter->setCurrentWorkPathAt(this->m_Parameter->currentGenerationNum(),pop_index);
     getRelaxedGeometryCoord();
     std::string out_filename("dmol_");

     if(IsNormalComplete()==true){
        isNormalExist=true;
        res=readFinalEnergy();
        out_filename = out_filename + std::to_string(currGeneration) + \
                   "_" + std::to_string(pop_index);
        Log::Info<<" Normally Finish DMol calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";
     }else{
        isNormalExist=false;
        res=9999999;
        out_filename = out_filename + std::to_string(currGeneration) + \
                   "_" + std::to_string(pop_index) + "_ERROR";
        Log::Info<<"InNormally Finish DMol calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";
     }

     CIOBase* tempIO = this->getIO(m_Parameter->output_struct_format,m_pCalcModeStruct->periodicFramework());
     out_filename = out_filename + m_Parameter->output_struct_format;
     this->m_Parameter->setCurrentWorkPathAt(CParameter::SCRATCH);
     tempIO->output(out_filename);
     delete tempIO;

     return res;
}
void CExeDMol::ConvOrigToRawScore(std::vector<double>& temporgValue)
{
    std::vector<double> tmpValue;
    tmpValue.assign(temporgValue.begin(),temporgValue.end());
    std::vector<double> ::iterator maxEnergy = std::max_element(tmpValue.begin(),tmpValue.end(),[](double a,double b){return a < b;});
    for(size_t i=0;i<tmpValue.size();i++)
         temporgValue[i] = *maxEnergy - tmpValue[i];
}
void CExeDMol::CheckInputFile()
{
    bool res=true;
    for(size_t i=0;i<this->m_pInputFile.size();i++)
        if(access(m_pInputFile[i]->c_str(),F_OK) != 0)
        {
           std::string error_info = *(m_pInputFile[i]) + std::string(" file is not exist!");
           Log::Error<<error_info<<" CheckInputFile_CExeDMol\n";
           res=false;
        }
    if(res==false)
    {
         Log::Error<<"Input files in DMOL is error! CheckInputFile_CExeDMol"<<std::endl;
         boost::throw_exception(std::runtime_error("Input files in DMOL is error! Check the file."));
    }
}
bool CExeDMol::IsNormalComplete()
{
      std::string NormalLabel("DMol3 job finished successfully");
      std::string* filename=new std::string(*dmol_inputfile);
      *filename = (*filename) + ".outmol";
      bool res=false;
      if(access(filename->c_str(),F_OK) != 0 )
      {
           Log::Error<<(*filename) <<" file is no exist! IsNormalComplete_CExeDMol"<<std::endl;
           return res;
      }
      std::ifstream *in;
      try{
          in= new std::ifstream(filename->c_str(),std::ifstream::in);
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
          Log::Error<<e.what()<<" IsNormalComplete_CExeDMol!"<<std::endl;
          exit(-1);
      }
      in->close();
      delete filename;
      return res;
}
double CExeDMol::readFinalEnergy()
{
      std::string NormalLabel("DMol3 job finished successfully");
      std::vector<double> *energy_Vect = new (std::vector<double>);
      std::ifstream *in;
      std::string* filename=new std::string(*dmol_inputfile);
      *filename = (*filename) + ".outmol";
      try{
          in= new std::ifstream(filename->c_str(),std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );
          std::string str;
          std::vector<std::string> strVec;
          while(!in->eof())
          {
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str, "opt==")>0)
             {
                  boost::algorithm::split(strVec,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
                  energy_Vect->push_back(std::stod(strVec[2]));
             }
             if(boost::algorithm::contains(str, NormalLabel)>0)
                  break;
          }
      }catch(const std::ifstream::failure& e){
          Log::Error<<e.what()<<" readFinalEnergy_CExeDMol"<<std::endl;
          exit(-1);
      }
      in->close();
      double temEnergy = energy_Vect->at(energy_Vect->size()-1);

      delete energy_Vect;
      delete filename;

      return temEnergy;
}
void CExeDMol::getRelaxedGeometryCoord()
{
     m_pIO->m_pPeriodicFramework->clear();
     std::string* filename=new std::string(*dmol_inputfile);
      *filename = (*filename) + ".car";
     m_pIO->input(*filename);

     delete filename;
}


}
