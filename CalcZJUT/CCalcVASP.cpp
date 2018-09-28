#include "unistd.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "../GaZJUT/GaUtilityFunction.h"
#include "CCalcVASP.h"
#include "CCalc2DSupport.h"
#include "CIOPoscar.h"
#include "CCalcCluster.h"
#include "CCalcClusterSupport.h"
#include "../CataZJUT/CFragment.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../Util/Point-Vector.h"
#include "../GaZJUT/CGaparameter.h"
#include "../Util/log.hpp"
#include "../CataZJUT/CBondTolerance.h"

using util::Log;

namespace CATAZJUT{
  class CFragment;
  class CPeriodicFramework;
}

using util::Point3;

namespace CALCZJUT{


CCalcVASP::CCalcVASP(CParameter* mpara):CCalcFitnessInterface(mpara)
{
    m_pInputFile.push_back(new std::string("INCAR"));
    m_pInputFile.push_back(new std::string("KPOINTS"));
    m_pInputFile.push_back(new std::string("POSCAR"));
    m_pInputFile.push_back(new std::string("POTCAR"));
}

CCalcVASP::~CCalcVASP()
{
    for(size_t i=0;i<m_pInputFile.size();i++)
        delete m_pInputFile[i];
    m_pInputFile.clear();
}
void CCalcVASP::init()
{
     if(m_Parameter->simulationMode !=CParameter::CLUSTER)
     {
        if(m_Parameter->simulationMode ==CParameter::MOL_2DMATERIAL){
             this->m_pCalcModeStruct=new CCalc2DSupport(this->m_Parameter);
        }else{
             this->m_pCalcModeStruct=new CCalcClusterSupport(this->m_Parameter);
        }

        m_pIO = new CIOPoscar( m_pCalcModeStruct->periodicFramework());

        if(m_Parameter->adso_supp_Struct != ""){
           //read coordinate of mixed adso_supp_Struct, add the pointer of m_pCalcModeStruct;
           //
           Cios* tempIO = this->getIO(m_Parameter->adso_supp_Struct,m_pCalcModeStruct->periodicFramework());
           tempIO->input(m_Parameter->adso_supp_Struct);
           delete tempIO;
           //construct chemical bond
           m_pCalcModeStruct->periodicFramework()->perceiveBonds();
           //construct chemical fragments to identify support and adsorbent
           m_pCalcModeStruct->periodicFramework()->perceiveFragments();
           //
           //if(m2dSupport->periodicFramework()->fragments().size())
           CATAZJUT::CFragment *m1,*m2;
           m1 = m_pCalcModeStruct->periodicFramework()->fragment(0);
           m2 = m_pCalcModeStruct->periodicFramework()->fragment(1);
           if(m1== nullptr || m2 ==nullptr){
               Log::Error<<"Support and adsorbent setting of input file are error! input_CCalcVASP.\n";
               boost::throw_exception(std::runtime_error("Support and adsorbent setting are error! Check the file: Error_information.txt."));
           }
           if(m1->atomCount() > m2->atomCount())
           {
                m_pCalcModeStruct->createSupport(m1->bitSet());
                m_pCalcModeStruct->createMoleAdsorb(m2->bitSet());
           }else{
                m_pCalcModeStruct->createSupport(m2->bitSet());
                m_pCalcModeStruct->createMoleAdsorb(m1->bitSet());
           }
         }else{ // treat cluster-support
            //read coordinate of mixed adso_supp_Struct, add the pointer of m_pCalcModeStruct;
           //
           Cios* tempIO = this->getIO(m_Parameter->supportStructFile,m_pCalcModeStruct->periodicFramework());
           tempIO->input(m_Parameter->supportStructFile);
           delete tempIO;

           Cios* tempIO_1 = this->getIO(m_Parameter->adsorbentStructFile,m_pCalcModeStruct->periodicFramework());
           Bitset Bit_adsorbent = tempIO_1->input(m_Parameter->adsorbentStructFile,CParameter::MOL_CLUSTER);
           delete tempIO_1;
           // get opposite bit for support
           m_pCalcModeStruct->createSupport(~Bit_adsorbent);
           //set molecular adsorbent bits
           m_pCalcModeStruct->createMoleAdsorb(Bit_adsorbent);
         }
     }else{    //pure cluster model
        CCalcCluster* cluster = new CCalcCluster(this->m_Parameter);
        m_pIO = new CIOPoscar( cluster->periodicFramework());
        this->m_pCalcModeStruct = cluster;
        m_pIO->input("POSCAR");
     }

    //initialize gene varible
    //
    m_Parameter->GaParameter()->setGeneVar( m_pCalcModeStruct->GeneVarRange());

    m_pCalcModeStruct->periodicFramework()->m_pBondEvaluator->setExcludeBond(m_Parameter->excludeBond);

    m_pCalcModeStruct->periodicFramework()->m_pBondEvaluator->setTolerancefactor(m_Parameter->bondToleranceFactor);

    if(m_Parameter->output_struct_format=="")
       m_Parameter->output_struct_format="poscar";   //default value;
}
void CCalcVASP::ConvOrigToRawScore(std::vector<double>* temporgValue)
{
    std::vector<double> tmpValue;
    tmpValue.assign(temporgValue->begin(),temporgValue->end());
    std::vector<double> ::iterator maxEnergy = std::max_element(tmpValue.begin(),tmpValue.end(),[](double a,double b){return a < b;});
    for(size_t i=0;i<tmpValue.size();i++)
         temporgValue->at(i) = *maxEnergy - tmpValue[i];
}
double CCalcVASP::CalcuRawFit(std::vector<double>* RealValueOfGenome,size_t& pop_index, bool& isNormalexist)
{
     pid_t pid;
     double res;
     //
     //construct new object of CPeriodicFramework class
     if( pop_index < m_Parameter->GaParameter()->PopNum() )
     {
        m_pCalcModeStruct->createStructureAtGene();
     }
     //transfer gene value to POSCAR file
     m_pCalcModeStruct->setGeneValueToStruct(*RealValueOfGenome);
     m_pIO->setConfiguration(m_pCalcModeStruct->periodicFramework(pop_index));
     m_pIO->output("POSCAR");
     //
     //Check whether input files is OK?
     CheckInputFile();

    // pid=fork();
     if (pid<0){
         Log::Error<<"Building new process is error! CalcuRawFit_CCalcVASP"<<std::endl;
         boost::throw_exception(std::runtime_error("Building new process is error! Check the file: Error_information.txt."));//ERROR TREATMENT;
     }else if(pid==0){
         //run VASP program
     }else
         wait(NULL);
     //
     //read CONTCAR file;
     getRelaxedGeometryCoord();
     std::string out_filename("CONTCAR_");

     if(IsNormalComplete()==true){
        isNormalExist=true;
        res=readFinalEnergy();
        out_filename = out_filename + std::to_string(m_Parameter->m_pGAParameter->Curr_Generation) + \
                   "_" + std::to_string(pop_index);
     }else{
        isNormalExist=false;
        res=9999999;
        out_filename = out_filename + std::to_string(m_Parameter->m_pGAParameter->Curr_Generation) + \
                   "_" + std::to_string(pop_index) + "_ERROR";
     }

     // Use the output file format to output structure;
     Cios* tempIO = this->getIO(m_Parameter->output_struct_format,m_pCalcModeStruct->periodicFramework());
     out_filename = out_filename + m_Parameter->output_struct_format;
     tempIO->output(out_filename);
     delete tempIO;

     // monitor whether all of tasks is completed, especically for parallel running
     pop_run_state[pop_index]=1;
     if( pop_run_state.flip().none()==true ) // all of gene is complete.
     {
         m_pCalcModeStruct->removeStructureOfGene();
         pop_run_state.reset();  //set to 0;
     }
     return res;
}
void CCalcVASP::CheckInputFile()
{
    bool res=true;
    for(unsigned int i=0;i<this->m_pInputFile->size();i++)
        if(access((m_pInputFile->at(i)).c_str(),F_OK) != 0)
        {
           std::string error_info = m_pInputFile->at(i) + std::string(" file is not exist!");
           Log::Error<<error_info <<" CheckInputFile_CCalcVASP\n";
           res=false;
        }
    if(res==false)
    {
         Log::Error<<"Input files in VASP is error!CheckInputFile_CCalcVASP"<<std::endl;
         boost::throw_exception(std::runtime_error("Input files in VASP is error! Check the file: Error_information.txt."));
    }
}
bool CCalcVASP::IsNormalComplete()
{
      std::string NormalLabel("General timing and accounting informations");

      bool res=false;
      if(access("OUTCAR",F_OK) != 0 )
      {
           Log::Error<<"OUTCAR file is no exist! IsNormalComplete_CCalcVASP\n";
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
          Log::Error<<e.what()<<" IsNormalComplete_CCalcVASP\n";
          exit(-1);
      }
      in->close();
      return res;
}
double CCalcVASP::readFinalEnergy()
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
          Log::Error<<e.what()<<"readFinalEnergy_CCalcVASP\n";
          exit(-1);
      }
      in->close();
      double temEnergy = energy_Vect->at(energy_Vect->size()-1);

      delete energy_Vect;

      return temEnergy;
}
void CCalcVASP::getRelaxedGeometryCoord()
{
     m_pIO->m_pPeriodicFramework->clear();
     m_pIO->input("CONTCAR");
}



}
