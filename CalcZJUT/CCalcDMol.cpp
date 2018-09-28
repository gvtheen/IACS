#include "unistd.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "../GaZJUT/GaUtilityFunction.h"
#include "CCalcDMol.h"
#include "CParameter.h"
#include "CIOCar.h"
#include "CCalc2DSupport.h"
#include "CCalcCluster.h"
#include "CCalcClusterSupport.h"
#include "../CataZJUT/CFragment.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../Util/Point-Vector.h"
#include "../Util/log.hpp"
#include "CGaparameter.h"

using util::Log;

namespace CALCZJUT{

CCalcDMol::CCalcDMol(CParameter* mpara)
:CCalcFitnessInterface(mpara)
{
    dmol_inputfile =new std::string(m_Parameter->runCmd[m_Parameter->runCmd.size()-1]);

    m_pInputFile.push_back(new std::string(*dmol_inputfile+".car"));
    m_pInputFile.push_back(new std::string(*dmol_inputfile+".input"));
    m_pInputFile.push_back(new std::string(*dmol_inputfile+".mdf"));
}
CCalcDMol::~CCalcDMol()
{
    if(dmol_inputfile!=nullptr)
        delete dmol_inputfile;

    for(size_t i=0;i<m_pInputFile.size();i++)
        delete m_pInputFile[i];
    m_pInputFile.clear();
}
void CCalcDMol::init()
{
     #ifdef DEBUG
       Log::Debug<<" Initialize CCalcDMol object!" <<std::endl;
     #endif // DEBUG

     if(m_Parameter->simulationMode !=CParameter::CLUSTER)
     {
        if(m_Parameter->simulationMode ==CParameter::MOL_2DMATERIAL){
             this->m_pCalcModeStruct=new CCalc2DSupport(this->m_Parameter);
        }else{
             this->m_pCalcModeStruct=new CCalcClusterSupport(this->m_Parameter);
        }

        m_pIO = new CIOCar( m_pCalcModeStruct->periodicFramework());

        if(m_Parameter->adso_supp_Struct != ""){
           //read coordinate of mixed adso_supp_Struct
           Cios* tempIO = this->getIO(m_Parameter->supportStructFile,m_pCalcModeStruct->periodicFramework());
           tempIO->input(m_Parameter->supportStructFile);
           delete tempIO;

           Cios* tempIO_1 = this->getIO(m_Parameter->adsorbentStructFile,m_pCalcModeStruct->periodicFramework());
           Bitset Bit_adsorbent = tempIO_1->input(m_Parameter->adsorbentStructFile,CParameter::MOL_CLUSTER);
           delete tempIO_1;
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
               Log::Error<<"Support and adsorbent setting of input file are error! input_CCalcVASP"<<std::endl;
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
           m_pIO->input(m_Parameter->supportStructFile);

           Bitset Bit_adsorbent = m_pIO->input(m_Parameter->adsorbentStructFile,CParameter::MOL_CLUSTER)
           // get opposite bit for support
           m_pCalcModeStruct->createSupport(~Bit_adsorbent);
           //set molecular adsorbent bits
           m_pCalcModeStruct->createMoleAdsorb(Bit_adsorbent);
         }
     }else{    //pure cluster model
        CCalcCluster* cluster = new CCalcCluster(this->m_Parameter);
        m_pIO = new CIOCar( cluster->periodicFramework());
        this->m_pCalcModeStruct = cluster;
        m_pIO->input("POSCAR");
     }

    //initialize gene varible
    //
    m_Parameter->GaParameter()->setGeneVar( m_pCalcModeStruct->GeneVarRange());

    m_pCalcModeStruct->periodicFramework()->m_pBondEvaluator->setExcludeBond(m_Parameter->excludeBond);

    m_pCalcModeStruct->periodicFramework()->m_pBondEvaluator->updateTolerancefactor(m_Parameter->bondToleranceFactor);

    if(m_Parameter->output_struct_format=="")
       m_Parameter->output_struct_format="mol";   //default value;
}
double CCalcDMol::CalcuRawFit(std::vector<double>* RealValueOfGenome,size_t& pop_index, bool& isNormalexist)
{
     pid_t pid;
     double res;
     size_t currGeneration =m_Parameter->m_pGAParameter->Curr_Generation;

     Log::Info<<" Run DMol calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";

     //construct new object of CPeriodicFramework class
     if( pop_index < m_Parameter->m_pGAParameter->PopNum() )
     {
        m_pCalcModeStruct->createStructureAtGene();
     }
     //transfer gene value to structure file
     m_pCalcModeStruct->setGeneValueToStruct(RealValueOfGenome);
     //
     m_pIO->setConfiguration(m_pCalcModeStruct->periodicFramework(pop_index));
     //write structure in the car file
     m_pIO->output(m_pInputFile[0]+".car"); //.car file
     //
     //Check whether input files is OK?
     CheckInputFile();

     pid=fork();
     if (pid<0){
         Log::Error<<"Building new process is error! CalcuRawFit_CCalcDMol"<<std::endl;
         boost::throw_exception(std::runtime_error("Building new process is error! Check the file: Error_information.txt."));//ERROR TREATMENT;
     }else if(pid==0){
         //run DMOL program
     }else
         wait(NULL);
     //
     //read CONTCAR file;
     getRelaxedGeometryCoord();
     std::string out_filename("relaxed_");

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
void CCalcDMol::ConvOrigToRawScore(std::vector<double>* OrigRawScore)
{
    std::vector<double> tmpValue;
    tmpValue.insert(temporgValue->begin(),temporgValue->end());
    double maxEle = std::max_element(tmpValue.begin(),tmpValue.end(),[](double a,double b){return a<b;})
    for(size_t i=0;i<tmpValue.size();i++);
         temporgValue->at(i)=maxEle-tmpValue[i];
}
void CCalcDMol::CheckInputFile()
{
    bool res=true;
    for(unsigned int i=0;i<this->m_pInputFile->size();i++)
        if(access((m_pInputFile->at(i)).c_str(),F_OK) != 0)
        {
           std::string error_info = m_pInputFile->at(i) + std::string(" file is not exist!");
           Log::Error<<error_info<<" CheckInputFile_CCalcDMol\n";
           res=false;
        }
    if(res==false)
    {
         Log::Error<<"Input files in DMOL is error! CheckInputFile_CCalcDMol"<<std::endl;
         boost::throw_exception(std::runtime_error("Input files in DMOL is error! Check the file."));
    }
}
bool CCalcDMol::IsNormalComplete()
{
      std::string NormalLabel("DMol3 job finished successfully");
      std::string* filename=new std::string(*dmol_inputfile);
      *filename = (*filename) + ".outmol";
      bool res=false;
      if(access(filename->c_str(),F_OK) != 0 )
      {
           Log::Error<<(*filename) <<" file is no exist! IsNormalComplete_CCalcDMol"<<std::endl;
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
          Log::Error<<e.what()<<" IsNormalComplete_CCalcDMol!"<<std::endl;
          exit(-1);
      }
      in->close();
      delete filename;
      return res;
}
double CCalcDMol::readFinalEnergy()
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
          Log::Error<<e.what()<<" readFinalEnergy_CCalcDMol"<<std::endl;
          exit(-1);
      }
      in->close();
      double temEnergy = energy_Vect->at(energy_Vect->size()-1);

      delete energy_Vect;
      delete filename;

      return temEnergy;
}
void CCalcDMol::getRelaxedGeometryCoord()
{
     m_pIO->m_pPeriodicFramework->clear();
     std::string* filename=new std::string(*dmol_inputfile);
      *filename = (*filename) + ".car";
     m_pIO->input(*filename);

     delete filename;
}


}
