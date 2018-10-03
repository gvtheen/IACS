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
#include "CCalcGaussian.h"
#include "../CataZJUT/CFragment.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../Util/Point-Vector.h"
#include "../GaZJUT/CGaparameter.h"
#include "../Util/log.hpp"
#include "../CataZJUT/CBondTolerance.h"

using util::Log;

namespace CALCZJUT{

//
CCalcGaussian::CCalcGaussian(CParameter* mpara)
:CCalcFitnessInterface(mpara)
{
    m_pInputFile =new std::string(m_Parameter->runCmd[m_Parameter->runCmd.size()-1]);
}
//
CCalcGaussian::~CCalcGaussian()
{
    if(m_pInputFile!=nullptr)
      delete m_pInputFile;
}
void CCalcGaussian::init()
{
     if(m_Parameter->simulationMode ==CParameter::MOL_CLUSTER)
     {
         this->m_pCalcModeStruct=new CCalcClusterSupport(this->m_Parameter);

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
               Log::Error<<"Support and adsorbent setting of input file are error! input_CCalcVASP\n";
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
         }else{
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
     }else if (m_Parameter->simulationMode ==CParameter::CLUSTER){    //pure cluster model
        CCalcCluster* cluster = new CCalcCluster(this->m_Parameter);
        m_pIO = new CIOPoscar( cluster->periodicFramework());
        this->m_pCalcModeStruct = cluster;
        m_pIO->input("POSCAR");
     }

    //initialize gene varible
    //
    m_Parameter->GaParameter()->setGeneVAR( m_pCalcModeStruct->GeneVARRange());

    m_pCalcModeStruct->periodicFramework()->m_pBondEvaluator->setExcludeBond(m_Parameter->excludeBond);

    m_pCalcModeStruct->periodicFramework()->m_pBondEvaluator->setTolerancefactor(m_Parameter->bondToleranceFactor);

    if(m_Parameter->output_struct_format=="")
       m_Parameter->output_struct_format="gjf";   //default value;
}
double CCalcGaussian::CalcuRawFit(std::vector<double>* RealValueOfGenome,size_t& pop_index, bool& isNormalExist)
{
         pid_t pid;
     double res;
     //
     size_t currGeneration =m_Parameter->GaParameter()->Curr_Generation;

     Log::Info<<" Run Gaussian calculation of the "<< pop_index<< "th Genome in "<< currGeneration <<"th generation!\n";
     //construct new object of CPeriodicFramework class
     if( pop_index < m_Parameter->GaParameter()->PopNum() )
     {
        m_pCalcModeStruct->createStructureAtGene();
     }
     //transfer gene value to POSCAR file
     m_pCalcModeStruct->setGeneValueToStruct(*RealValueOfGenome);
     m_pIO->setConfiguration(m_pCalcModeStruct->periodicFramework(pop_index));
     m_pIO->output(*m_pInputFile+".gjf"); //.mol file
     //
     //Check whether input files is OK?
     CheckInputFile();

     //pid=fork();
     if (pid<0){
         Log::Error<<"Building new process is error! CalcuRawFit_CCalcDMol\n";
         boost::throw_exception(std::runtime_error("Building new process is error! Check the file: Error_information.txt."));//ERROR TREATMENT;
     }else if(pid==0){
         //run VASP program
     }else
        ;// wait(NULL);
     //
     //read CONTCAR file;
     getRelaxedGeometryCoord();
     std::string out_filename("relaxed_");

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
void CCalcGaussian::ConvOrigToRawScore(std::vector<double>* temporgValue)
{

}
void CCalcGaussian::CheckInputFile()
{
   if(access(m_pInputFile->c_str(),F_OK) != 0)
    {
         std::string error_info = *(m_pInputFile) + std::string(" file is not exist!");
         Log::Error<<error_info<<" CheckInputFile_CCalcGaussian\n";
         boost::throw_exception(std::runtime_error("Input files in Gaussian isnot exist! Check the file."));
    }
}
bool CCalcGaussian::IsNormalComplete()
{

}
double CCalcGaussian::readFinalEnergy()
{

}
void CCalcGaussian::getRelaxedGeometryCoord()
{

}

}
