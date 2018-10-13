#include "unistd.h"
#include <fstream>
#include <string.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include "../GaZJUT/GaUtilityFunction.h"
#include "CParameter.h"
#include "../GaZJUT/CGaparameter.h"
#include "../Util/log.hpp"

using util::Log;

namespace CALCZJUT{

CParameter::CParameter(char* file)
:m_pfile(new std::string(file))
{
    //building the connection between cmd and corresponding function
    m_mapCmdFunc["[System_Name]"] = &CParameter::setSysName;
    m_mapCmdFunc["[Evaluator_Code]"] = &CParameter::setEvaluator_Code;
    m_mapCmdFunc["[Evaluator_Benchmark]"] = &CParameter::setEvaluator_Criterion;
    m_mapCmdFunc["[Running_Command]"]=&CParameter::setRunCmd;
    m_mapCmdFunc["[Evaluator_Criterion]"]=&CParameter::setEvaluator_Criterion;
    // fur supported catalyst
    m_mapCmdFunc["[Support_Structure]"]=&CParameter::setSupport_Structure;
    m_mapCmdFunc["[Adsorbent_Structure]"]=&CParameter::setAdsorbent_Structure;
    m_mapCmdFunc["[Adsorbent_Support_Structure]"]=&CParameter::setAdsorbent_Support_Structure;
    // for cluster
    m_mapCmdFunc["[Cluster_Formula]"]=&CParameter::setCluster_Formula;
    m_mapCmdFunc["[Cluster_Structure]"]=&CParameter::setCluster_Input_File;
    // for bond Tolerance
    m_mapCmdFunc["[Bond_Tolerance_Factor]"]=&CParameter::setBond_Tolerance_Factor;
    m_mapCmdFunc["[Exclude_Bond]"]=&CParameter::setExclude_Bond;
    // for GA
    m_mapCmdFunc["[Generation_Number]"]=&CParameter::setGenNum;
    m_mapCmdFunc["[Population_Size]"]=&CParameter::setPopSize;
    m_mapCmdFunc["[Cross_Probability]"]=&CParameter::setPc;
    m_mapCmdFunc["[Mutation_Probability]"]=&CParameter::setPm;
    m_mapCmdFunc["[Scaling_Mode]"]=&CParameter::setScaling_Mode;
    m_mapCmdFunc["[Mutation_Mode]"]=&CParameter::setMutation_Mode;
    m_mapCmdFunc["[Gene_Code]"]=&CParameter::setGene_Code;
    m_mapCmdFunc["[Cross_Mode]"]=&CParameter::setCross_Mode;
    m_mapCmdFunc["[Select_Mode]"]=&CParameter::setSelect_Mode;
    m_mapCmdFunc["[Gene_Formation_Mode]"]=&CParameter::setGene_Formation_Mode;

    m_pGAParameter=new GAZJUT::CGaparameter();
}
void CParameter::input()
{
    if(m_pfile==nullptr)
       m_pfile = new std::string("input.dat");
    if(access(m_pfile->c_str(),F_OK) != 0 )
      {
         Log::Error<<*m_pfile <<" file is no exist! input_CParameter\n";
         boost::throw_exception(std::runtime_error(*m_pfile +"input file is no exist! Check the file: Error_information.txt."));
      }
      std::ifstream *in;
      std::string str;
      std::vector<std::string> cmd_Str;
      std::vector<std::string> keyValue;
      try{
          in= new std::ifstream(m_pfile->c_str(),std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );

          Log::Output<<"*********** Read input file of "<<m_pfile->c_str()<<"***********"<< std::endl;

          while(!in->eof()){
            std::getline(*in,str,'\n');
            boost::algorithm::trim(str);
            if(str=="") continue;
            boost::algorithm::split(cmd_Str,str,boost::algorithm::is_any_of("#"),boost::algorithm::token_compress_on);
            for(size_t i=0;i<cmd_Str.size();i++){
                boost::algorithm::trim(cmd_Str[i]);
                if(cmd_Str[i]!="" && i%2==0){
                     Log::Output<< ">>>" <<cmd_Str[i]<< std::endl;
                     boost::algorithm::split(keyValue,cmd_Str[i],boost::algorithm::is_any_of("="),boost::algorithm::token_compress_on);

                     if(keyValue.size()!=2){
                        Log::Error<<cmd_Str[i] <<" command is wrong!" <<std::endl;
                        boost::throw_exception(std::runtime_error(cmd_Str[i] + " command is wrong! Check the file: Error_information.txt."));
                     }

                     boost::algorithm::trim(keyValue[0]);
                     boost::algorithm::trim(keyValue[1]);
                     if(m_mapCmdFunc.find(keyValue[0]) != m_mapCmdFunc.end()){
                         (this->*m_mapCmdFunc[keyValue[0]])(keyValue[1]);
                     }else{
                         Log::Error<< cmd_Str[i] <<" command is wrong!" <<std::endl;
                         boost::throw_exception(std::runtime_error(cmd_Str[i] + " command is wrong! Check the file: Error_information.txt."));
                     }
                 }
             }
         }
      }catch(const std::ifstream::failure& e){
          Log::Error<< e.what() <<"Input_CParameter\n";
          boost::throw_exception(std::runtime_error(std::string(e.what()) + "Check the file: Error_information.txt."));
      }
      Log::Output<<"*********** End read "<<"***********"<< std::endl;
}
bool CParameter::checkIsValidParameters()
{
     if( m_pGAParameter->EvaluateEXE()==GAZJUT::GAUSSIAN &&
        this->simulationMode == CParameter::MOL_2DMATERIAL )
            return false;

     if( simulationMode != CLUSTER && adso_supp_Struct==""
        && (adso_supp_Struct=="" || adsorbentStructFile==""))
            return false;

     return true;
}
void CParameter::setSysName(std::string mtr)
{
    this->sysName = mtr;
}
void CParameter::setPopSize(std::string mtr)
{
    (*m_pGAParameter)["[Population_Size]"]=mtr;
}
void CParameter::setPm(std::string mtr)
{
    (*m_pGAParameter)["[Mutation_Probability]"]=mtr;
}
void CParameter::setPc(std::string mtr)
{
    (*m_pGAParameter)["[Mutation_Probability]"]=mtr;
}
void CParameter::setGenNum(std::string mtr)
{
    (*m_pGAParameter)["[Cross_Probability]"]=mtr;
}
void CParameter::setRunCmd(std::string mtr)
{
    boost::algorithm::split(runCmd,mtr,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
    for(size_t i=0;i<runCmd.size();i++)
        boost::algorithm::trim(runCmd[0]);
}
void CParameter::setBond_Tolerance_Factor(std::string mtr)
{
    std::vector<std::string> vectStr;
    boost::algorithm::split(vectStr,mtr,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
    if(vectStr.size()<2){
        Log::Error<<mtr << " command is wrong! setBond_Tolerance_Factor_CParameter!\n";
        boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: Error_information.txt."));
    }
    bondToleranceFactor.first =std::stod(vectStr[0]);
    bondToleranceFactor.second=std::stod(vectStr[1]);
    if(bondToleranceFactor.first > bondToleranceFactor.second )
    {
        double temp;
        temp = bondToleranceFactor.first;
        bondToleranceFactor.first = bondToleranceFactor.second;
        bondToleranceFactor.second = temp;
    }
}
void CParameter::setExclude_Bond(std::string mtr)
{
    std::vector<std::string> vectStr;
    boost::algorithm::split(vectStr,mtr,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
    if(vectStr.size()<2){
        Log::Error<<mtr << " command is wrong! setExclude_Bond_CParameter!\n";
        boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: Error_information.txt."));
    }
    excludeBond.push_back(std::pair<std::string*,std::string*>(new std::string(vectStr[0]),\
                                                               new std::string(vectStr[1])));
}
void CParameter::setEvaluator_Code(std::string mtr)
{
/*
# 1: VASP
# 2: GAUSSIAN
# 3: DMOL
# 4: CASTEP
# 5: LAMMPS   */
      if(mtr=="VASP" || std::stoi(mtr)==1)
        (*m_pGAParameter)["[Evaluator_Code]"]="VASP";
      else if(mtr=="GAUSSIAN" || std::stoi(mtr)==2)
        (*m_pGAParameter)["[Evaluator_Code]"]="GAUSSIAN";
      else if(mtr=="DMOL" || std::stoi(mtr)==3)
        (*m_pGAParameter)["[Evaluator_Code]"]="DMOL";
      else if(mtr=="CASTEP" || std::stoi(mtr)==4)
        (*m_pGAParameter)["[Evaluator_Code]"]="CASTEP";
      else if(mtr=="LAMMPS" || std::stoi(mtr)==5)
        (*m_pGAParameter)["[Evaluator_Code]"]="LAMMPS";
      else{
         Log::Error << mtr <<" command is wrong! setEvaluator_Code_CParameter\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: Error_information.txt."));
      }
}
void CParameter::setOutput_struct_format(std::string mtr)
{
/*
## 1: poscar   //.poscar
## 2: mol      //.mol
## 3: cif      //.cif
## 4: car      //.car
## 5: gjf      //.gjf
## 6: cell     //.cell
*/
     if(mtr=="poscar" || std::stoi(mtr)==1)
        this->output_struct_format = "poscar";
      else if(mtr=="mol" || std::stoi(mtr)==2)
        this->output_struct_format = "mol";
      else if(mtr=="cif" || std::stoi(mtr)==3)
        this->output_struct_format = "cif";
      else if(mtr=="car" || std::stoi(mtr)==4)
        this->output_struct_format = "car";
      else if(mtr=="gjf" || std::stoi(mtr)==5)
        this->output_struct_format = "gjf";
      else if(mtr=="cell" || std::stoi(mtr)==6)
        this->output_struct_format = "cell";
      else{
          Log::Error<< mtr << " format isnot supported! setOutput_struct_format_CParameter!\n";
          boost::throw_exception(std::runtime_error(mtr+ "  format isnot supported! Check the file: Error_information.txt."));
      }
}
void CParameter::setEvaluator_Criterion(std::string mtr)
{
     if( mtr=="Energy" || std::stoi(mtr)==1 )
         this->evaluatorCriterion = CParameter::ENERGY;
     else if( mtr=="Force" || std::stoi(mtr)==2 )
         this->evaluatorCriterion = CParameter::FORCE;
     else if( mtr=="Band_gap" || std::stoi(mtr)==3 )
         this->evaluatorCriterion = CParameter::BAND_GAP;
     else{
         Log::Error<< mtr << " command is wrong! setEvaluator_Criterion_CParameter!\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: Error_information.txt."));
     }
}
void CParameter::setMutation_Mode(std::string mtr)
{
     if(mtr=="UNIFORM_M" || std::stoi(mtr)==1 )
        (*m_pGAParameter)["[Mutation_Mode]"]="UNIFORM_M";
     else if(mtr=="BOUNDARY" || std::stoi(mtr)==2)
        (*m_pGAParameter)["[Mutation_Mode]"]="BOUNDARY";
     else if(mtr=="NOUNIFORM" || std::stoi(mtr)==3)
        (*m_pGAParameter)["[Mutation_Mode]"]="NOUNIFORM";
     else if(mtr=="GAUSSIAN_M" || std::stoi(mtr)==4)
        (*m_pGAParameter)["[Mutation_Mode]"]="GAUSSIAN_M";
     else{
         Log::Error<<mtr << " command is wrong! setMutation_Mode_CParameter\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: Error_information.txt."));
     }
}
void CParameter::setGene_Code(std::string mtr)
{
/*
1: BINARY ; 2:GRAY;  3: REAL
*/
     if(mtr=="BINARY " || std::stoi(mtr)==1)
        (*m_pGAParameter)["[Gene_Code]"]="BINARY ";
     else if(mtr=="GRAY" || std::stoi(mtr)==2)
        (*m_pGAParameter)["[Gene_Code]"]="GRAY";
     else if(mtr=="REAL" || std::stoi(mtr)==3)
        (*m_pGAParameter)["[Gene_Code]"]="REAL";
     else{
         Log::Error<<mtr << " command is wrong! setGene_Code_CParameter!\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: Error_information.txt."));
     }
}
void CParameter::setCross_Mode(std::string mtr)
{
/*
#1: SINGLE
#2: MULTIPLE
#3: UNIFORM_C
#4: ARITHMETIC
#5: UNARITHMETIC
*/
     if(mtr=="SINGLE" || std::stoi(mtr)==1 )
        (*m_pGAParameter)["[Cross_Mode]"]="SINGLE";
     else if(mtr=="MULTIPLE" || std::stoi(mtr)==2)
        (*m_pGAParameter)["[Cross_Mode]"]="MULTIPLE";
     else if(mtr=="UNIFORM_C" || std::stoi(mtr)==3)
        (*m_pGAParameter)["[Cross_Mode]"]="UNIFORM_C";
     else if(mtr=="UNARITHMETIC" || std::stoi(mtr)==4)
        (*m_pGAParameter)["[Cross_Mode]"]="UNARITHMETIC";
     else{
         Log::Error<<mtr << " command is wrong! setCross_Mode_CParameter!\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: Error_information.txt."));
     }
}
void CParameter::setSelect_Mode(std::string mtr)
{
/*##
## 1:  ROULETTE_WHEEL
## 2:  TOURNAMENT
## 3:  RANDOM
## 4:  mixed
*/
     if(mtr=="ROULETTE_WHEEL" || std::stoi(mtr)==1 )
        (*m_pGAParameter)["[Search_Mode]"] = "ROULETTE_WHEEL";
     else if(mtr=="TOURNAMENT" || std::stoi(mtr)==2)
        (*m_pGAParameter)["[Search_Mode]"] = "TOURNAMENT";
     else if(mtr=="RANDOM" || std::stoi(mtr)==3)
        (*m_pGAParameter)["[Search_Mode]"] = "RANDOM";
     else if(mtr=="MIXED" || std::stoi(mtr)==4)
        (*m_pGAParameter)["[Search_Mode]"] = "MIXED";
     else{
         Log::Error<< mtr << " command is wrong! setSelect_Mode_CParameter!\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: Error_information.txt."));
     }
}
void CParameter::setGene_Formation_Mode(std::string mtr)
{
/*
#1: RANDOM
#2: FILE
*/
     if(mtr=="RANDOM" || std::stoi(mtr)==1)
        (*m_pGAParameter)["[Gene_Formation_Mode]"]="RANDOM";
     else if(mtr=="FILE" || std::stoi(mtr)==2)
        (*m_pGAParameter)["[Gene_Formation_Mode]"]="FILE";
     else{
         Log::Error<< mtr << " command is wrong! setGene_Formation_Mode_CParameter!\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: Error_information.txt."));
     }

}
void CParameter::setSupport_Structure(std::string str)
{
    this->supportStructFile=str;
}
void CParameter::setAdsorbent_Structure(std::string str)
{
    this->adsorbentStructFile=str;
}
void CParameter::setAdsorbent_Support_Structure(std::string str)
{
    std::vector<std::string> file_vect;
    boost::algorithm::trim(str);

    boost::algorithm::split(file_vect,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
    for(size_t i=0;i<file_vect.size();i++)
       this->adso_supp_Input_File.push_back(new std::string(file_vect[i]));

}
void CParameter::setCluster_Formula(std::string str)
{
    this->cluster_Formula=str;
}
void CParameter::setCluster_Input_File(std::string str)
{
    std::vector<std::string> file_vect;
    boost::algorithm::trim(str);

    boost::algorithm::split(file_vect,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);

    for(size_t i=0;i<file_vect.size();i++)
        this->cluster_Input_File.push_back(new std::string(file_vect[i]));
}
GAZJUT::CGaparameter* CParameter::GaParameter()
{
    return m_pGAParameter;
}
size_t CParameter::currentGenerationNum()
{
    return this->m_pGAParameter->GenerationNum();
}
size_t CParameter::popNum()
{
    return this->m_pGAParameter->PopNum();
}
CParameter::~CParameter()
{
    if(m_pfile!=nullptr)
       delete m_pfile;
}




}
