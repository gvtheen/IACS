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
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "CParameter.h"
#include "../GaZJUT/CGaparameter.h"
#include "../Util/log.hpp"
#include "../Util/utilFunction.h"
#include "../IACS.h"

using util::Log;
using util::strcasecmp;

namespace CALCZJUT{

CParameter::CParameter(std::string file)
{
    //building the connection between cmd and corresponding function
    m_mapCmdFunc["[System_Name]"] = &CParameter::setSysName;
    m_mapCmdFunc["[Evaluator_Code]"] = &CParameter::setEvaluator_Code;
    m_mapCmdFunc["[Simulation_Mode]"] = &CParameter::setSimulationMode;
    m_mapCmdFunc["[Running_Command]"]=&CParameter::setRunCmd;
    m_mapCmdFunc["[Evaluator_Criterion]"]=&CParameter::setEvaluator_Criterion;
    m_mapCmdFunc["[Search_Mode]"]=&CParameter::setSearch_Mode;
    m_mapCmdFunc["[Output_Structure_Format]"]=&CParameter::setOutputStructureFormat;
    // fur supported catalyst
    m_mapCmdFunc["[Support_Structure]"]=&CParameter::setSupport_Structure;
    m_mapCmdFunc["[Adsorbent_Structure]"]=&CParameter::setAdsorbent_Structure;
    m_mapCmdFunc["[Adsorbent_Support_Structure]"]=&CParameter::setAdsorbent_Support_Structure;
    m_mapCmdFunc["[Distance_Adsorbent_On_Support]"]=&CParameter::setDistanceRangeAdsorbentOnSupport;
    // for cluster
    m_mapCmdFunc["[Cluster_Formula]"]=&CParameter::setCluster_Formula;
    m_mapCmdFunc["[Cluster_Structure]"]=&CParameter::setCluster_Input_File;
    // for bond Tolerance
    m_mapCmdFunc["[Bond_Tolerance_Factor]"]=&CParameter::setBond_Tolerance_Factor;
    m_mapCmdFunc["[Exclude_Bond]"]=&CParameter::setExclude_Bond;
    m_mapCmdFunc["[Bond_Tolerance]"]=&CParameter::setBond_Tolerance;
    // for GA
    m_mapCmdFunc["[Generation_Number]"]=&CParameter::setGenNum;
    m_mapCmdFunc["[Population_Size]"]=&CParameter::setPopSize;
    m_mapCmdFunc["[Cross_Probability]"]=&CParameter::setPc;
    m_mapCmdFunc["[Mutation_Probability]"]=&CParameter::setPm;
    m_mapCmdFunc["[Scaling_Mode]"]=&CParameter::setScaling_Mode;
    m_mapCmdFunc["[Mutation_Mode]"]=&CParameter::setMutation_Mode;
    m_mapCmdFunc["[Gene_Code]"]=&CParameter::setGene_Code;
    m_mapCmdFunc["[Cross_Mode]"]=&CParameter::setCross_Mode;
    m_mapCmdFunc["[Cross_Number]"]=&CParameter::setCross_Num;
    m_mapCmdFunc["[Select_Mode]"]=&CParameter::setSelect_Mode;
    m_mapCmdFunc["[Gene_Formation_Mode]"]=&CParameter::setGene_Formation_Mode;

    m_pGAParameter=new GAZJUT::CGaparameter();
    this->m_pfile=new std::string(file);
}
CParameter::~CParameter()
{
    if(m_pfile!=nullptr)
       delete m_pfile;
    m_mapCmdFunc.clear();
}

//**********Public function*********/
void CParameter::input()
{
    Log::Info<<"Reading input file..."<<std::endl;

    if(m_pfile==nullptr)
       m_pfile = new std::string("input.dat");

    if(access(m_pfile->c_str(),F_OK) != 0 )
      {
         Log::Error<<*m_pfile <<" file is no exist! input_CParameter\n";
         boost::throw_exception(std::runtime_error(*m_pfile +"input file is no exist! CParameter::input()."));
      }
      std::ifstream *in=nullptr;
      std::string str;
      std::vector<std::string> cmd_Str;
      std::vector<std::string> keyValue;
      try{
          in= new std::ifstream(m_pfile->c_str(),std::ifstream::in);
          in->exceptions (std::ifstream::badbit );
          #ifdef DEBUG
            Log::Debug<<"*********** Read input file of "<<*m_pfile<<"***********"<< std::endl;
          #endif
          if(in==nullptr || !in->is_open()){
            Log::Error<<*m_pfile <<" file open with an error! input_CParameter\n";
            boost::throw_exception(std::runtime_error(*m_pfile +" file open with an error! Check the file: CParameter::input()."));
          }
          while( !in->eof()){
            std::getline(*in,str,'\n');
            boost::algorithm::trim(str);
            if(str=="") continue;

            cmd_Str.clear();

            boost::algorithm::split(cmd_Str,str,boost::algorithm::is_any_of("#"),boost::algorithm::token_compress_on);
            for(size_t i=0;i<cmd_Str.size();i++){
                boost::algorithm::trim(cmd_Str[i]);
                if(cmd_Str[i]!="" && i%2==0){
                     #ifdef DEBUG
                        Log::Debug<<cmd_Str[i]<< std::endl;
                     #endif // DEBUG
                     boost::algorithm::split(keyValue,cmd_Str[i],boost::algorithm::is_any_of("="),boost::algorithm::token_compress_on);

                     if(keyValue.size()!=2){
                        Log::Error<<cmd_Str[i] <<" command is wrong! see code of CParameter::input()!" <<std::endl;
                        boost::throw_exception(std::runtime_error(cmd_Str[i] + " command is wrong! see code of CParameter::input()!"));
                     }

                     boost::algorithm::trim(keyValue[0]);
                     boost::algorithm::trim(keyValue[1]);
                     if(m_mapCmdFunc.find(keyValue[0]) != m_mapCmdFunc.end()){
                         // call pointer function
                         (this->*m_mapCmdFunc[keyValue[0]])(keyValue[1]);
                     }else{
                         Log::Error<< cmd_Str[i] <<" command is wrong!" <<std::endl;
                         boost::throw_exception(std::runtime_error(cmd_Str[i] + " command is wrong! see code of CParameter::input()!"));
                     }
                 }
             }
         }
      }catch(const std::ifstream::failure& e){
          Log::Error<< e.what() <<"CParameter::input()\n";
          boost::throw_exception(std::runtime_error(std::string(e.what()) + "Check the file: CParameter::input()."));
      }
      in->close();
      Log::Info<<"END Reading input file"<<std::endl;

}
void CParameter::output()
{
    Log::Output<<"***********Parameters of GA***********"<<std::endl;

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

// set work environment
void CParameter::initWorkEnvironment()
{
    /** \brief  create working environment
    * \path1: scratch
    * \path2: work          create  pop_num  folder
    * \path3: parameter     initial files
    */
    boost::filesystem::path working_path,scratch_path,parameter_path;
      working_path = m_root_WorkingPath/"work";
    if(!boost::filesystem::exists(working_path))
        boost::filesystem::create_directory(working_path);

      scratch_path = m_root_WorkingPath/"scratch";
    if(!boost::filesystem::exists(scratch_path))
        boost::filesystem::create_directory(scratch_path);

    parameter_path = m_root_WorkingPath/"parameter";
    if(!boost::filesystem::exists(parameter_path))
        boost::filesystem::create_directory(parameter_path);

}
void CParameter::setCurrentWorkPathAt(size_t generation,size_t population)
{
    std::stringstream str;
    str<<"Gen_"<<generation<<"Pop_"<<population;
    std::string pathstr;
    str>>pathstr;
    boost::filesystem::path working_path;
    working_path = m_root_WorkingPath/"work";

    working_path=working_path/pathstr.c_str();

    if(! boost::filesystem::exists(working_path))
         boost::filesystem::create_directory(working_path);

    boost::filesystem::current_path(working_path);
}
void CParameter::setCurrentWorkPathAt(CParameter::PATH_TYPE mType)
{
    boost::filesystem::path working_path;
    switch(int(mType)){
       case CParameter::WORK:
           working_path=m_root_WorkingPath/"work";
           break;
       case CParameter::PARAMETER:
           working_path=m_root_WorkingPath/"parameter";
           break;
       case CParameter::SCRATCH:
           working_path=m_root_WorkingPath/"scratch";
           break;
       case CParameter::MAINPATH:
       default:
           working_path=m_root_WorkingPath;
           break;
    }
    boost::filesystem::current_path(working_path);
}
std::string CParameter::rootPath()
{
    return this->m_root_WorkingPath.string();
}
std::string CParameter::scratchPath()
{
    boost::filesystem::path scratch_path = m_root_WorkingPath/"scratch";
    return scratch_path.string();
}
std::string CParameter::parameterPath()
{
    boost::filesystem::path  parameter_path;
    parameter_path = m_root_WorkingPath/"parameter";
    return parameter_path.string();
}
std::string CParameter::currentWorkPath()
{
    return boost::filesystem::current_path().string();
}
void CParameter::checkExeNecessaryFiles(std::vector<std::string*>& files)
{
    boost::filesystem::path tempPath,tempPath_2;
    for(size_t i=0;i<files.size();i++){
       tempPath=m_root_WorkingPath/files[i]->c_str();

       tempPath_2 = m_root_WorkingPath/"parameter";
       tempPath_2 = tempPath_2/files[i]->c_str();
       if(!boost::filesystem::is_regular_file(tempPath)||!boost::filesystem::is_regular_file(tempPath_2)
         ||boost::filesystem::file_size(tempPath) ==0 || boost::filesystem::file_size(tempPath_2) ==0){
           Log::Error<<*files[i]<<" parameter file is not exist or empty in current director and parameter folder!"<<std::endl;
           boost::throw_exception(std::runtime_error(*files[i] +" parameter file is no exist or empty! CParameter::checkExeNecessaryFiles!\n"));
       }else{
           boost::filesystem::path parameter_path = m_root_WorkingPath/"parameter";
           this->copyFileToPath(tempPath.string(),parameter_path.string());
       }
    }
}
void CParameter::moveFileToPath(const std::string& file, const std::string& dir)
{
    boost::filesystem::path oldFilePath(file);

    if(!boost::filesystem::is_regular_file(oldFilePath)){
        Log::Error<<file<<" parameter file is not exist or empty!"<<std::endl;
        boost::throw_exception(std::runtime_error(file +" parameter file is no exist or empty! CParameter::moveFileToPath!\n"));
    }

    boost::filesystem::path newDirPath(dir);

    if(! boost::filesystem::exists(newDirPath))
         boost::filesystem::create_directory(newDirPath);

    boost::filesystem::path newFilePath(dir + "/" + oldFilePath.leaf().string());
    if(!boost::filesystem::is_regular_file(newFilePath) ||
       boost::filesystem::file_size(newFilePath) ==0 )
       boost::filesystem::rename(oldFilePath,newFilePath);
}

void CParameter::copyFileToPath(const std::string& file, const std::string& dir)
{
    boost::filesystem::path oldFilePath(file);

    if(!boost::filesystem::is_regular_file(oldFilePath)){
        Log::Error<<file<<" parameter file is not exist or empty!"<<std::endl;
        boost::throw_exception(std::runtime_error(file +" parameter file is no exist or empty! CParameter::moveFileToPath!\n"));
    }
    boost::filesystem::path newDirPath(dir);

    if(! boost::filesystem::exists(newDirPath))
         boost::filesystem::create_directory(newDirPath);

    boost::filesystem::path newFilePath(dir + "/" + oldFilePath.leaf().string());

    boost::filesystem::copy_file(oldFilePath,newFilePath);
}

//**********End Public function*********/

//**********Private function*********/

void CParameter::setSysName(std::string mtr)
{
    this->sysName = mtr;
}

void CParameter::setRunCmd(std::string mtr)
{
    boost::algorithm::split(runCmd,mtr,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
    for(size_t i=0;i<runCmd.size();i++)
        boost::algorithm::trim(runCmd[0]);
}
void CParameter::setSimulationMode(std::string mtr)
{
     int res=std::stoi(mtr);
     if (res<1 || res>5 ){
        Log::Error<<" Simulation mode setting is error!"<<std::endl;
        boost::throw_exception(std::runtime_error(" Simulation mode setting is error!! CParameter::setSimulationMode!\n"));
     }
     switch (res){
       case 1:
           this->simulationMode=CParameter::CLUSTER;
           break;
       case 2:
           this->simulationMode=CParameter::MOL_CLUSTER;
           break;
       case 3:
           this->simulationMode=CParameter::MOL_2DMATERIAL;
           break;
       case 4:
           this->simulationMode=CParameter::MOL_CLUSTER2DMATERIAL;
           break;
       case 5:
           this->simulationMode=CParameter::PERIODIC;
           break;
       default:
           this->simulationMode=CParameter::CLUSTER;
           break;
     }
}
void CParameter::setBond_Tolerance_Factor(std::string mtr)
{
    std::vector<std::string> vectStr;
    boost::algorithm::split(vectStr,mtr,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
    if(vectStr.size()<2){
        Log::Error<<mtr << " command is wrong! setBond_Tolerance_Factor_CParameter!\n";
        boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: CParameter::setBond_Tolerance_Factor."));
    }
    bondToleranceFactor.first =std::stod(vectStr[0]);
    bondToleranceFactor.second=std::stod(vectStr[1]);
    if(bondToleranceFactor.first > bondToleranceFactor.second)
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
        boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: CParameter::setExclude_Bond_CParameter."));
    }
    this->excludeBond.push_back(std::pair<std::string*,std::string*>(new std::string(vectStr[0]),\
                                                               new std::string(vectStr[1])));
}
void CParameter::setBond_Tolerance(std::string mtr)
{
    std::vector<std::string> vectStr;
    boost::algorithm::split(vectStr,mtr,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
    if(vectStr.size()!=4){
       Log::Error<<mtr << " command is wrong! setBond_Toleranc!\n";
       boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: CParameter::setBond_Tolerance."));
    }
    double d1=std::stof(vectStr[2]);
    double d2=std::stof(vectStr[3]);
    CATAZJUT::CBondPrivate* tmpBondTol=new CATAZJUT::CBondPrivate(vectStr[0],vectStr[1],
                                                                  d1,d2);
    this->bondTolerance.push_back(tmpBondTol);
}
void CParameter::setEvaluator_Code(std::string mtr)
{
/*
# 1: VASP
# 2: GAUSSIAN
# 3: DMOL
# 4: CASTEP
# 5: LAMMPS   */
      if(strcasecmp(mtr,"VASP") || std::stoi(mtr)==1)
        m_pGAParameter->setKeyValue("[Evaluator_Code]","VASP");
      else if(strcasecmp(mtr,"GAUSSIAN") || std::stoi(mtr)==2)
        m_pGAParameter->setKeyValue("[Evaluator_Code]","GAUSSIAN");
      else if(strcasecmp(mtr,"DMOL") || std::stoi(mtr)==3)
        m_pGAParameter->setKeyValue("[Evaluator_Code]","DMOL");
      else if(strcasecmp(mtr,"CASTEP") || std::stoi(mtr)==4)
        m_pGAParameter->setKeyValue("[Evaluator_Code]","CASTEP");
      else if(strcasecmp(mtr,"LAMMPS") || std::stoi(mtr)==5)
        m_pGAParameter->setKeyValue("[Evaluator_Code]","LAMMPS");
      else if(strcasecmp(mtr,"DFTB") || std::stoi(mtr)==6)
        m_pGAParameter->setKeyValue("[Evaluator_Code]","DFTB");
      else{
         Log::Error << mtr <<" command is wrong! setEvaluator_Code_CParameter\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: setEvaluator_Code_CParameter!"));
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
     if(strcasecmp(mtr,"poscar") || std::stoi(mtr)==1)
        this->output_struct_format = "poscar";
      else if( strcasecmp(mtr,"mol") || std::stoi(mtr)==2)
        this->output_struct_format = "mol";
      else if(strcasecmp(mtr,"cif") || std::stoi(mtr)==3)
        this->output_struct_format = "cif";
      else if(strcasecmp(mtr,"car") || std::stoi(mtr)==4)
        this->output_struct_format = "car";
      else if(strcasecmp(mtr,"gjf") || std::stoi(mtr)==5)
        this->output_struct_format = "gjf";
      else if(strcasecmp(mtr,"cell") || std::stoi(mtr)==6)
        this->output_struct_format = "cell";
      else if(strcasecmp(mtr,"xyz") || std::stoi(mtr)==7)
        this->output_struct_format = "xyz";
      else{
          Log::Error<< mtr << " format isnot supported! setOutput_struct_format_CParameter!\n";
          boost::throw_exception(std::runtime_error(mtr+ "  format is not supported! Check the file: CParameter::setOutput_struct_format!"));
      }
}
void CParameter::setEvaluator_Criterion(std::string mtr)
{
     if( strcasecmp(mtr,"Energy") || std::stoi(mtr)==1 )
         this->evaluatorCriterion = CParameter::ENERGY;
     else if( strcasecmp(mtr,"Force") || std::stoi(mtr)==2 )
         this->evaluatorCriterion = CParameter::FORCE;
     else if( strcasecmp(mtr,"Band_gap") || std::stoi(mtr)==3 )
         this->evaluatorCriterion = CParameter::BAND_GAP;
     else{
         Log::Error<< mtr << " command is wrong! setEvaluator_Criterion_CParameter!\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: CParameter::setEvaluator_Criterion."));
     }
}
void CParameter::setSearch_Mode(std::string mtr)
{
    int res=std::stoi(mtr);
    if(res !=1 && res !=2 ){
       Log::Error<< mtr << " command is wrong! CParameter::setSearch_Mode!\n";
       boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: CParameter::setSearch_Mode."));
    }
    if(res==1)
        this->searchMode=CParameter::MIN;
    else
        this->searchMode=CParameter::MAX;
}
// Commands for GA Part
void CParameter::setPopSize(std::string mtr)
{
    m_pGAParameter->setKeyValue("[Population_Size]",mtr);
}
void CParameter::setPm(std::string mtr)
{
    m_pGAParameter->setKeyValue("[Mutation_Probability]",mtr);
}
void CParameter::setPc(std::string mtr)
{
    m_pGAParameter->setKeyValue("[Mutation_Probability]",mtr);
}
void CParameter::setGenNum(std::string mtr)
{
    m_pGAParameter->setKeyValue("[Cross_Probability]",mtr);
}
void CParameter::setScaling_Mode(std::string mtr)
{
    m_pGAParameter->setKeyValue("[Scaling_Mode]",mtr);
}
void CParameter::setMutation_Mode(std::string mtr)
{
     if(strcasecmp(mtr,"UNIFORM_M") || std::stoi(mtr)==1 )
        m_pGAParameter->setKeyValue("[Mutation_Mode]","UNIFORM_M");
     else if(strcasecmp(mtr,"BOUNDARY") || std::stoi(mtr)==2)
        m_pGAParameter->setKeyValue("[Mutation_Mode]","BOUNDARY");
     else if( strcasecmp(mtr,"NOUNIFORM") || std::stoi(mtr)==3)
        m_pGAParameter->setKeyValue("[Mutation_Mode]","NOUNIFORM");
     else if(strcasecmp(mtr,"GAUSSIAN_M") || std::stoi(mtr)==4)
        m_pGAParameter->setKeyValue("[Mutation_Mode]","GAUSSIAN_M");
     else{
         Log::Error<<mtr << " command is wrong! setMutation_Mode_CParameter\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: setMutation_Mode_CParameter."));
     }
}
void CParameter::setGene_Code(std::string mtr)
{
/*
1: BINARY ; 2:GRAY;  3: REAL
*/
     if(strcasecmp(mtr,"BINARY") || std::stoi(mtr)==1)
        m_pGAParameter->setKeyValue("[Gene_Code]","BINARY");
     else if(strcasecmp(mtr,"GRAY") || std::stoi(mtr)==2)
        m_pGAParameter->setKeyValue("[Gene_Code]","GRAY");
     else if(strcasecmp(mtr,"REAL") || std::stoi(mtr)==3)
        m_pGAParameter->setKeyValue("[Gene_Code]","REAL");
     else{
         Log::Error<<mtr << " command is wrong! setGene_Code_CParameter!\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: CParameter::setGene_Code."));
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
     if(strcasecmp(mtr,"SINGLE") || std::stoi(mtr)==1 )
        m_pGAParameter->setKeyValue("[Cross_Mode]","SINGLE");
     else if(strcasecmp(mtr,"MULTIPLE") || std::stoi(mtr)==2)
        m_pGAParameter->setKeyValue("[Cross_Mode]","MULTIPLE");
     else if(strcasecmp(mtr,"UNIFORM_C") || std::stoi(mtr)==3)
        m_pGAParameter->setKeyValue("[Cross_Mode]","UNIFORM_C");
     else if(strcasecmp(mtr,"UNARITHMETIC") || std::stoi(mtr)==4)
        m_pGAParameter->setKeyValue("[Cross_Mode]","UNARITHMETIC");
     else{
         Log::Error<<mtr << " command is wrong! setCross_Mode_CParameter!\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: CParameter::setCross_Mode!"));
     }
}
void CParameter::setCross_Num(std::string mtr)
{
   m_pGAParameter->setKeyValue("[Cross_Number]",mtr);
}
void CParameter::setSelect_Mode(std::string mtr)
{
/*##
## 1:  ROULETTE_WHEEL
## 2:  TOURNAMENT
## 3:  RANDOM
## 4:  mixed
*/
     if(strcasecmp(mtr,"ROULETTE_WHEEL") || std::stoi(mtr)==1 )
        m_pGAParameter->setKeyValue("[Search_Mode]","ROULETTE_WHEEL");
     else if(strcasecmp(mtr,"TOURNAMENT") || std::stoi(mtr)==2)
        m_pGAParameter->setKeyValue("[Search_Mode]", "TOURNAMENT");
     else if(strcasecmp(mtr,"RANDOM") || std::stoi(mtr)==3)
        m_pGAParameter->setKeyValue("[Search_Mode]", "RANDOM");
     else if(strcasecmp(mtr,"MIXED") || std::stoi(mtr)==4)
        m_pGAParameter->setKeyValue("[Search_Mode]", "MIXED");
     else{
         Log::Error<< mtr << " command is wrong! CParameter::setSelect_Mode!\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: CParameter::setSelect_Mode."));
     }
}
void CParameter::setGene_Formation_Mode(std::string mtr)
{
/*
#1: RANDOM
#2: FILE
*/
     if(strcasecmp(mtr,"RANDOM") || std::stoi(mtr)==1)
        m_pGAParameter->setKeyValue("[Gene_Formation_Mode]","RANDOM");
     else if(strcasecmp(mtr,"FILE") || std::stoi(mtr)==2)
        m_pGAParameter->setKeyValue("[Gene_Formation_Mode]","FILE");
     else{
         Log::Error<< mtr << " command is wrong! CParameter::setGene_Formation_Mode!\n";
         boost::throw_exception(std::runtime_error(mtr+ " value format is wrong! Check the file: CParameter::setGene_Formation_Mode."));
     }

}

/*Setting command for Initial structure*/
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
void CParameter::setDistanceRangeAdsorbentOnSupport(std::string str)
{
    std::vector<std::string> str_vect;
    boost::algorithm::trim(str);

    boost::algorithm::split(str_vect,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
    if(str_vect.size()!=2){
      Log::Error<< str << " command is wrong! CParameter::setDistanceRangeAdsorbentOnSupport!\n";
      boost::throw_exception(std::runtime_error(str+ " value format is wrong! Check the file: CParameter::setDistanceRangeAdsorbentOnSupport."));
    }
    this->DistRange_Adsorbent_Support.first=std::stof(str_vect[0]);
    this->DistRange_Adsorbent_Support.second=std::stof(str_vect[1]);

    if(DistRange_Adsorbent_Support.first>DistRange_Adsorbent_Support.second){
       double tmp= DistRange_Adsorbent_Support.first;
       DistRange_Adsorbent_Support.first=DistRange_Adsorbent_Support.second;
       DistRange_Adsorbent_Support.second=tmp;
    }
}
void CParameter::setOutputStructureFormat(std::string mtr)
{
   int res=std::stoi(mtr);
   if(res<1 || res>7){
     Log::Error<< mtr << " command is wrong! CParameter::setOutputStructureFormat!\n";
     boost::throw_exception(std::runtime_error(mtr+ " command format is wrong! Check the file: CParameter::setOutputStructureFormat."));
   }
   switch(res){
     case 1:
         this->output_struct_format="poscar";
         break;
     case 2:
         this->output_struct_format="mol";
         break;
     case 3:
         this->output_struct_format="cif";
         break;
     case 4:
         this->output_struct_format="car";
         break;
     case 5:
         this->output_struct_format="gjf";
         break;
     case 6:
         this->output_struct_format="cell";
         break;
     case 7:
         this->output_struct_format="xyz";
         break;
     default:
         this->output_struct_format="poscar";
         break;
   }
}
bool CParameter::checkIsValidParameters()
{
     if( m_pGAParameter->EvaluateEXE()==GAZJUT::GAUSSIAN &&
        this->simulationMode == CParameter::MOL_2DMATERIAL )
            return false;

     if(( simulationMode != CLUSTER && simulationMode != PERIODIC )&&
        ( adso_supp_Input_File.size()==0 || ( supportStructFile=="" && adsorbentStructFile=="" )))
            return false;

     return true;
}


}
