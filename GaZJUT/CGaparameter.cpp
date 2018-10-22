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
#include "CGaparameter.h"
#include "GaUtilityFunction.h"
#include "../Util/log.hpp"
#include "../GACatalyst.h"
//# - Scale methods defaults
/*CDefScaleLinearMultiplier     = 1.2
CDefScaleSigmaTruncMultiplier = 2.0
CDefScalePowerLawFactor       = 1.0005
CDefScaleBoltzMinTemp         = 1.0
CDefScaleBoltzFactor          = 0.05
# 40 temp. = 500 generations
CDefScaleBoltzStart           = 40.0
*/

using util::Log;
//using namespace GAZJUT;

namespace GAZJUT{

CGaparameter::CGaparameter()
{
    this->m_GeneVARofPopulation.clear();
    this->defaultInit();
}

CGaparameter::CGaparameter(std::vector <GAZJUT::GeneVAR>& myVar)
{
    this->m_GeneVARofPopulation.assign(myVar.begin(),myVar.end());
    this->defaultInit();
}
CGaparameter::~CGaparameter()
{
    delete m_mapCmdString;
    m_GeneVARofPopulation.clear();
}

CGaparameter::CGaparameter( CGaparameter& other)
{
    assert(&other);
    assert(other.m_mapCmdString);

    this->m_GeneVARofPopulation.assign(other.GeneVAR().begin(), other.GeneVAR().end());

    this->defaultInit();

    this->Curr_Generation = other.Curr_Generation;

    this->m_mapCmdString = new (std::map <std::string, std::string>);
    this->m_mapCmdString->insert(other.m_mapCmdString->begin(),other.m_mapCmdString->end());

}
CGaparameter& CGaparameter::operator=( CGaparameter& other)
{
    if (this == &other)
        return *this; // handle self assignment
    //assignment operator
    this->Curr_Generation = other.Curr_Generation;

    this->m_mapCmdString = new (std::map <std::string, std::string>);
    this->m_mapCmdString->insert(other.m_mapCmdString->begin(),other.m_mapCmdString->end());

    this->m_GeneVARofPopulation.assign(other.GeneVAR().begin(), other.GeneVAR().end());

    return *this;
}

void CGaparameter::defaultInit()
{
    m_mapCmdString = new (std::map <std::string, std::string>);

    m_mapCmdString->insert(std::pair<std::string,std::string>("[Generation_Number]",    "30"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Population_Size]",      "20"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Cross_Probability]",    "0.80"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Mutation_Probability]", "0.20"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Initial_Gene_File]",     "gene.txt"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Select_Mode]",   "ROULETTE_WHEEL"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Search_Mode]",   "MIN"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Cross_Mode]",    "SINGLE"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Cross_Number]",      "2"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Mutation_Mode]",   "UNIFORM_M"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Evaluator_Code]",     "VASP"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Scaling_Mode]",  "LINEAR"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Gene_Code]",     "BINARY"));
    m_mapCmdString->insert(std::pair<std::string,std::string>("[Gene_Formation_Mode]","FILE_INPUT"));

    this->Curr_Generation = 0;
}
void CGaparameter::add_Curr_Generation()
{
    this->Curr_Generation = this->Curr_Generation + 1;
}
std::vector <GeneVAR>& CGaparameter::GeneVAR()
{
    return this->m_GeneVARofPopulation;
}
size_t CGaparameter::PopNum()
{
     std::string keyValue=(*m_mapCmdString)["[Population_Size]"];
     int d_value;
     try{
        d_value=stoi(keyValue);
     }catch(const std::exception& e){
        Log::Error<< "key value is error:" << e.what() <<" CGaparameter_PopNum!\n";
        boost::throw_exception(std::runtime_error("key value is error! CGaparameter_PopNum!\n"));
     }
     return d_value;
}
size_t CGaparameter::GenerationNum()
{
     std::string keyValue=(*m_mapCmdString)["[Generation_Number]"];
     int d_value;
     try{
        d_value=stoi(keyValue);
     }catch(const std::exception& e){
        Log::Error<< "key value is error:" << e.what() <<" CGaparameter_GenerationNum!\n";
        boost::throw_exception(std::runtime_error("key value is error! CGaparameter_GenerationNum!\n"));
     }
     return d_value;
}
size_t CGaparameter::CrossNum()
{
     std::string keyValue=(*m_mapCmdString)["[Cross_Number]"];
     int d_value;
     try{
        d_value=stoi(keyValue);
     }catch(const std::exception& e){
        Log::Error<< "key value is error:" << e.what() <<" CGaparameter_CrossNum!\n";
        boost::throw_exception(std::runtime_error("key value is error! CGaparameter_CrossNum!\n"));
     }
     return d_value;
}
void CGaparameter::setGeneVAR(std::vector <GAZJUT::GeneVAR>& myVar)
{
    assert(&myVar);

    if( m_GeneVARofPopulation.size()>0 )
       m_GeneVARofPopulation.clear();

    for(size_t i=0;i<myVar.size();i++)
        m_GeneVARofPopulation.push_back(myVar[i]);

    checkGeneVAR();
}
void CGaparameter::checkGeneVAR()
{
    std::vector<GAZJUT::GeneVAR>::iterator it;
    for(it=m_GeneVARofPopulation.begin();it<m_GeneVARofPopulation.end();it++)
    {
        if((it->min) > (it->max))
           std::swap(it->min,it->max);     //call std::swap(x,y) function

        if((it->accuracy)==0.0)
        {
	       Log::Error<< "checkGeneVAR! CGaparameter_checkGeneVAR!\n";
           boost::throw_exception(std::runtime_error("Variable value is error! CGaparameter_checkGeneVAR!\n"));
        }
    }
}
double CGaparameter::CrossProb()
{
     std::string keyValue=(*m_mapCmdString)["[Cross_Probability]"];
     double d_value;
     try{
        d_value=stof(keyValue);
     }catch(const std::exception& e){
        Log::Error<< "key value is error:" << e.what() <<" CGaparameter_CrossProb!\n";
        boost::throw_exception(std::runtime_error("key value is error! CGaparameter_CrossProb!\n"));
     }
     return d_value;
}
double CGaparameter::MutaProb()
{
     std::string keyValue=(*m_mapCmdString)["[Mutation_Probability]"];
     double d_value;
     try{
        d_value=stof(keyValue);
     }catch(const std::exception& e){
        Log::Error<< "key value is error:" << e.what() <<" CGaparameter_MutaProb!\n";
        boost::throw_exception(std::runtime_error("key value is error! CGaparameter_MutaProb!\n"));
     }
     return d_value;
}
std::string CGaparameter::GeneFile()
{
    return (*m_mapCmdString)["[Initial_Gene_File]"];
}
E_GA_TYPE CGaparameter::SearchType()
{
    std::map <std::string, E_GA_TYPE> *mySearch=new (std::map <std::string, E_GA_TYPE>);
    E_GA_TYPE res;
    mySearch->insert(std::pair<std::string,E_GA_TYPE>("MIN",    MIN));
    mySearch->insert(std::pair<std::string,E_GA_TYPE>("MAX",    MAX));
    res = (*mySearch)[(*m_mapCmdString)["[Search_Mode]"]];
    delete mySearch;
    return res;
}
E_SELECT_OPERATOR  CGaparameter::SelectMode()
{
      std::map <std::string, E_SELECT_OPERATOR> *mySearch=new (std::map <std::string, E_SELECT_OPERATOR>);
      E_SELECT_OPERATOR res;
      mySearch->insert(std::pair<std::string,E_SELECT_OPERATOR>("RANDOM",        RANDOM));
      mySearch->insert(std::pair<std::string,E_SELECT_OPERATOR>("TOURNAMENT",    TOURNAMENT));
      mySearch->insert(std::pair<std::string,E_SELECT_OPERATOR>("ROULETTE_WHEEL",ROULETTE_WHEEL));
      mySearch->insert(std::pair<std::string,E_SELECT_OPERATOR>("MIXED",         MIXED));
      res = (*mySearch)[(*m_mapCmdString)["[Select_Mode]"]];
      delete mySearch;
      return res;
}
E_CROSSOVER_OPERATOR CGaparameter::CrossMode()
{
    std::map <std::string, E_CROSSOVER_OPERATOR> *mySearch=new (std::map <std::string, E_CROSSOVER_OPERATOR>);
    E_CROSSOVER_OPERATOR res;
    mySearch->insert(std::pair<std::string,E_CROSSOVER_OPERATOR>("SINGLE",        SINGLE));
    mySearch->insert(std::pair<std::string,E_CROSSOVER_OPERATOR>("MULTIPLE",      MULTIPLE));
    mySearch->insert(std::pair<std::string,E_CROSSOVER_OPERATOR>("UNIFORM_C",     UNIFORM_C));
    mySearch->insert(std::pair<std::string,E_CROSSOVER_OPERATOR>("ARITHMETIC",    ARITHMETIC));
    mySearch->insert(std::pair<std::string,E_CROSSOVER_OPERATOR>("UNARITHMETIC",  UNARITHMETIC));
    res = (*mySearch)[(*m_mapCmdString)["[Cross_Mode]"]];
    delete mySearch;
    return res;
}
E_MUTATE_OPERATOR CGaparameter::MutateMode()
{
    std::map <std::string, E_MUTATE_OPERATOR> *mySearch=new (std::map <std::string,E_MUTATE_OPERATOR>);
    E_MUTATE_OPERATOR res;

    mySearch->insert(std::pair<std::string,E_MUTATE_OPERATOR>("UNIFORM_M",    UNIFORM_M));
    mySearch->insert(std::pair<std::string,E_MUTATE_OPERATOR>("BOUNDARY",     BOUNDARY));
    mySearch->insert(std::pair<std::string,E_MUTATE_OPERATOR>("NOUNIFORM",    NOUNIFORM));
    mySearch->insert(std::pair<std::string,E_MUTATE_OPERATOR>("GAUSSIAN_M",   GAUSSIAN_M));
    res = (*mySearch)[(*m_mapCmdString)["[Mutation_Mode]"]];
    delete mySearch;
    return res;
}
E_GENEFORMATION_TYPE CGaparameter::InitGenMode()
{
    std::map <std::string, E_GENEFORMATION_TYPE> *mySearch=new (std::map <std::string,E_GENEFORMATION_TYPE>);
    E_GENEFORMATION_TYPE res;
    mySearch->insert(std::pair<std::string,E_GENEFORMATION_TYPE>("RANDOM_FORMATION",    RANDOM_FORMATION));
    mySearch->insert(std::pair<std::string,E_GENEFORMATION_TYPE>("FILE_INPUT",FILE_INPUT));
    res=(*mySearch)[(*m_mapCmdString)["[Gene_Formation_Mode]"]];
    delete mySearch;
    return res;
}
E_EVALUATOR_EXE CGaparameter::EvaluateEXE()
{
    std::map <std::string, E_EVALUATOR_EXE> *mySearch=new (std::map <std::string,E_EVALUATOR_EXE>);
    E_EVALUATOR_EXE res;
    mySearch->insert(std::pair<std::string,E_EVALUATOR_EXE>("VASP",    VASP));
    mySearch->insert(std::pair<std::string,E_EVALUATOR_EXE>("GAUSSIAN",GAUSSIAN));
    mySearch->insert(std::pair<std::string,E_EVALUATOR_EXE>("DMOL",    DMOL));
    mySearch->insert(std::pair<std::string,E_EVALUATOR_EXE>("CASTEP",  CASTEP));
    mySearch->insert(std::pair<std::string,E_EVALUATOR_EXE>("LAMMPS",  LAMMPS));
    res =(*mySearch)[(*m_mapCmdString)["[Evaluator_Code]"]];
    delete mySearch;
    return res;
}
E_SCALING_TYPE CGaparameter::ScalingMode()
{
    std::map <std::string, E_SCALING_TYPE> *mySearch=new (std::map <std::string,E_SCALING_TYPE>);
    E_SCALING_TYPE res;
    mySearch->insert(std::pair<std::string,E_SCALING_TYPE>("LINEAR",  LINEAR));
    mySearch->insert(std::pair<std::string,E_SCALING_TYPE>("SIGMA",   SIGMA));
    mySearch->insert(std::pair<std::string,E_SCALING_TYPE>("POWER",   POWER));
    res = (*mySearch)[(*m_mapCmdString)["[Scaling_Mode]"]];
    delete mySearch;
    return res;
}
E_CODE_TYPE CGaparameter::CodeMode()
{
    std::map <std::string, E_CODE_TYPE> *mySearch=new (std::map <std::string,E_CODE_TYPE>);
    E_CODE_TYPE res;

    mySearch->insert(std::pair<std::string,E_CODE_TYPE>("BINARY",  BINARY));
    mySearch->insert(std::pair<std::string,E_CODE_TYPE>("GRAY",    GRAY));
    mySearch->insert(std::pair<std::string,E_CODE_TYPE>("REAL",    REAL));
    res = (*mySearch)[(*m_mapCmdString)["[Gene_Code]"]];
    delete mySearch;
    return res;
}
std::string& CGaparameter::operator[](std::string key_name)
{
    if(m_mapCmdString->find(key_name) == m_mapCmdString->end())
    {
        Log::Error<< "keyname is error: CGaparameter_operator[]!\n";
        boost::throw_exception(std::runtime_error("keyname is error! CGaparameter_operator[]!\n"));
    }
    return (*m_mapCmdString)[key_name];
}
void CGaparameter::setKeyValue(const std::string key,const std::string value)
{

}




}
