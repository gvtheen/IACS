#include "CGenome.h"
#include "GaUtilityFunction.h"
#include "../GACatalyst.h"
#include "../Util/log.hpp"

namespace GAZJUT{

CGenome::CGenome()
{
	this->m_geneNum=0;
	this->m_varNumofGenome=0;
	this->m_totalbitNum=0;
	//pointer
	this->m_pGeneVARofGenome=nullptr;
	this->m_pGenome=nullptr;
	this->m_ptotalRealofGenome=nullptr;
}
CGenome::CGenome(std::vector <GeneVAR>*vec_genVal,E_CODE_TYPE codeType)
{
	this->m_codeType = codeType;
	this->init(vec_genVal);
}
CGenome::CGenome(CGenome& myGenome)
{
	this->init(myGenome.m_pGeneVARofGenome);
	std::vector <double> realValue;
	myGenome.getDecValue(realValue);
    this->updateDecValueGene(realValue);
}
//copy object
CGenome* CGenome::clone()
{
	return new CGenome(*this);
}

//init
void  CGenome::init()
{

}
///
void CGenome::updateDecValueGene(std::vector <double>& myDecValue)
{
	assert(m_ptotalRealofGenome);
	assert(m_ptotalgeneofGenome);

    int index=0;
	for(int i=0;i<this->m_geneNum;i++)
    {
        m_pGenome->at(i)->updatecode(myDecValue[i]);
        //collect gene into total gen of genome
        if( this->m_codeType!=REAL ){
             std::vector<unsigned int>* tempbitGene=m_pGenome->at(i)->bitGene();
             if(m_ptotalgeneofGenome->empty()==true)
                m_ptotalgeneofGenome->insert(m_ptotalgeneofGenome->end(),tempbitGene->begin(), \
                                         tempbitGene->end());
             else
                for(unsigned int j=0;j<tempbitGene->size();j++)
                    m_ptotalgeneofGenome->at(index++) = tempbitGene->at(j);
        }else{
             if(m_ptotalRealofGenome->empty()==true)
                  m_ptotalRealofGenome->push_back(m_pGenome->at(i)->realGene());
             else
                  m_ptotalRealofGenome->at(i) = m_pGenome->at(i)->realGene();
        }

    }
}

void CGenome::updateTotalGeneToIndividualGene()
{
    if(this->m_codeType == REAL)
    {
        for(int i=0;i<this->m_geneNum;i++)
            (*m_pGenome)[i]->realGene()=this->m_ptotalRealofGenome->at(i);
    }else{
        int index=0;
        for(int i=0;i<this->m_geneNum;i++)
        {
            int bit_Num = this->m_pGenome->at(i)->bitNum();
            std::vector<unsigned int> *bitGene = this->m_pGenome->at(i)->bitGene();
            for(int j=0;j<bit_Num;j++)
                bitGene->at(j)=this->m_ptotalgeneofGenome->at(index++);
            this->m_pGenome->at(i)->decode();
        }
    }
}
void CGenome::getDecValue(std::vector <double>& realValue)
{
	std::vector<CGenebase*>::iterator it;
	for(it=this->m_pGenome->begin(); it<this->m_pGenome->end(); it++)
        realValue.push_back((*it)->value());
}
//
void CGenome::setFitness(const double fit)
{
	this->m_fitness=fit;
}
double CGenome::fitness()
{
	return this->m_fitness;
}
//
void CGenome::setRawScore(const double rawcore)
{
	this->m_rawscore=rawcore;
}
double CGenome::rawscore()
{
	return this->m_rawscore;
}
//
void CGenome::setOrigValue(const double oriValue)
{
    this->m_origValue = oriValue;
}
double CGenome::origValue()
{
    return this->m_origValue;
}
//
void  CGenome::setCumufitness(const double value)
{
	this->m_cumufitness=value;
}
double CGenome::cumufitness()
{
	return this->m_cumufitness;
}
//
double CGenome::relativefitness()
{
	return this->m_relativefitness;
}
void  CGenome::setRelativefitness(const double value)
{
	this->m_relativefitness=value;
}
//return two type of gene
Bitset& CGenome::totalbitGene()
{

	return this->m_totalgeneofGenome;
}
void CGenome::setTotalbitGene(Bitset my_totalbitgene)
{
    this->m_totalgeneofGenome = my_totalbitgene;
}

//
std::vector <double>* CGenome::totalrealGene()
{
    return this->m_ptotalRealofGenome;
}
void CGenome::setTotalrealGene(std::vector <double>& my_totalrealgene)
{
    assert(this->m_ptotalRealofGenome);
    this->m_ptotalRealofGenome->assign(my_totalrealgene.begin(),my_totalrealgene.end());
}

//
int CGenome::geneNum()
{
	return this->m_geneNum;
}
int CGenome::totalbitNum()
{
	return this->m_totalbitNum;
}
std::vector <GeneVAR>* CGenome::GeneVARiable()
{
    return this->m_pGeneVARofGenome;
}
CGenebase* CGenome::extractGene(int i)
{
	if(i<1)
	{
	    util::Log::Error<<"i<1. CGenome_extractGene"<<std::endl;
	    std::exit(-1);
	}
	return this->m_pGenome->at(i-1);
}
void  CGenome::insertGeneToGenome(CGenebase* mygene)
{
	this->m_geneNum=this->m_geneNum+1;
	this->m_totalbitNum=this->m_geneNum + mygene->bitNum();
	this->m_pGenome->push_back(mygene);
	this->m_pGeneVARofGenome->push_back(*(mygene->m_GeneVAR));
	for(int i=0;i<mygene->bitNum();i++)
	   this->m_ptotalgeneofGenome->push_back(mygene->bitNum());
}
bool CGenome:isNormalFinish()const
{
    return this->FinishState;
}
void CGenome:setFinishState(bool stat)
{
    this->FinishState = stat;
}
//

//
void CGenome::init(std::vector <GeneVAR>* vec_genVal)
{
	assert(vec_genVal);

	this->m_varNumofGenome=vec_genVal->size();
	this->m_geneNum=m_varNumofGenome;
	this->m_totalbitNum=0;
	// pointer
	//General_GA_Parameter.CodeMode=REAL;
	this->m_pGeneVARofGenome=new (std::vector<GeneVAR>)();   //contruct empty vector<GeneVAR> pointer
    this->m_pGeneVARofGenome->assign(vec_genVal->begin(),vec_genVal->end());

    this->m_pGenome = new (std::vector<CGenebase*>)(m_geneNum);

    //for 3 type of gene,initialization.
    if( General_GA_Parameter.CodeMode != REAL ){
        m_ptotalgeneofGenome = new (std::vector <unsigned int>);
    }else
        m_ptotalRealofGenome = new (std::vector <double>);

    for(int i=0;i<this->m_varNumofGenome;i++)
	{
	    if (General_GA_Parameter.CodeMode==BINARY)
	         m_pGenome->at(i) = new CBinarygene(vec_genVal->at(i));
        else if (General_GA_Parameter.CodeMode==GRAY)
             m_pGenome->at(i) = new CGraygene(vec_genVal->at(i));
        else
             m_pGenome->at(i) = new CRealgene(vec_genVal->at(i));
        if( General_GA_Parameter.CodeMode!=REAL ){
            std::vector<unsigned int>* tempbitGene=m_pGenome->at(i)->bitGene();
            m_ptotalgeneofGenome->insert(m_ptotalgeneofGenome->end(),tempbitGene->begin(), \
                                         tempbitGene->end());
            this->m_totalbitNum=this->m_totalbitNum + m_pGenome->at(i)->bitNum();
        }else
            m_ptotalRealofGenome->push_back(m_pGenome->at(i)->realGene());
	}
}
// operator
bool CGenome::operator==(CGenome& myGenome)
{
	bool res=true;
	if (this->m_geneNum!=myGenome.geneNum() || (General_GA_Parameter.CodeMode!=REAL && \
                                             this->m_totalbitNum!=myGenome.totalbitNum()) )
	   return false;
	for(unsigned int i=0; i<m_ptotalgeneofGenome->size(); i++)
	   if(m_ptotalgeneofGenome->at(i)!=(myGenome.totalbitGene())->at(i))
	   {
	   	    res=false;
	   	    break;
	   }
	return res;
}
bool CGenome::operator ^= (CGenome& myGenome)
{
	bool res=true;
	std::vector <double> *V_left, *V_right;
	double tmp;
	double COMPARE_RADIO=0.05;
	if (this->m_geneNum!=myGenome.geneNum() || (General_GA_Parameter.CodeMode!=REAL && \
                                             this->m_totalbitNum!=myGenome.totalbitNum()) )
		return false;

	V_left=this->getDecValue();
	V_right=myGenome.getDecValue();

	for(unsigned int i=0;i<V_left->size();i++)
	{
		 tmp= ((V_left->at(i)) > (V_right->at(i)))? V_left->at(i) : V_right->at(i);
		 tmp=fabs(V_left->at(i) - V_right->at(i))/tmp;
		 if(tmp > COMPARE_RADIO){
		 	res=false;
		 	break;
		 }
	}
	return res;
}
double CGenome::operator[](std::string index_name)
{
	if(index_name=="rawscore")
	   return this->m_rawscore;
	else if(index_name=="fitness")
	   return this->m_fitness;
	else if(index_name=="origvalue")
	   return this->m_origValue;
	else{
	   ERROR_OUTPUT("Error: str of [] can not match","CGenome","[]operator");
	   return 0;
	}
}
CGenome::~CGenome()
{
	for(size_t i=0;i<m_pGenome->size();i++)
        delete (*m_pGenome)[i];
	delete m_pGenome;
	delete m_pGeneVARofGenome;
}

}
