#include "CGenome.h"
#include "GaUtilityFunction.h"
#include "../GACatalyst.h"
#include "../Util/log.hpp"

using util::Log;

namespace GAZJUT{

CGenome::CGenome()
{
	this->m_geneNum=0;
	this->m_varNumofGenome=0;
	this->m_totalbitNum=0;
	//pointer
}
CGenome::CGenome(std::vector <GeneVAR>*vec_genVal,E_CODE_TYPE codeType)
{
	this->m_codeType = codeType;
	this->init(*vec_genVal);
}
CGenome::CGenome(CGenome& myGenome)
{
	this->init(myGenome.m_GeneVARofGenome);
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

    int index=0;
	for(int i=0;i<this->m_geneNum;i++)
    {
        m_Genome[i]->updatecode(myDecValue[i]);
        //collect gene into total gen of genome
        if( this->m_codeType!=REAL ){
             std::vector<unsigned int>* tempbitGene=m_Genome[i]->bitGene();
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
	if(i<0)
	{
	    Log::Error<<"i<1. CGenome_extractGene"<<std::endl;
	    boost::throw_exception(std::runtime_error("i<1. CGenome_extractGene!\n"));
	}
	return this->m_Genome[i];
}
void  CGenome::insertGeneToGenome(CGenebase* mygene)
{
	this->m_geneNum=this->m_geneNum+1;
	this->m_totalbitNum=this->m_geneNum + mygene->bitNum();
	this->m_Genome.push_back(mygene);
	this->m_GeneVARofGenome.push_back(*(mygene->m_GeneVAR));

	if(mygene->codeType!= REAL)
	  this->concatenateBitsets(this->m_totalgeneofGenome,mygene->bitGene());
    else
      this->m_totalRealofGenome.push_back(mygene->m_value);

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
void CGenome::init(std::vector <GeneVAR>& vec_genVal)
{
	assert(vec_genVal);

	this->m_varNumofGenome=vec_genVal.size();
	this->m_geneNum=m_varNumofGenome;
	this->m_totalbitNum=0;

    this->m_GeneVARofGenome.assign(vec_genVal.begin(),vec_genVal.end());

    for(int i=0;i<this->m_varNumofGenome;i++)
	{
	    if (General_GA_Parameter.CodeMode==BINARY)
	         m_Genome.push_back( new CBinarygene(vec_genVal[i] ));
        else if (General_GA_Parameter.CodeMode==GRAY)
             m_Genome.push_back( new CGraygene(vec_genVal[i])) );
        else
             m_Genome.push_back( new CRealgene(vec_genVal[i]) );
        if( General_GA_Parameter.CodeMode!=REAL ){
              this->concatenateBitsets(this->m_totalgeneofGenome,m_Genome[i]->bitGene());
        }else
               m_totalRealofGenome.push_back(m_Genome[i]->realGene());
	}
	this->m_totalbitNum = m_totalgeneofGenome.size();
}
void CGenome::concatenateBitsets(boost::dynamic_bitset<>& first, const boost::dynamic_bitset<>& second)
{
    Bitset secondCopy(second);

    secondCopy.resize(first.size() + second.size());
    //Increase the size of the bit buffer to fit the data being placed in it
    first.resize(first.size() + second.size());


    //shift the bits in the firstCopy to the left
    first <<= second.size();

    //"copy" the bits from the secondCopy into the firstCopy
    first |= secondCopy;

}
// operator
bool CGenome::operator==(CGenome& myGenome)
{
	if (this->m_geneNum!=myGenome.geneNum() || (General_GA_Parameter.CodeMode!=REAL && \
                                             this->m_totalbitNum!=myGenome.totalbitNum()) )
	   return false;
    if(m_totalgeneofGenome == myGenome.totalbitGene())
       return true;
    else
	   return false;
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
	   Log::Error<<(" str of [] can not match! CGenome_[]operator!\n";
	   boost::throw_exception(std::runtime_error("str of [] can not match! CGenome_[]operator!\n"));
	}
}
CGenome::~CGenome()
{
	for(size_t i=0;i<m_Genome.size();i++)
        delete m_pGenome[i];
	m_Genome.clear();
	m_GeneVARofGenome.clear();
}

}
