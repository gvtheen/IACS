#include "CGenome.h"
#include "CRealgene.h"
#include "CGraygene.h"
#include "CBinarygene.h"
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
CGenome::CGenome(std::vector <GeneVAR>&vec_genVal,E_CODE_TYPE codeType,size_t index)
{
	this->m_codeType = codeType;
	this->init(vec_genVal);
	this->m_index = index;
}
CGenome::CGenome(CGenome& myGenome)
{
	this->init(myGenome.m_GeneVARofGenome);
	std::vector <double> realValue;
	myGenome.getDecValue(realValue);
    this->updateDecValueGene(realValue);
    this->m_index = myGenome.index();
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
    if( this->m_codeType!=REAL )   //clear all data.
        this->m_totalgeneofGenome.clear();
    else
        this->m_totalRealofGenome.clear();

	for(size_t i=0;i<this->m_geneNum;i++)
    {
        m_Genome[i]->updatecode(myDecValue[i]);
        //collect gene into total gen of genome
        if( this->m_codeType!=REAL )
            this->concatenateBitsets(this->m_totalgeneofGenome,m_Genome[i]->bitGene());
        else
            this->m_totalRealofGenome.push_back(m_Genome[i]->realGene());
    }
}

void CGenome::updateTotalGeneToIndividualGene()
{
    if(this->m_codeType == REAL)
    {
        for(size_t i=0;i<this->m_geneNum;i++)
            m_Genome[i]->realGene()=this->m_totalRealofGenome[i];
    }else{
/*
m_totalgeneofGenome = gene1 gene2 gene3 gene4....geneN;
m_totalgeneofGenome[0]= geneN[0]              // Warnning:    Importantly!!!!!
*/
        size_t index=0;
        for(size_t i=this->m_geneNum-1;i>=0;i--)
        {
             size_t bit_Num = this->m_Genome[i]->bitNum();
            Bitset& bitGene = this->m_Genome[i]->bitGene();
            for(size_t j=0;j<bit_Num;j++)
                bitGene[j]=this->m_totalgeneofGenome[index++];
            this->m_Genome[i]->decode();
        }
    }
}
void CGenome::getDecValue(std::vector <double>& realValue)
{
	std::vector<CGenebase*>::iterator it;
	for(it=this->m_Genome.begin(); it<this->m_Genome.end(); it++)
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
void CGenome::setTotalbitGene(const Bitset my_totalbitgene)
{
    this->m_totalgeneofGenome = my_totalbitgene;
}

//
std::vector <double>& CGenome::totalrealGene()
{
    return this->m_totalRealofGenome;
}
void CGenome::setTotalrealGene(std::vector <double>& my_totalrealgene)
{
    this->m_totalRealofGenome.assign(my_totalrealgene.begin(),my_totalrealgene.end());
}

//
size_t CGenome::geneNum()
{
	return this->m_geneNum;
}
size_t CGenome::totalbitNum()
{
	return this->m_totalbitNum;
}
std::vector <GeneVAR>& CGenome::GeneVARiable()
{
    return this->m_GeneVARofGenome;
}
CGenebase* CGenome::extractGene(int i)
{
	if(i<0)
	{
	    Log::Error<<i <<" is error! CGenome_extractGene"<<std::endl;
	    boost::throw_exception(std::runtime_error("i<0. CGenome_extractGene!\n"));
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
bool CGenome::isNormalFinish()const
{
    return this->FinishState;
}
void CGenome::setFinishState(bool stat)
{
    this->FinishState = stat;
}
//
size_t CGenome::index()
{
   return this->m_index;
}
//
/*
according to the gene-variable array, all of gene objects were initialized
*/
void CGenome::init(std::vector <GeneVAR>& vec_genVal)
{
	this->m_varNumofGenome=vec_genVal.size();
	this->m_geneNum=m_varNumofGenome;
	this->m_totalbitNum=0;

    // set gene-variable array to vector of gene-var
    this->m_GeneVARofGenome.assign(vec_genVal.begin(),vec_genVal.end());

    for(size_t i=0;i<this->m_varNumofGenome;i++)
	{
	    if (this->m_codeType==BINARY)
	         m_Genome.push_back( new CBinarygene(vec_genVal[i]));
        else if (this->m_codeType==GRAY)
             m_Genome.push_back( new CGraygene(vec_genVal[i]));
        else
             m_Genome.push_back( new CRealgene(vec_genVal[i]));
        if( this->m_codeType!=REAL ){
              this->concatenateBitsets(this->m_totalgeneofGenome,m_Genome[i]->bitGene());
        }else
              m_totalRealofGenome.push_back(m_Genome[i]->realGene());
	}
	this->m_totalbitNum = m_totalgeneofGenome.size();
}
void CGenome::concatenateBitsets(Bitset& first, const Bitset& second)
{
    Bitset secondCopy(second);
    secondCopy.resize(first.size() + second.size());
    //Increase the size of the bit buffer to fit the data being placed in it
    first.resize(first.size() + second.size());
    //shift the bits in the firstCopy to the left
    first <<= second.size();
    //"copy" the bits from the secondCopy into the firstCopy
    first |= secondCopy;
/*
  gene1=1010
  gene2=0111
  after doing concatenateBitsets(gene1,gene2)
  gene1=10100111  // gene1 =  original gene1 + gene2;
*/

}
// operator
bool CGenome::operator==(CGenome& myGenome)
{
	if (this->m_geneNum!=myGenome.geneNum() || ( this->m_codeType!=REAL && \
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
	std::vector <double> V_left, V_right;
	double tmp;
	double COMPARE_RADIO=0.05;
	if (this->m_geneNum!=myGenome.geneNum() || ( this->m_codeType!=REAL && \
                                             this->m_totalbitNum!=myGenome.totalbitNum()) )
		return false;

	this->getDecValue(V_left);
	myGenome.getDecValue(V_right);

	for(unsigned int i=0;i<V_left.size();i++)
	{
		 tmp= ( V_left[i] > V_right[i] )? V_left[i] : V_right[i];
		 tmp=fabs(V_left[i] - V_right[i])/tmp;
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
	   Log::Error<<" str of [] can not match! CGenome_[]operator!\n";
	   boost::throw_exception(std::runtime_error("str of [] can not match! CGenome_[]operator!\n"));
	}
}
CGenome::~CGenome()
{
	for(size_t i=0;i<m_Genome.size();i++)
        delete m_Genome[i];
	m_Genome.clear();
	m_GeneVARofGenome.clear();
	m_totalRealofGenome.clear();
	m_totalgeneofGenome.clear();
}

}
