#include<cmath>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include "CGpopulation.h"
namespace GAZJUT{

CGpopulation::CGpopulation()
{
    m_minRawGenome  = nullptr;
    m_maxRawGenome  = nullptr;
    m_ObjGaparameter= nullptr;
}
CGpopulation::CGpopulation(CGaparameter* my_Gaparameter)
{
	assert(my_Gaparameter);
    assert(my_Gaparameter->PopNum());

	//firstly built empty vector<GENEVAR>
	this->m_pObjGaparameter=my_Gaparameter;

	for(size_t i=0;i<my_Gaparameter->PopNum();i++)
    {
        CGenome *it = new CGenome(m_ObjGaparameter->m_pGeneVarofPopulation,m_ObjGaparameter->CodeMode());
        m_Gpopulation.push_back(it);
    }

    //until this time, both best genome and worst genome are unknown.
    m_minRawGenome = nullptr;
    m_maxRawGenome = nullptr;
}
CGpopulation::CGpopulation(CGpopulation& mypopulation)
{
	//assert(mypopulation);
    assert(mypopulation.m_pObjGaparameter);

    // pointer array
    this->m_pObjGaparameter= new CGaparameter(*(mypopulation.m_pObjGaparameter));
	CGenome* tempCGenome=nullptr;
	for(i=0;i<m_ObjGaparameter->PopNum();i++)
    {
        tempCGenome = new CGenome(*(mypopulation[i]));
        m_Gpopulation.push_back(tempCGenome);
    }
//    this->m_pMinRawGenome=new CGenome(*(mypopulation.m_minRawGenome));
//    this->m_pMaxRawGenome=new CGenome(*(mypopulation.m_maxRawGenome));
}
CGpopulation::~CGpopulation()
{
    for(size_t i=0;i<m_Gpopulation.size();i++)
	    delete m_Gpopulation[i];
}
CGpopulation* CGpopulation::clone()
{
	return new CGpopulation(*this);
}
//operator function
CGenome* CGpopulation::operator[](const int index)
{
	if(index>=m_pObjGaparameter->PopNum())
    {
        ERROR_OUTPUT("Error: index is out of range","CGpopulation","operator[]");
        exit(-1);
    }
	return m_Gpopulation[index];
}
double CGpopulation::operator[](std::string mystr)
{
	 std::string indexName = boost::to_lower_copy<std::string>(myStr);

	 if (indexName=="minraw")
	     return this->m_MinRawScore;
	 else if(indexName=="minfit")
	     return this->m_MinFitness;
	 else if(indexName=="maxfit")
	     return this->m_MaxFitness;
	 else if(indexName=="maxraw")
	     return this->m_MaxRawScore;
	 else if(indexName=="avgraw")
	     return this->m_AvgRawScore;
	 else if(indexName=="avgfit")
	     return this->m_AvgFitness;
     else if(indexName=="rawdev")
         return this->m_DevrawScore;
     else if(indexName=="minori")
         return this->m_MinOriValue;
     else if(indexName=="maxori")
         return this->m_MaxOriValue;
	 else{
	 	 ERROR_OUTPUT("Error: indexStr cannot match with them","CGpopulation","operator[]");
		 abort();
	 }
}

// property function
//void CGpopulation::init()
//{
//
//}
void CGpopulation::asscendSort()
{
	std::sort(m_Gpopulation.begin(),m_Gpopulation.end(),\
      [](CGenome* A,CGenome* B){ return (*A)["fitness"] < (*B)["fitness"];});   // use lamba
}
void CGpopulation::descendSort()
{
	std::sort(m_Gpopulation.begin(),m_Gpopulation.end(),\
      [](CGenome* A,CGenome* B){ return (*A)["fitness"] > (*B)["fitness"];});   // use lamba
}
//after genetic operatoration, updatepopulation is necessary.
void CGpopulation::updatePopulation()
{
	std::vector <CGenome*>::iterator it;
	for(it=this->m_Gpopulation.begin();it<this->m_Gpopulation.end();it++)
		(*it)->updateDecValueGene();
}
void CGpopulation::fitness_statistic()
{
    assert(m_Gpopulation.size()>0);

    std::vector <CGenome*>::iterator it;
	double sumFit=0.0,tmp=0.0;

    for(it=m_Gpopulation.begin();it<m_Gpopulation.end();it++)
        sumFit = sumFit + (*it)->fitness();

    this->m_AvgFitness  = sumFit/this->popNum();

	for(it=m_Gpopulation.begin();it<m_Gpopulation.end();it++)
	{
		(*it)->setRelativefitness( (*it)->fitness() / sumFit );
		tmp = tmp + (*it)->relativefitness();
		(*it)->setCumufitness(tmp);
	}
	this->m_pMaxFitGenome = std::max_element(m_Gpopulation.begin(),m_Gpopulation.end(),\
                                    [](CGenome* A,CGenome* B){ return (*A)["fitness"] < (*B)["fitness"];});
    this->m_pMinFitGenome = std::min_element(m_Gpopulation.begin(),m_Gpopulation.end(),\
                                    [](CGenome* A,CGenome* B){ return (*A)["fitness"] < (*B)["fitness"];});
    this->m_MinFitness = this->m_pMinFitGenome->fitness();
    this->m_MaxFitness = this->m_pMaxFitGenome->fitness();
}
void CGpopulation::raw_statistic()
{
	assert(m_Gpopulation.size()>0);

	std::vector <CGenome*>::iterator it;
    double sumRawScore=0.0,s,tmp=0;

	for(it=m_Gpopulation.begin();it<m_Gpopulation.end();it++)
		sumRawScore = sumRawScore +(*it)->rawscore();

    this->m_AvgRawScore = sumRawScore/this->popNum();

	for(it=m_Gpopulation.begin();it<m_Gpopulation.end();it++)
	{
		s=(*it)->rawscore() - this->avgRawScore;
		tmp = tmp + s*s;
	}
    this->m_DevrawScore = sqrt(tmp);

    this->m_pMaxRawGenome = std::max_element(m_Gpopulation.begin(),m_Gpopulation.end(),\
                                    [](CGenome* A,CGenome* B){ return (*A)["rawscore"] < (*B)["rawscore"];});
    this->m_pMinRawGenome = std::min_element(m_Gpopulation.begin(),m_Gpopulation.end(),\
                                    [](CGenome* A,CGenome* B){ return (*A)["rawscore"] < (*B)["rawscore"];});
    this->m_MinRawScore = this->m_pMinRawGenome->rawscore();
    this->m_MaxRawScore = this->m_pMaxRawGenome->rawscore();
	//obtained best fitness

	this->m_pMinOriGenome = std::min_element(m_Gpopulation.begin(),m_Gpopulation.end(),\
                                    [](CGenome* A,CGenome* B){ return (*A)["origvalue"] < (*B)["origvalue"];});
    this->m_pMaxOriGenome = std::max_element(m_Gpopulation.begin(),m_Gpopulation.end(),\
                                    [](CGenome* A,CGenome* B){ return (*A)["origvalue"] < (*B)["origvalue"];});
    this->m_MinOriValue = this->m_pMinOriGenome->origValue();
    this->m_MaxOriValue = this->m_pMaxOriGenome->origValue();
}
//input output
int CGpopulation::popNum()
{
	return m_pObjGaparameter->PopNum();
}
void CGpopulation::setPopNum(int popnum)
{
    if(popnum<=0)
    {
        ERROR_OUTPUT("Error: pop num <=0","CGpopulation","setPopNum");
        exit(-1);
    }
    (*m_pObjGaparameter)["Population"]=std::to_string(popnum);

}
std::vector <GENEVAR>* CGpopulation::geneVarArray()
{
    return m_pObjGaparameter->GeneVar();
}
void CGpopulation::setGeneVarArray(std::vector <GENEVAR>* my_geneVar)
{
    if(my_geneVar==nullptr)
    {
        ERROR_OUTPUT("Error: Pointer of my_geneVar is null","CGpopulation","setGeneVarArray");
        exit(-1);
    }
    m_pObjGaparameter->setGeneVar(my_geneVar);
}
void CGpopulation::modifyPopulation(std::vector <CGenome*> *newGenome)
{
    assert(newGenome);
}


}


