#include<math.h>
#include "CBinarygene.h"
#include "../Util/CRandomgenerator.h"
#include "../Util/log.hpp"
#include "GaUtilityFunction.h"
namespace GAZJUT{

CBinarygene::CBinarygene():CGenebase()
{
     this->m_bitdata=nullptr;
}
CBinarygene::CBinarygene(GENEVAR myVal):CGenebase(myVal)
{//ctor
//    CGenebase::CGenebase(myVal);
    this->init(myVal);
}

CBinarygene::~CBinarygene()
{
    delete m_bitdata;
}
void CBinarygene::updatecode(double value)
{
    double low,high;
	int i, m;

	long int TenValue;
	low=this->m_geneVar->min;
    high=this->m_geneVar->max;
    if (value<low || value>high)
    {
    	util::Log::OutputToFile<<"ERROR: The value exceeds the variable range!"<<"ref: CGene_encode"<<std::endl;
    	abort();
	}
	m=this->m_bitdata->size();
	TenValue=(int)((value-low)*(std::pow(2.0,m)-1)/(high-low));
	//clear the value of m_bitdata.
	std::fill(m_bitdata->begin(),m_bitdata->end(),0);
	for(i=m-1;i>=0;i--)
	{
	    this->m_bitdata->at(i)=TenValue%2;
        TenValue=TenValue/2;
        if (TenValue==0) break;
	}
//	for(i=0;i<m;i++)
//        cout<<m_bitdata->at(i);
//    cout<<endl;
	this->decode();
}
double CBinarygene::decode()
{
     int i,m;
     double low,high,sum;
   //  std::cout<<"value of CBinarygene:"<<endl;
   //  std::cout<<(CGenebase::m_geneVar)->max<<endl;
     low= this->m_geneVar->min;
     high=this->m_geneVar->max;
     sum=0.0;
     m=this->m_bitdata->size();
   //  std::cout<<m<<endl;
     for(i=0;i<m;i++)
          sum=sum+(double)(this->m_bitdata->at(i))*std::pow(2.0,m-i-1);

     this->m_value=low+(high-low)*sum/(std::pow(2.0,m)-1);
    //std::cout<<"value of CBinarygene:"<<this->m_value<<endl;
     return  this->m_value;
}
void CBinarygene::init(GENEVAR myVal)
{
    double rndNum;
    assert(this->m_geneVar);

    if(this->m_geneVar==nullptr)
    {
        this->m_geneVar=new GENEVAR();
        *m_geneVar=myVal;
    }
   // using CGeneBase::init(myVal);
    this->m_bitNum = calcBitNum(myVal);
    this->m_bitdata= new (std::vector<unsigned int>);
    util::CRandomgenerator *rndgenerator=new util::CRandomgenerator();

    for(int i=0;i<this->m_bitNum;i++)
    {
        rndNum=rndgenerator->uniformRandom01(i+100);
        if (rndNum>=0.5)
	       this->m_bitdata->push_back(1);
	    else
	       this->m_bitdata->push_back(0);
    }
    this->m_value=decode();
}
std::vector<unsigned int>* CBinarygene::bitGene()
{
    return this->m_bitdata;
}
int CBinarygene::bitNum()
{
    return this->m_bitNum;
}





}
