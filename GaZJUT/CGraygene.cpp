#include "CGraygene.h"
#include "../Util/CRandomgenerator.h"
#include "GaUtilityFunction.h"
#include "../Util/log.hpp"
namespace GAZJUT{
CGraygene::CGraygene():CGenebase()
{
    //nothing is done.
    this->codeType=GRAY;
}
CGraygene::CGraygene(GeneVAR myVal):CGenebase(myVal)
{
    this->init(myVal);
    this->codeType=GRAY;
}

CGraygene::~CGraygene()
{
    //nothing is done.
}

double CGraygene::decode()
{
     int i,m;
     double low,high,sum;

     low= this->m_GeneVAR->min;
     high=this->m_GeneVAR->max;
     sum=0.0;
     m=this->m_bitdata.size();

     Bitset tempData=this->m_bitdata;
     grayTobit(tempData);

     for(i=0;i<m;i++)
          sum=sum+(double)(tempData[i])*pow(2.0,m-i-1);

     this->m_value=low+(high-low)*sum/(pow(2.0,m)-1);
     tempData.clear();
     return  this->m_value;
}
void CGraygene::init(GeneVAR myVal)
{
    double rndNum;
    assert(this->m_GeneVAR);
    if(this->m_GeneVAR==nullptr)
    {
        this->m_GeneVAR=new GeneVAR();
        *m_GeneVAR=myVal;
    }
    this->m_bitNum = calcBitNum(myVal);

    util::CRandomgenerator *rndgenerator=new util::CRandomgenerator();
    for(int i=0;i<this->m_bitNum;i++)
    {
        rndNum=rndgenerator->uniformRandom01(i+100);
        if (rndNum>=0.5)
	       this->m_bitdata.push_back(1);
	    else
	       this->m_bitdata.push_back(0);
    }
    this->m_value=decode();
}
void CGraygene::updatecode(double value)
{
    double low,high;
	int i, m;

	long int TenValue;
	low=this->m_GeneVAR->min;
    high=this->m_GeneVAR->max;
    if (value<low || value>high)
    {
    	util::Log::Error<<"The value exceeds the valable range! CGraygene::updatecode!\n";
    	boost::throw_exception(std::runtime_error("The value exceeds the variable range! CGraygene::updatecode!\n"));
	}
	m=this->m_bitdata.size();
	TenValue=(int)((value-low)*(pow(2.0,m)-1)/(high-low));
	//clear the value of m_bitdata;
	this->m_bitdata.set(0);
	for(i=m-1;i>=0;i--)
	{
	    this->m_bitdata[i]=TenValue%2;
        TenValue=TenValue/2;
        if (TenValue==0) break;
	}
    bitTogray(this->m_bitdata);
    this->m_value=decode();
}
Bitset& CGraygene::bitGene()
{
    return this->m_bitdata;
}
int CGraygene::bitNum()
{
    return this->m_bitNum;
}
}
