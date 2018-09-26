#include "CGraygene.h"
#include "GaUtilityFunction.h"
namespace GAZJUT{
CGraygene::CGraygene():CGenebase()
{
    this->m_bitdata=nullptr;
}
CGraygene::CGraygene(GENEVAR myVal):CGenebase(myVal)
{
    this->init(myVal);
}

CGraygene::~CGraygene()
{
    delete m_bitdata;
}

double CGraygene::decode()
{
     int i,m;
     double low,high,sum;

     low= this->m_geneVar->min;
     high=this->m_geneVar->max;
     sum=0.0;
     m=this->m_bitdata->size();

     std::vector<unsigned int>  *tempData=new (std::vector<unsigned int>)(m);
     tempData->assign(this->m_bitdata->begin(),this->m_bitdata->end());
     grayTobit(tempData);

     for(i=0;i<m;i++)
          sum=sum+(double)(tempData->at(i))*pow(2.0,m-i-1);

     this->m_value=low+(high-low)*sum/(pow(2.0,m)-1);
     return  this->m_value;
}
void CGraygene::init(GENEVAR myVal)
{
    double rndNum;
    assert(this->m_geneVar);
    if(this->m_geneVar==nullptr)
    {
        this->m_geneVar=new GENEVAR();
        *m_geneVar=myVal;
    }
    this->m_bitNum = calcBitNum(myVal);
    this->m_bitdata=new (std::vector<unsigned int>)(this->m_bitNum);

    CRandomgenerator *rndgenerator=new CRandomgenerator();
    for(int i=0;i<this->m_bitNum;i++)
    {
        rndNum=rndgenerator->uniformRandom01(i+100);
        if (rndNum>=0.5)
	       this->m_bitdata->at(i)=1;
	    else
	       this->m_bitdata->at(i)=0;
    }
    this->m_value=decode();
}
void CGraygene::updatecode(double value)
{
    double low,high;
	int i, m;

	long int TenValue;
	low=this->m_geneVar->min;
    high=this->m_geneVar->max;
    if (value<low || value>high)
    {
    	ERROR_OUTPUT("Error:The value exceeds the valable range!","CGene","encode");
    	abort();
	}
	m=this->m_bitdata->size();
	TenValue=(int)((value-low)*(pow(2.0,m)-1)/(high-low));
	//clear the value of m_bitdata;
	std::fill(m_bitdata->begin(),m_bitdata->end(),0);
	for(i=m-1;i>=0;i--)
	{
	    this->m_bitdata->at(i)=TenValue%2;
        TenValue=TenValue/2;
        if (TenValue==0) break;
	}
    bitTogray(this->m_bitdata);
    this->m_value=decode();
}
std::vector<unsigned int>* CGraygene::bitGene()
{
    return this->m_bitdata;
}
int CGraygene::bitNum()
{
    return this->m_bitNum;
}
}
