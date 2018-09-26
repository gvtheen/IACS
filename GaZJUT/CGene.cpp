#include "CGene.h"
#include <math.h>
#include <stddef.h>

namespace GAZJUT{
CGene::CGene()
{
	GENEVAR iniVal={0,1.0,0.01};
	this->m_geneVar=iniVal;
	this->m_gene=NULL;
	this->m_bitNum=0;
	this->m_value=0;
}
CGene::CGene(GENEVAR myval)
{
	this->init(myval);
}
CGene::CGene(const CGene& oldCGene)
{
	this->m_geneVar=oldCGene.m_geneVar;
    this->m_bitNum=oldCGene.m_bitNum;
	this->m_gene=new (std::vector<unsigned int>);
	for(int i=0;i<m_bitNum;i++)
	   this->m_gene->push_back(oldCGene.m_gene->at(i));
	this->m_value=oldCGene.m_value;
}
CGene* CGene::clone()
{
	return new CGene(*this);
}
void CGene::randomInitial()
{
	for(int i;i<this->m_bitNum;i++)
	    if (genrand_real1()>=0.5)
	       this->m_gene->at(i)=1;
	    else
	       this->m_gene->at(i)=0;
}
int CGene::encode(double value)
{
	double low,high;
	int i, m,start;

	long int TenValue;
	low=this->m_geneVar.min;
    high=this->m_geneVar.max;
    if (value<low || value>high)
    {
    	ERROR_OUTPUT("Error:The value exceeds the valable range!","CGene","encode");
    	abort();
	}
	m=this->BinaryBitNum();
	TenValue=(int)((value-low)*(pow(2.0,m)-1)/(high-low));
	for(i=m-1;i>=0;i--)
	{
	    this->m_gene->at(i)=TenValue%2;
        TenValue=TenValue/2;
        if (TenValue==0) break;
	}
}
double CGene::decode()
{
	 int i,k,m;
     double low,high,sum;

     low= this->m_geneVar.min;
     high=this->m_geneVar.max;
     sum=0.0;
     m=this->BinaryBitNum();
     for(i=0;i<m;i++)
          sum=sum+(double)(this->m_gene->at(i))*pow(2.0,m-i-1);

     this->m_value=low+(high-low)*sum/(pow(2.0,m)-1);

    return  this->m_value;
}

double CGene::value()
{
	return this->m_value;
}


int CGene::BinaryBitNum()
{
    double c,m;
    m=(this->m_geneVar.max - this->m_geneVar.min)/this->m_geneVar.accuracy;
    c=log10(m)/log10(2.0);
    //printf("%f\n",c);
    return (int)c+1;
}
void CGene::init(GENEVAR myval)
{
	this->m_geneVar=myVal;
	this->m_bitNum=this->BinaryBitNum();
	this->m_gene=new (vector <int>)(this->m_bitNum);
}
void CGene::updatebit(vector <unsigned int>* mygene)
{
	if(mygene->size()!=this->m_gene->size())
	{
    	ERROR_OUTPUT("Error:The mygene bitsize exceeds this of m_gene!","CGene","updatebit");
    	exit(0);
	}
	this->m_gene->assign(mygene->begin(),mygene->end());
	this->decode();
}
void CGene::init()
{

}
CGene::~CGene()
{
	delete m_gene;
}


} ///namespace
