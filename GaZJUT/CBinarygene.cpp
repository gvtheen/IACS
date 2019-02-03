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
#include<math.h>
#include "CBinarygene.h"
#include "../Util/CRandomgenerator.h"
#include "../Util/log.hpp"
#include "../Util/utilFunction.h"
namespace GAZJUT{

CBinarygene::CBinarygene():CGenebase()
{
    this->codeType=BINARY;
}
CBinarygene::CBinarygene(GeneVAR *myVal):CGenebase(myVal)
{//ctor
//    CGenebase::CGenebase(myVal);
    this->init(myVal);
    this->codeType=BINARY;
}

CBinarygene::~CBinarygene()
{
    m_bitdata.clear();
}
void CBinarygene::updatecode(double value)
{
    double low,high;
	int i, m;

	long int TenValue;
	low=this->m_GeneVAR->min;
    high=this->m_GeneVAR->max;
    if (value<low || value>high)
    {
    	util::Log::Error<<"The value exceeds the variable range!"<<"ref: CGene_encode"<<std::endl;
        boost::throw_exception(std::runtime_error("The value exceeds the variable range! CBinarygene::updatecode!\n"));
	}
	m=this->m_bitdata.size();
	TenValue=(int)((value-low)*(std::pow(2.0,m)-1)/(high-low));
	//clear the value of m_bitdata.
	this->m_bitdata.set(0);
	for(i=m-1;i>=0;i--)
	{
	    this->m_bitdata[i]=TenValue%2;
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
   //  std::cout<<(CGenebase::m_GeneVAR)->max<<endl;
     low= this->m_GeneVAR->min;
     high=this->m_GeneVAR->max;
     sum=0.0;
     m=this->m_bitdata.size();
   //  std::cout<<m<<endl;
     for(i=0;i<m;i++)
          sum=sum+(double)(this->m_bitdata[i])*std::pow(2.0,m-i-1);

     this->m_value=low+(high-low)*sum/(std::pow(2.0,m)-1);
    //std::cout<<"value of CBinarygene:"<<this->m_value<<endl;
     return  this->m_value;
}
void CBinarygene::init(GeneVAR*myVal)
{
    double rndNum;
    assert(myVal);

    this->m_GeneVAR=myVal;
   // using CGeneBase::init(myVal);
    this->m_bitNum = util::calcBitNum(*myVal);

    util::CRandomgenerator *rndgenerator=new util::CRandomgenerator();

    this->m_bitdata.clear();
    for(size_t i=0;i<this->m_bitNum;i++)
    {
        rndNum=rndgenerator->uniformRandom01(i+100);
        if (rndNum>=0.5)
	       this->m_bitdata.push_back(1);
	    else
	       this->m_bitdata.push_back(0);
    }
    this->m_value=decode();

    delete rndgenerator;
}
Bitset& CBinarygene::bitGene()
{
    return this->m_bitdata;
}
size_t CBinarygene::bitNum()
{
    return this->m_bitNum;
}





}
