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
#include "CGraygene.h"
#include "../Util/CRandomgenerator.h"
#include "../Util/utilFunction.h"
#include "../Util/log.hpp"

namespace GAZJUT{
CGraygene::CGraygene():CGenebase()
{
    //nothing is done.
    this->codeType=GRAY;
}
CGraygene::CGraygene(VarRangeStruct*myVal):CGenebase(myVal)
{
    this->init(myVal);
    this->codeType=GRAY;
}

CGraygene::~CGraygene()
{
    m_bitdata.clear();
    //nothing is done.
}

double CGraygene::decode()
{
     double low,high,sum;

     low= this->m_VarRangeStruct->min;
     high=this->m_VarRangeStruct->max;
     sum=0.0;
     size_t m=this->m_bitdata.size();

     Bitset tempData=this->m_bitdata;
     util::grayTobit(tempData);

     for(size_t i=0;i<m;i++)
          sum=sum+(double)(tempData[i])*pow(2.0,m-i-1);

     this->m_value=low+(high-low)*sum/(pow(2.0,m)-1);
     tempData.clear();
     return  this->m_value;
}
void CGraygene::init(VarRangeStruct* myVal)
{
    double rndNum;
    assert(myVal);

    this->m_VarRangeStruct=myVal;

    this->m_bitNum = util::calcBitNum(*myVal);

    util::CRandomgenerator *rndgenerator=new util::CRandomgenerator();
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
void CGraygene::updatecode(double value)
{
    double low,high;
	int i, m;

	long int TenValue;
	low=this->m_VarRangeStruct->min;
    high=this->m_VarRangeStruct->max;
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
    util::bitTogray(this->m_bitdata);
    this->m_value=decode();
}
Bitset& CGraygene::bitGene()
{
    return this->m_bitdata;
}
size_t CGraygene::bitNum()
{
    return this->m_bitNum;
}
}
