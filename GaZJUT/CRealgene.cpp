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
#include "CRealgene.h"
#include "../Util/CRandomgenerator.h"

namespace GAZJUT{

CRealgene::CRealgene()
:CGenebase()
{
    this->codeType=REAL;
}
CRealgene::CRealgene(VarRangeStruct*myVal)
:CGenebase(myVal)
{
    this->init(myVal);
    this->codeType=REAL;
}

CRealgene::~CRealgene()
{
    //dtor
}
double CRealgene::decode()
{
    return this->m_value;
}
void CRealgene::init(VarRangeStruct *myVal)
{
    double low,high;
    assert(myVal);

    this->m_VarRangeStruct=myVal;

    util::CRandomgenerator *rndgenerator=new util::CRandomgenerator();
    low =this->m_VarRangeStruct->min;
    high=this->m_VarRangeStruct->max;
    this->m_value=low + (high-low)*rndgenerator->uniformRandom01(100);
}
void CRealgene::updatecode(double value)
{
    this->m_value=value;
}
double& CRealgene::realGene()
{
    return this->m_value;
}
size_t CRealgene::bitNum()
{
    return 1;
}

}
