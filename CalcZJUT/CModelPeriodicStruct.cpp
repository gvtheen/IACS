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
#include "CModelPeriodicStruct.h"


namespace CALCZJUT{

CModelPeriodicStruct::CModelPeriodicStruct(CParameter* mPara,
                                           size_t index)
:CModelBase(mPara,index)
{
     this->m_pPeriodicFramework = new CATAZJUT::CPeriodicFramework(mPara);
}

CModelPeriodicStruct::~CModelPeriodicStruct()
{
    delete this->m_pPeriodicFramework;
}

CModelBase* CModelPeriodicStruct::clone()
{
    CModelPeriodicStruct* res =  new CModelPeriodicStruct (this->m_pParameter,0);

    res->m_pGeneVAR->assign(this->m_pGeneVAR->begin(),this->m_pGeneVAR->end());
    res->m_IsNeedRandomInit = this->RandomInitState();
    res->setChemicalFormula(this->chemicalFormula());

    return res;
}
void CModelPeriodicStruct::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{

}
void CModelPeriodicStruct::getGeneValuefromStruct(std::vector<double>& currentGeneRealValue)
{

}
void CModelPeriodicStruct::GeneVARRange(std::vector<GeneVAR>& currentGeneVARible)
{

}
void CModelPeriodicStruct::eliminateCloseContacts(CATAZJUT::CPeriodicFramework* strut,
                                                  double distanceCutOff=1.0)
{


}
void eliminateFragment(CATAZJUT::CPeriodicFramework*);


}
