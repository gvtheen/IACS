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
#include <string>
#include <sstream>
#include "CModelBase.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "CIOCar.h"
namespace CALCZJUT{

CModelBase::CModelBase(CParameter* temParameter,size_t index)
:m_pParameter(temParameter),m_index(index)
{
    m_IsNeedRandomInit=true;
}

CModelBase::~CModelBase()
{
    //dtor
}
void CModelBase::init()
{

}
CATAZJUT::CPeriodicFramework* CModelBase::periodicFramework()
{
    return m_pPeriodicFramework;
}
void CModelBase::setPeriodicFramekwork(CATAZJUT::CPeriodicFramework* mbf)
{
    m_pPeriodicFramework = mbf;
}
void CModelBase::createSupport(const Bitset & mt)
{
    // no doing
}
void CModelBase::createMoleAdsorb(const Bitset & mt)
{
    // no doing
}
//void CModelBase::createStructureAtGene()
//{
//
//}
void CModelBase::setRandomInitState(const bool& mht)
{
    this->m_IsNeedRandomInit=mht;
}
bool CModelBase::RandomInitState()
{
    return this->m_IsNeedRandomInit;
}
//void CModelBase::removeStructureOfGene()
//{
//   for(size_t i=0;i<m_PopuPeriodicFramework.size();i++)
//       delete m_PopuPeriodicFramework[i];
//   m_PopuPeriodicFramework.clear();
//}

std::vector<std::pair<std::string,size_t>>& CModelBase::chemicalFormula()
{
     std::vector<std::pair<std::string,size_t>> none_res;
     return none_res;
}
void CModelBase::setChemicalFormula(const std::vector<std::pair<std::string,size_t>>& mth)
{

}

size_t CModelBase::index()
{
    return this->m_index;
}
void CModelBase::setIndex(size_t index)
{
     this->m_index=index;
}
void CModelBase::outputStructureToFile()
{
   CIOCar* io = new CIOCar(this->m_pPeriodicFramework);
   std::stringstream str ;
   str<<"Struct_Gene_" << m_pParameter->currentGenerationNum()<<"_Pop_" << this->index();
   std::string filename;
   str>>filename;
   io->output(filename);
   delete io;
}

}
