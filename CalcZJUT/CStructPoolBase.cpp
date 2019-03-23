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
#include "../IACS.h"
#include "CStructPoolBase.h"
#include "CModelBase.h"
#include "../Util/log.hpp"

using util::Log;

namespace CALCZJUT{

CStructPoolBase::CStructPoolBase(CParameter* othr)
:m_pParameter(othr)
{
    //ctor
}

CStructPoolBase::~CStructPoolBase()
{
     for(size_t i=0;i<m_CalcStructPool.size();i++)
        if(m_CalcStructPool[i]!=nullptr)
           delete m_CalcStructPool[i];
     m_CalcStructPool.clear();
}
CModelBase* CStructPoolBase::operator[](size_t index)
{
    if(index>=this->m_CalcStructPool.size())
    {
         Log::Error<<"Index is out of range! CStructPoolBase::operator[]!\n";
         boost::throw_exception(std::runtime_error("Index is out of range! CStructPoolBase::operator[]!"));
    }
    return this->m_CalcStructPool[index];
}
void CStructPoolBase::VarRangeStructRange(std::vector<VarRangeStruct>& mth)
{

}
void CStructPoolBase::init()
{

}




}
