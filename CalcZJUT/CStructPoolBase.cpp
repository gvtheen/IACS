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
#include "CCalcStructBasePool.h"



namespace CALCZJUT{

CStructPoolBase::CStructPoolBase(CParameter* othr)
:m_Parameter(othr)
{
    //ctor
}

CStructPoolBase::~CStructPoolBase()
{
    //dtor
}
CModelBase* CStructPoolBase::operator[](size index)
{
    if(index>=this->m_CalcStructPool.size())
    {
         Log::Error<<"Index is out of range! CStructPoolBase::operator[]!\n";
         boost::throw_exception(std::runtime_error("Index is out of range! CStructPoolBase::operator[]!"));
    }
    return this->m_CalcStructPool[index];
}
void CStructPoolBase::GeneVARRange(std::vector<GeneVAR>& mth)
{

}
void CStructPoolBase::init()
{

}




}
