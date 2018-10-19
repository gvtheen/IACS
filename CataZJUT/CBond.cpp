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
#include "CBond.h"
#include "CAtom.h"
#include "CConfigurationBase.h"
#include "CConfigurationPrivateData.h"

namespace CATAZJUT{

CBond::CBond(CConfigurationBase* mConfig,size_t index)
 :m_pConfiguration(mConfig),m_index(index)
{

}
CBond::~CBond()
{
    //dtor
}
CAtom* CBond::atom(size_t index) const
{
   return index == 0 ? atom1() : atom2();
}
CAtom* CBond::atom1() const
{
   return m_pConfiguration->m_pData->bondAtoms[this->m_index].first;
}
CAtom* CBond::atom2() const
{
   return m_pConfiguration->m_pData->bondAtoms[this->m_index].second;
}
CAtom* CBond::otherAtom(const CAtom *atom) const
{
    const std::pair<CAtom*,CAtom*> &mpair = m_pConfiguration->m_pData->bondAtoms[this->m_index];
    return mpair.first == atom? mpair.second: mpair.first;
}
bool CBond::contains(const CAtom *atom) const
{
    return atom1()== atom || atom2()==atom;
}
bool CBond::contains(const CElement &element) const
{
    return atom1()->IsElement(element) || atom2()->IsElement(element);
}
bool CBond::containsBoth(const CAtom *a, const CAtom *b) const
{
    return contains(a) && contains(b);
}
bool CBond::containsBoth(const CElement &a, const CElement &b) const
{
    return (atom1()->IsElement(a)&& atom2()->IsElement(b)) ||(atom1()->IsElement(b)&& atom2()->IsElement(a)) ;
}
bool CBond::isTerminal() const
{
    return atom1()->isTerminal()&& atom2()->isTerminal();
}
void CBond::setBondOrder(CBond::BondOrderType order)
{
    m_pConfiguration->m_pData->bondOrders[this->m_index]=order;
}
CBond::BondOrderType CBond::bondOrder() const
{
    return m_pConfiguration->m_pData->bondOrders[this->m_index];
}
bool CBond::is(CBond::BondOrderType order) const
{
    if(m_pConfiguration->m_pData->bondOrders[this->m_index] == order)
        return true;
    else
        return false;
}
//other atom <------- atom
Vector3 CBond::bondVector(const CAtom* atom)const
{
    if(this->contains(atom)==false)
        return Vector3(0,0,0);
    return otherAtom(atom)->position() - atom->position();
}



}
