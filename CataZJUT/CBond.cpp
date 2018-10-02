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
void CBond::setOrder(CBond::BondOrderType order)
{
    m_pConfiguration->m_pData->bondOrders[this->m_index]=order;
}
CBond::BondOrderType CBond::order() const
{
    return m_pConfiguration->m_pData->bondOrders[this->m_index];
}
bool is(CBond::BondOrderType order) const
{
    if(m_pConfiguration->m_pData->bondOrders[this->m_index] == order)
        return true;
    else
        return false;
}

}
