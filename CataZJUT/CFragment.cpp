#include "CFragment.h"
#include <cassert>
#include "CAtom.h"
#include "foreach.h"
#include "Cbond.h"
#include "CConfigurationBase.h"

namespace CATAZJUT{

CFragment::CFragment(CConfigurationBase* mConfig,const Bitset& bit_set)
:m_pConfiguration(mConfig),m_bitset(bit_set)
{
    //ctor
}
CFragment::~CFragment()
{
    //dtor
}
std::vector<CBond*> CFragment::bonds() const
{
    std::vector<CBond*> bonds;
    foreach(CAtom* atom, atoms()){
        foreach(CBond* bond, atom->bonds()){
           if(std::find(bonds.begin(),bonds.end(),bond)==bonds.end())
               bonds.push_back(bond);
        }
    }
    return bonds;
}
size_t CFragment::bondCount() const
{
   return bonds().size();
}
bool CFragment::contains(const CBond *bond) const
{
   return bond->atom1()->fragment()==this;
}
Bitset CFragment::bitSet()const
{
    return this->m_bitset;
}
}
