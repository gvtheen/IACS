#ifndef CFRAGMENT_SUB_H
#define CFRAGMENT_SUB_H
#include "CFragment.h"
#include "CAtom.h"
#include "foreach.h"
#include "Cbond.h"
#include "CConfigurationBase.h"
#include "../Util/log.hpp"
using util::Log;

namespace CATAZJUT{

inline size_t CFragment::size() const
{
    return m_bitset.count();
}
inline CConfigurationBase* CFragment::Configuration() const
{
    return this->m_pConfiguration;
}
inline CAtom* CFragment::atom(size_t index) const
{
     size_t pos = m_bitset.find_first();
     while(index--)
     {
         pos = m_bitset.find_next(pos);
         if( pos == Bitset::npos )
            break;
     }
     return m_pConfiguration->m_Atom[pos];
}
inline std::vector<CAtom*> CFragment::atoms() const
{
      std::vector<CAtom*> atoms;
      #ifdef DEBUG
           Log::Debug<<"CFragment::atoms()size:"<<m_bitset.size()<< std::endl;
      #endif
      for(size_t i=0;i<m_bitset.size();i++)
        if(m_bitset.test(i)==1)   // judge whether m_bitset[i]==1
            atoms.push_back(m_pConfiguration->m_Atom[i]);
      return atoms;
}
inline size_t CFragment::atomCount() const
{
     return m_bitset.count();
}
inline bool CFragment::contains(const CAtom *atom) const
{
     return m_bitset.test(atom->index());
}
inline bool CFragment::containMetal()const
{
    foreach(CAtom* atom, atoms())
    {
        if(atom->element().isMetal()==true)
            return true;
    }
    return false;
}
inline std::vector<std::pair<std::string,size_t>>* CFragment::composition()
{
       return nullptr;
}
inline std::string CFragment::formula()
{
     return nullptr;
}


}
#endif

