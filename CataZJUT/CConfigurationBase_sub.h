#ifndef CCONFIGURATIONBASE_SUB_H
#define CCONFIGURATIONBASE_SUB_H
#include <boost/foreach.hpp>
#include "CConfigurationBase.h"
#include "../GaZJUT/GaUtilityFunction.h"
#include "../Util/log.hpp"
#include "CElement.h"

namespace CATAZJUT{

inline size_t CConfigurationBase::size() const
{
    return atomCount();
}
inline bool CConfigurationBase::isEmpty() const
{
    return size()==0;
}
inline CConfigurationBase::AtomRange CConfigurationBase::atoms() const
{
    return boost::make_iterator_range(m_Atom.begin(),m_Atom.end());
}
inline size_t CConfigurationBase::atomCount() const
{
    return m_Atom.size();
}
inline CAtom* CConfigurationBase::atom(size_t index) const
{
    return m_Atom[index];
}
//
template<typename Range>
inline void CConfigurationBase::removeAtoms(Range range)
{
   this->removeAtoms(std::vector<CAtom*>(range.begin(),range.end()));
}
//
template<typename Range>
inline void CConfigurationBase::removeBonds(Range range)
{
   this->removeBonds(std::vector<CBond *>(range.begin(), range.end()));
}
//
inline void CConfigurationBase::checkElement(CElement& currentElement)
{
   if(!currentElement.isValid())
   {
       util::Log::Error<<"Atomic label is error! CConfigurationBase_checkElement"<<std::endl;
       boost::throw_exception(std::runtime_error("Atomic label is error!Check the file: Error_information.txt."));
   }
}

}

#endif // CCONFIGURATIONBASE_SUB_H
