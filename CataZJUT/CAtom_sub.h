#ifndef CATOM_SUB_H
#define CATOM_SUB_H
#include "CAtom.h"
#include "CConfigurationBase.h"


namespace CATAZJUT{

inline CElement CAtom::element() const
{
     return m_pConfiguration->m_Element[m_index];
}
inline bool CAtom::IsElement(const CElement &element) const
{
    return element == this->element();
}
inline CConfigurationBase* CAtom::Configuration() const
{
    return this->m_pConfiguration;
}
inline size_t CAtom::index() const
{
    return this->m_index;
}
inline void CAtom::setIndex(const size_t index)
{
    this->m_index = index;
}

}

#endif
