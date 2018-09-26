#ifndef CBOND_SUB_H
#define CBOND_SUB_H

#include "CBond.h"
#include "CConfigurationBase.h"

namespace CATAZJUT{

inline CConfigurationBase* CBond::Configuration() const
{
    return m_pConfiguration;
}

/// Returns the bond's index in the molecule.
inline size_t CBond::index() const
{
    return m_index;
}


}

#endif
