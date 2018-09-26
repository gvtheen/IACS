#ifndef CFRAGMENT_H
#define CFRAGMENT_H

//#include "gacata.h"
#include <vector>

#include "../Util/Bitset.h"

using util::Bitset;

namespace CATAZJUT{

class CConfigurationBase;
class CAtom;
class CBond;

class CFragment
{
    public:

        CFragment(CConfigurationBase*,const Bitset&);
        virtual ~CFragment();

            // properties
       inline size_t size() const;
       inline CConfigurationBase* Configuration() const;

    // structure
       inline CAtom* atom(size_t index) const;
       inline std::vector<CAtom *> atoms() const;
       inline size_t atomCount() const;
       inline bool contains(const CAtom *atom) const;
       inline bool containMetal()const;
       Bitset bitSet()const;
       std::vector<CBond *> bonds() const;
       size_t bondCount() const;
       bool contains(const CBond *bond) const;

    protected:

    private:
        CConfigurationBase*  m_pConfiguration;
        Bitset  m_bitset;
};
}
#include "CFragment_sub.h"
#endif // CFRAGMENT_H
