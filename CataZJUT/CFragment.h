#ifndef CFRAGMENT_H
#define CFRAGMENT_H

//#include "gacata.h"
#include <vector>

#include "../Util/Bitset.h"
#include "../Util/Point-Vector.h"

using util::Bitset;
using util::Point3;
using util::Vector3;

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

       inline std::vector<std::pair<std::string,size_t>>* composition();
       inline std::string formula();

       Bitset bitSet()const;
       std::vector<CBond *> bonds() const;
       size_t bondCount() const;
       bool contains(const CBond *bond) const;
       bool isBondTo(const CFragment* otherFragment);
       std::vector<Point3> coordinates() const;
       void move(const Vector3&);

    protected:

    private:
        CConfigurationBase*  m_pConfiguration;
        Bitset  m_bitset;
};

}
#include "CFragment_sub.h"
#endif // CFRAGMENT_H
