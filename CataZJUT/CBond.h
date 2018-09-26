#ifndef CBOND_H
#define CBOND_H
#include "Point-Vector.h"
#include "CElement.h"

namespace CATAZJUT{

class CConfigurationBase;
class CAtom;

class CBond
{
    public:
        CBond(CConfigurationBase*,size_t);
        virtual ~CBond();

        Point3 Center() const;
        double Length() const;
        CAtom* atom(size_t index) const;
        CAtom* atom1() const;
        CAtom* atom2() const;
        CAtom* otherAtom(const CAtom *atom) const;

        bool contains(const CAtom *atom) const;
        bool contains(const CElement &element) const;
        bool containsBoth(const CAtom *a, const CAtom *b) const;
        bool containsBoth(const CElement &a, const CElement &b) const;
        bool isTerminal() const;

        inline CConfigurationBase* Configuration() const;
        inline size_t index() const;
    private:
        friend class CConfigurationBase;
    private:
        CConfigurationBase *m_pConfiguration;
        size_t m_index;

};

}

#include "CBond_sub.h"
#endif // CBOND_H
