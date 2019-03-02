#ifndef CBOND_H
#define CBOND_H
#include "../Util/Point-Vector.h"
#include "CElement.h"


using util::Point3;
using util::Vector3;

namespace CATAZJUT{

class CConfigurationBase;
class CAtom;

class CBond
{
    public:
        enum BondOrderType{
             Single = 1,
             Double = 2,
             Triple = 3,
             Quadruple = 4
             };
        CBond(CConfigurationBase*,size_t);
        virtual ~CBond();

        Point3 Center() const;
        double length() const;
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
        // bond order
        void setBondOrder(BondOrderType order);
        BondOrderType bondOrder() const;
        size_t bondOrderValue()const;
        bool is(BondOrderType order) const;
        //bond vector

        Vector3 bondVector(const CAtom*)const;

    private:
        friend class CConfigurationBase;
    private:
        CConfigurationBase *m_pConfiguration;
        size_t m_index;

};

}

#include "CBond_sub.h"
#endif // CBOND_H
