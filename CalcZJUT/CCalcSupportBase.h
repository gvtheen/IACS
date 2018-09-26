#ifndef CCALCSUPPORTBASE_H
#define CCALCSUPPORTBASE_H
#include <vector>
#include <boost/range/iterator_range.hpp>
//#
#include "../Util/Bitset.h"
//#include "CCartesianCoordinates.h"
using util::Bitset;

namespace CATAZJUT{
         class CConfigurationBase;
         class CAtom;
         class CCartesianCoordinates;
         class CPeriodicFramework;
}
namespace CALCZJUT{

class CCalcSupportBase
{
    public:
        typedef boost::iterator_range<std::vector<CATAZJUT::CAtom *>::const_iterator> AtomRange;

        //construction & deconstruction
        CCalcSupportBase(CATAZJUT::CPeriodicFramework* mpconf, Bitset& supportBit);
        CCalcSupportBase(CATAZJUT::CPeriodicFramework* mpconf, std::vector<CATAZJUT::CAtom*> support);

        virtual ~CCalcSupportBase();

        size_t size() const;
        bool isEmpty() const;

        CATAZJUT::CConfigurationBase* configuration() const;
        void setConfiguration(CATAZJUT::CPeriodicFramework*);
        CATAZJUT::CCartesianCoordinates* coordinates();

        //structure
        CATAZJUT::CAtom* atom(size_t index) const;
        Bitset bitSet()const;
        AtomRange atoms() const;
        size_t atomCount() const;

   public:
        CATAZJUT::CPeriodicFramework*      m_pPeriodicFramework;
        std::vector<CATAZJUT::CAtom*>      m_Atoms;
                     Bitset      m_AtomicBits;

};


}
#endif // CCalcSupportBase_H
