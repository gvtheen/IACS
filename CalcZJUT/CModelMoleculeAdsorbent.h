#ifndef CModelMoleculeAdsorbent_H
#define CModelMoleculeAdsorbent_H
/*
Noting:
   This class is used to treat the adsorbent on the support rather than pure cluster evolution.
*/
#include <vector>
#include <boost/range/iterator_range.hpp>
#include "../Util/Point-Vector.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../Util/Bitset.h"

namespace CATAZJUT{
 class CConfigurationBase;
 class CAtom;
}

using util::Point3;
using util::Bitset;

namespace CALCZJUT{

class CModelMoleculeAdsorbent
{
    public:
        typedef boost::iterator_range<std::vector<CATAZJUT::CAtom *>::const_iterator> AtomRange;

        CModelMoleculeAdsorbent(CATAZJUT::CConfigurationBase* mpConf,Bitset& molIndexBit);
        CModelMoleculeAdsorbent(CATAZJUT::CConfigurationBase* mpConf,std::vector<CATAZJUT::CAtom*>& molAtoms);
        virtual ~CModelMoleculeAdsorbent();

    size_t size() const;
    bool isEmpty() const;
    CATAZJUT::CConfigurationBase* configuration() const;
    void setConfiguration(CATAZJUT::CConfigurationBase*);
    //structure
    CATAZJUT::CAtom* atom(size_t index) const;
    AtomRange atoms() const;
    size_t atomCount() const;
    Bitset bitSet()const;

    //operators// only for the molecule;

    CModelMoleculeAdsorbent& operator=(const CModelMoleculeAdsorbent &moleculeadsorbent);

    Point3 gravityCentre();
    void moveBy(const Point3 vect);
    void moveBy(const double x,const double y,const double z);
    void rotate(const Point3& vect, const double& angle);

    private:
        CATAZJUT::CConfigurationBase*      m_pConfiguration;
        std::vector<CATAZJUT::CAtom*>      m_Atoms;
                               Bitset      m_AtomicBits;
};





}
#endif // CMOLECULE_H
