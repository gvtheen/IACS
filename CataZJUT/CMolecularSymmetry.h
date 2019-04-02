#ifndef CMOLECULARSYMMETRY_H
#define CMOLECULARSYMMETRY_H

#include "../Util/Bitset.h"

using util::Bitset;

namespace CATAZJUT{

class CConfigurationBase;

class CMolecularSymmetry
{
    public:
        CMolecularSymmetry(CConfigurationBase*);
        CMolecularSymmetry(CConfigurationBase*,Bitset);
        virtual ~CMolecularSymmetry();

        char* GetPointGroup();

    private:
        CConfigurationBase    *m_pConfiguration;
        char* point_group;
        Bitset      m_AtomicBits;
};


}
#endif // CMOLECULARSYMMETRY_H
