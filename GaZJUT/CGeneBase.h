#ifndef CODE_H
#define CODE_H
#include<vector>
#include "GaDeclaration.h"
#include "../Util/Bitset.h"

using util::Bitset;

namespace GAZJUT{

class CGenebase
{
    public:
        CGenebase();
        CGenebase(GeneVAR);
        virtual ~CGenebase();
        virtual double decode();
        virtual void init(GeneVAR);
        virtual void updatecode(double);
        virtual int bitNum();
        virtual Bitset& bitGene();
        virtual double realGene();
        double  value();
    public:
        double  m_value;
        GeneVAR *m_GeneVAR;

};


}
#endif // CODE_H
