#ifndef CBinarygene_H
#define CBinarygene_H
#include "CGenebase.h"
#include "GaDeclaration.h"
#include "../Util/Bitset.h"

using util::Bitset;

namespace GAZJUT{

class CBinarygene:public CGenebase
{
    public:
        CBinarygene();
        CBinarygene(GeneVAR*);
        virtual ~CBinarygene();
        double decode();
        void init(GeneVAR*);
        void updatecode(double);
        Bitset& bitGene();
        size_t bitNum();
    protected:
        Bitset m_bitdata;
        size_t m_bitNum;
    private:
};

}
#endif // CBinarygene_H
