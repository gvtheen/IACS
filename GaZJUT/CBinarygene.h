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
        CBinarygene(GeneVAR);
        virtual ~CBinarygene();
        virtual double decode();
        virtual void init(GeneVAR);
        virtual void updatecode(double);
        virtual Bitset& bitGene();
        virtual int bitNum();
    protected:
        Bitset m_bitdata;
        int m_bitNum;
    private:
};

}
#endif // CBinarygene_H
