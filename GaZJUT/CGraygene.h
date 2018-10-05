#ifndef CGraygene_H
#define CGraygene_H
#include "CGenebase.h"
#include "GaDeclaration.h"
#include "../Util/CRandomgenerator.h"

#include "../Util/Bitset.h"

using util::Bitset;

namespace GAZJUT{

class CGraygene:public CGenebase
{
    public:
        CGraygene();
        CGraygene(GeneVAR*);
        virtual ~CGraygene();
        virtual double decode();
        virtual void init(GeneVAR*);
        virtual void updatecode(double);
        virtual Bitset& bitGene();
        virtual size_t bitNum();
    protected:
        Bitset m_bitdata;
        size_t m_bitNum;
    private:

};


}
#endif // CGraygene_H
