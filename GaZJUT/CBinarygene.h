#ifndef CBinarygene_H
#define CBinarygene_H
#include "CGenebase.h"
#include "GaDeclaration.h"
#include "../Util/Bitset.h"
#include "../IACS.h"

using IACSZJUT::VarRangeStruct;
using util::Bitset;

namespace GAZJUT{

class CBinarygene:public CGenebase
{
    public:
        CBinarygene();
        CBinarygene(VarRangeStruct*);
        virtual ~CBinarygene();

        double decode();
        void init(VarRangeStruct*);
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
