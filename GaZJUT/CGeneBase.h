#ifndef CODE_H
#define CODE_H
#include<vector>
#include "GaDeclaration.h"
#include "../Util/Bitset.h"
#include "../IACS.h"


using IACSZJUT::VarRangeStruct;
using util::Bitset;

namespace GAZJUT{

class CGenebase
{
    public:
        CGenebase();
        CGenebase(VarRangeStruct*);
        virtual ~CGenebase();
        virtual double decode();
        virtual void init(VarRangeStruct*);
        virtual void updatecode(double);
        virtual size_t bitNum();
        virtual Bitset& bitGene();
        virtual double& realGene();
        double  value();
    public:
        double      m_value;
        VarRangeStruct    *m_VarRangeStruct;
        E_CODE_TYPE codeType;

};


}
#endif // CODE_H
