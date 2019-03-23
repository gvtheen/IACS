#ifndef READCODE_H
#define READCODE_H
#include "CGenebase.h"
#include "../IACS.h"

using IACSZJUT::VarRangeStruct;

namespace GAZJUT{

class CRealgene:public CGenebase
{
    public:
        CRealgene();
        CRealgene(VarRangeStruct*);
        virtual ~CRealgene();
        double decode();
        void init(VarRangeStruct*);
        void updatecode(double);
        double& realGene();
        size_t bitNum();
    protected:

    private:
};


}
#endif // READCODE_H
