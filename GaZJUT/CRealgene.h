#ifndef READCODE_H
#define READCODE_H
#include "CGenebase.h"
namespace GAZJUT{

class CRealgene:public CGenebase
{
    public:
        CRealgene();
        CRealgene(GeneVAR*);
        virtual ~CRealgene();
        double decode();
        void init(GeneVAR*);
        void updatecode(double);
        double& realGene();
        size_t bitNum();
    protected:

    private:
};


}
#endif // READCODE_H
