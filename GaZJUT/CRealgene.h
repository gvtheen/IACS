#ifndef READCODE_H
#define READCODE_H
#include "CGenebase.h"
namespace GAZJUT{

class CRealgene:public CGenebase
{
    public:
        CRealgene();
        CRealgene(GeneVAR);
        virtual ~CRealgene();
        virtual double decode();
        virtual void init(GeneVAR);
        virtual void updatecode(double);
        virtual double& realGene();
        virtual int bitNum();
    protected:

    private:
};


}
#endif // READCODE_H
