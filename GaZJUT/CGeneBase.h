#ifndef CODE_H
#define CODE_H
#include<vector>
#include "GaDeclaration.h"

namespace GAZJUT{

class CGenebase
{
    public:
        CGenebase();
        CGenebase(GENEVAR);
        virtual ~CGenebase();
        virtual double decode();
        virtual void init(GENEVAR);
        virtual void updatecode(double);
        virtual int bitNum();
        virtual std::vector<unsigned int>* bitGene();
        virtual double& realGene();
        double  value();
    public:
        double  m_value;
        GENEVAR *m_geneVar;

};


}
#endif // CODE_H
