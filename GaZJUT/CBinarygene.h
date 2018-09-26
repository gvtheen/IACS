#ifndef CBinarygene_H
#define CBinarygene_H
#include "CGenebase.h"
#include "GaDeclaration.h"

namespace GAZJUT{

class CBinarygene:public CGenebase
{
    public:
        CBinarygene();
        CBinarygene(GENEVAR);
        virtual ~CBinarygene();
        virtual double decode();
        virtual void init(GENEVAR);
        virtual void updatecode(double);
        virtual std::vector<unsigned int>* bitGene();
        virtual int bitNum();
    protected:
        std::vector<unsigned int> *m_bitdata;
        int m_bitNum;
    private:
};

}
#endif // CBinarygene_H
