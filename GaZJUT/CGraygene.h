#ifndef CGraygene_H
#define CGraygene_H
#include "CGenebase.h"
#include "GaDeclaration.h"
#include "CRandomgenerator.h"
namespace GAZJUT{

class CGraygene:public CGenebase
{
    public:
        CGraygene();
        CGraygene(GENEVAR);
        virtual ~CGraygene();
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
#endif // CGraygene_H
