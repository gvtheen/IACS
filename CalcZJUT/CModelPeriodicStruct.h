#ifndef CModelPeriodicStruct_H
#define CModelPeriodicStruct_H
#include "CModelBase.h"


namespace CALCZJUT{

class CModelPeriodicStruct:public CModelBase
{
    public:
        CModelPeriodicStruct();
        virtual ~CModelPeriodicStruct();

        CModelBase* clone();

        //virtual function derived from CModelBase
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        void getGeneValuefromStruct(std::vector<double>&);
        void GeneVARRange(std::vector<GeneVAR>&);

        std::vector<std::pair<std::string,size_t>>& chemicalFormula();
        void setChemicalFormula(const std::vector<std::pair<std::string,size_t>>&);

    protected:

        void eliminateCloseContacts(CATAZJUT::CPeriodicFramework* strut,double distanceCutOff=1.0);
        void eliminateFragment(CATAZJUT::CPeriodicFramework*);
    private:
};



}
#endif // CModelPeriodicStruct_H
