#ifndef CModelPeriodicStruct_H
#define CModelPeriodicStruct_H
#include "CModelBase.h"

namespace CATAZJUT{
   class CConfigurationBase;
}

namespace CALCZJUT{

class CModelPeriodicStruct:public CModelBase
{
    public:
        CModelPeriodicStruct(CParameter*,size_t);
        virtual ~CModelPeriodicStruct();

        CModelBase* clone();

        //virtual function derived from CModelBase
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        void getGeneValuefromStruct(std::vector<double>&);
        void VarRangeStructRange(std::vector<VarRangeStruct>&);

        std::vector<std::pair<std::string,size_t>>& chemicalFormula();
        void setChemicalFormula(const std::vector<std::pair<std::string,size_t>>&);

    protected:

        void eliminateCloseContacts(CATAZJUT::CConfigurationBase* strut,double distanceCutOff=1.0);
        void eliminateFragment(CATAZJUT::CConfigurationBase*);
    private:
};



}
#endif // CModelPeriodicStruct_H
