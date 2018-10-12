#ifndef CCALCCLUSTERSTRUCTPOOL_H
#define CCALCCLUSTERSTRUCTPOOL_H
#include <vector>
#include "CCalcStructBasePool.h"
#include "../GaZJUT/CGaDeclaration.h"

using GAZJUT::GeneVAR;

namespace CATAZJUT{
  class CPeriodicFramework;
  class CElement;
}

namespace CALCZJUT{

class CCalcClusterStructPool:public CCalcStructureBasePool
{
    public:
        CCalcClusterStructPool(CParameter* othr);
        virtual ~CCalcClusterStructPool();

        void init();
        void GeneVARRange(std::vector<GeneVAR>&);

        void Initialization(const std::string& chemicalformulaStr);          // initialize from chemical formula
        void Initialization(const char* chemicalformulaStr);                 // initialize from chemical formula
        void Initialization(const std::vector<std::string*>& inputfiles);     // initialize from exit structure

    protected:
        void RandomBuildFromChemicalFormula(CATAZJUT::CPeriodicFramework* strut);
        void eliminateCloseContacts(CATAZJUT::CPeriodicFramework* strut,double distanceCutOff=1.0);
        void eliminateFragment(CATAZJUT::CPeriodicFramework*);
        size_t ClusterType(std::vector<CATAZJUT::CElement*>&);

        //predicting methods: sphere,
        void metalClusterPredict(CATAZJUT::CPeriodicFramework*);
        void nonMetalClusterPredict(CATAZJUT::CPeriodicFramework*);
        void mixedClusterPredict(CATAZJUT::CPeriodicFramework*);

    private:

        std::vector<std::pair<std::string,size_t>>   chemicalFormula;
};


}
#endif // CCALCCLUSTERSTRUCTPOOL_H
