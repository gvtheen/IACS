#ifndef CModelClusterSTRUCTPOOL_H
#define CModelClusterSTRUCTPOOL_H
#include <vector>
#include "CStructPoolBase.h"
#include "../GaZJUT/GaDeclaration.h"

using GAZJUT::GeneVAR;

namespace CATAZJUT{
  class CPeriodicFramework;
  class CElement;
}

namespace CALCZJUT{

class CStructPoolCluster:public CStructPoolBase
{
    public:
        CStructPoolCluster(CParameter* othr);
        virtual ~CStructPoolCluster();

        //virtual function
        void init();
        void GeneVARRange(std::vector<GeneVAR>&);

        void Initialization(const std::string& chemicalformulaStr);          // initialize from chemical formula
        void Initialization(const char* chemicalformulaStr);                 // initialize from chemical formula
        void Initialization(const std::vector<std::string*>& inputfiles);     // initialize from exit structure

    protected:
        void RandomBuildFromChemicalFormula(CATAZJUT::CPeriodicFramework* strut,
                                            std::vector<std::pair<std::string,size_t>>&);
        void eliminateCloseContacts(CATAZJUT::CPeriodicFramework* strut,double distanceCutOff=1.0);
        void eliminateFragment(CATAZJUT::CPeriodicFramework*);
        size_t ClusterType(std::vector<CATAZJUT::CElement*>&);

        //predicting methods: sphere,
        void metalClusterPredict(CATAZJUT::CPeriodicFramework*,
                                 std::vector<std::pair<std::string,size_t>>&);
        void nonMetalClusterPredict(CATAZJUT::CPeriodicFramework*,
                                    std::vector<std::pair<std::string,size_t>>&);
        void mixedClusterPredict(CATAZJUT::CPeriodicFramework*,
                                 std::vector<std::pair<std::string,size_t>>&);

    private:

        std::vector<std::pair<std::string,size_t>>   chemicalFormula;
};


}
#endif // CModelClusterSTRUCTPOOL_H
