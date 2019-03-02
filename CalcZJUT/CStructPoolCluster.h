#ifndef CModelClusterSTRUCTPOOL_H
#define CModelClusterSTRUCTPOOL_H
#include <vector>
#include "CStructPoolBase.h"
#include "../GaZJUT/GaDeclaration.h"

using GAZJUT::GeneVAR;

namespace CATAZJUT{
  class CConfigurationBase;
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
        void RandomBuildFromChemicalFormula(CATAZJUT::CConfigurationBase* strut,
                                            std::vector<std::pair<std::string,size_t>>&);
        void eliminateCloseContacts(CATAZJUT::CConfigurationBase* strut,double distanceCutOff=1.0);
        void eliminateFragment(CATAZJUT::CConfigurationBase*);
        size_t ClusterType(std::vector<CATAZJUT::CElement*>&);

        //predicting methods: sphere,
        void metalClusterPredict(CATAZJUT::CConfigurationBase*,
                                 std::vector<std::pair<std::string,size_t>>&);
        void nonMetalClusterPredict(CATAZJUT::CConfigurationBase*,
                                    std::vector<std::pair<std::string,size_t>>&);
        void mixedClusterPredict(CATAZJUT::CConfigurationBase*,
                                 std::vector<std::pair<std::string,size_t>>&);
        void outPutStructure(CATAZJUT::CConfigurationBase*);

    private:

        std::vector<std::pair<std::string,size_t>>   chemicalFormula;
};


}
#endif // CModelClusterSTRUCTPOOL_H
