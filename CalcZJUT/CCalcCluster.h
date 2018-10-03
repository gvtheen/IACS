#ifndef CCALCCLUSTER_H
#define CCALCCLUSTER_H

/*
This class is used to GA evolution of pure cluster rather than molecular adsorbent.
*/
#include "../CataZJUT/CConfigurationBase.h"
#include "CCalcModeStruct.h"
#include "../GaZJUT/GaDeclaration.h"
#include "../Util/Point-Vector.h"
#include "CClusterGeneVar.h"

using util::Point3i;
using GAZJUT::GENEVAR;

namespace CATAZJUT{
  class CPeriodicFramework;
  class CConfigurationBase;
  class CAtom;
  class CElement;
}
namespace CALCZJUT{

class CCalcFitnessInterface;
class CCalcModeStruct;
class CParameter;

class CCalcCluster:public CCalcModeStruct
{
    public:
        CCalcCluster(CParameter* );
        virtual ~CCalcCluster();

        //virtual function from CCalcModeStruct
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        void getGeneValuefromStruct(std::vector<double>&);
        void GeneVarRange(std::vector<GENEVAR>&);

        void init();
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
        std::vector<Point3i>                         currentConnection;
        std::vector<std::vector<CClusterGeneVar*>>   GeneStorage;

        std::vector<std::pair<std::string,size_t>>   chemicalFormula;
        //CCalcFitnessInterface  *m_pCalcFitness;

};


}
#endif // CCalcCluster_H
