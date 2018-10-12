#ifndef CCALCCLUSTER_H
#define CCALCCLUSTER_H

/*
This class is used to GA evolution of pure cluster rather than molecular adsorbent.
*/
#include "../CataZJUT/CConfigurationBase.h"
#include "CCalcModeStruct.h"
#include "../GaZJUT/GaDeclaration.h"
#include "../Util/Point-Vector.h"
#include "CClusterGeneVAR.h"

using util::Point3i;
using GAZJUT::GeneVAR;

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

        CCalcModeStruct* clone();

        //virtual function from CCalcModeStruct
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        void getGeneValuefromStruct(std::vector<double>&);
        void GeneVARRange(std::vector<GeneVAR>&);

        std::vector<std::pair<std::string,size_t>>& chemicalFormula();
        void setChemicalFormula(const std::vector<std::pair<std::string,size_t>>&);

    protected:

        void eliminateCloseContacts(CATAZJUT::CPeriodicFramework* strut,double distanceCutOff=1.0);
        void eliminateFragment(CATAZJUT::CPeriodicFramework*);

    private:


        //CCalcFitnessInterface  *m_pCalcFitness;

};


}
#endif // CCalcCluster_H
