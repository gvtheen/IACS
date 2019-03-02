#ifndef CModelCluster_H
#define CModelCluster_H

/*
This class is used to GA evolution of pure cluster rather than molecular adsorbent.
*/
#include "../CataZJUT/CConfigurationBase.h"
#include "CModelBase.h"
#include "../GaZJUT/GaDeclaration.h"
#include "../Util/Point-Vector.h"
#include "CClusterGeneVAR.h"

using util::Point3i;
using GAZJUT::GeneVAR;

namespace CATAZJUT{
  class CConfigurationBase;
  class CConfigurationBase;
  class CAtom;
  class CElement;
}

namespace CALCZJUT{

class CExeFitnessInterface;
class CModelBase;
class CParameter;

class CModelCluster:public CModelBase
{
    public:
        CModelCluster(CParameter*,
                      size_t index);
        virtual ~CModelCluster();

        CModelBase* clone();

        //virtual function from CModelBase
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        void getGeneValuefromStruct(std::vector<double>&);
        void GeneVARRange(std::vector<GeneVAR>&);
        std::vector<std::pair<std::string,size_t>>& chemicalFormula();
        void setChemicalFormula(const std::vector<std::pair<std::string,size_t>>&);
        void outputStructureToFile();
    protected:

        void eliminateCloseContacts(CATAZJUT::CConfigurationBase* strut,double distanceCutOff=1.0);
        void eliminateFragment(CATAZJUT::CConfigurationBase*);

    private:


        //NO new members

};


}
#endif // CModelCluster_H
