#ifndef CCALCCLUSTER_H
#define CCALCCLUSTER_H

/*
This class is used to GA evolution of pure cluster rather than molecular adsorbent.
*/
#include "CConfigurationBase.h"
#include "CCalcModeStruct.h"
#include "GaDeclaration.h"
#include "Point-Vector.h"

using CATAZJUT::Point3i;
using GAZJUT::GENEVAR;
namespace CATAZJUT{
  class CPeriodicFramework;
  class CConfigurationBase;
  class CAtom;
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
        std::vector<double>*  getGeneValuefromStruct()const;
        std::vector<GENEVAR>* GeneVarRange();

        void Initialization(std::string&);        // initialize from chemical formula
        void Initialization(char*);               // initialize from chemical formula
        void Initialization(CConfigurationBase&); // initialize from exit structure

    protected:

    private:
        std::vector<Point3i> currentConnection;
        CCalcFitnessInterface *m_pCalcFitness;


};


}
#endif // CCalcCluster_H
