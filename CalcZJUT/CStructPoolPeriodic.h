#ifndef CSTRUCTPOOLPERIODIC_H
#define CSTRUCTPOOLPERIODIC_H
#include <vector>
#include "CStructPoolBase.h"
#include "../GaZJUT/GaDeclaration.h"

namespace CALCZJUT{

class CStructPoolPeriodic:public CStructPoolBase
{
    public:
        CStructPoolPeriodic(CParameter* othr);
        virtual ~CStructPoolPeriodic();
        //virtual function
        void init();
        void GeneVARRange(std::vector<GeneVAR>&);
       //normal function
        void Initialization(const std::string& chemicalformulaStr);          // initialize from chemical formula
        void Initialization(const char* chemicalformulaStr);                 // initialize from chemical formula
        void Initialization(const std::vector<std::string*>& inputfiles);     // initialize from exit structure

    protected:


    private:
        void RandomBuildFromChemicalFormula(CATAZJUT::CPeriodicFramework* strut,
                                            std::vector<std::pair<std::string,size_t>>&);
};


}
#endif // CSTRUCTPOOLPERIODIC_H
