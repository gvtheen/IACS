#ifndef CGAENGINE_H
#define CGAENGINE_H
#include<iostream>
#include<vector>
#include "GaDeclaration.h"
#include "GaUtilityFunction.h"

namespace CALCZJUT{
    class CCalcFitnessInterface;
    class CParameter;
}

namespace GAZJUT{

class CGpopulation;
class CGaoperatorBase;
class CGaparameter;

class CGAEngine
{
	public:
		CGAEngine(CALCZJUT::CParameter*);
		~CGAEngine();

		void SetFitCalculator(const CalcFitnessInterface*);
        CALCZJUT::CalcFitnessInterface* FitCalculator()const;



        void SetGeneVar(const std::vector<GENEVAR>*);
        std::vector<GENEVAR>* GeneVar()const;

		void init();
		void evolve();

	protected:
	    std::vector<CGaoperatorBase*>  m_GeneticOperator;
        CGpopulation                  *m_pCurrentPopulation;
CALCZJUT::CCalcFitnessInterface       *m_pFitnessCalculator;
	    CGaparameter                  *m_pGaparameter;
	    CALCZJUT::CParameter          *m_pParameter;
	    std::vector<GENEVAR>           m_GeneVar;

};


}
#endif
