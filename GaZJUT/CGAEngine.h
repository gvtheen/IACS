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
        CALCZJUT::CalcFitnessInterface* FitnessCalculator()const;



        void SetGeneVAR(const std::vector<GeneVAR>*);
        std::vector<GeneVAR>* GeneVAR()const;

		void init();
		void evolve();

	protected:
	    std::vector<CGaoperatorBase*>  m_GeneticOperator;
        CGpopulation                  *m_pCurrentPopulation;
CALCZJUT::CCalcFitnessInterface       *m_pFitnessCalculator;
	    CGaparameter                  *m_pGaparameter;
	    CALCZJUT::CParameter          *m_pParameter;
	    std::vector<GeneVAR>           m_GeneVAR;

};


}
#endif
