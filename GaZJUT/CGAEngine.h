#ifndef CGAENGINE_H
#define CGAENGINE_H
#include<iostream>
#include<vector>
#include "GaDeclaration.h"
#include "GaUtilityFunction.h"

namespace CALCZJUT{
    class CExeFitnessInterface;
    class CParameter;
    class CStructPoolBase;
}

namespace GAZJUT{

class CGpopulation;
class CGaOperatorBase;
class CGaparameter;

class CGAEngine
{
	public:
		CGAEngine(CALCZJUT::CParameter*);
		~CGAEngine();

//		void SetFitCalculator(CALCZJUT::CExeFitnessInterface*);
//        CALCZJUT::CExeFitnessInterface* FitnessCalculator()const;



        void SetGeneVAR(const std::vector<GeneVAR>&);
        std::vector<GeneVAR>& GeneVAR()const;

		void init();
		void evolve();

	protected:
	    std::vector<CGaOperatorBase*>                 m_GeneticOperator;
        CGpopulation                                 *m_pCurrentPopulation;

        std::vector<CALCZJUT::CExeFitnessInterface*>  m_FitnessCalculator;
        CALCZJUT::CStructPoolBase                    *m_pStructurePool;
	    CALCZJUT::CParameter                         *m_pParameter;
	    std::vector<GAZJUT::GeneVAR>                  m_GeneVAR;

};


}
#endif
