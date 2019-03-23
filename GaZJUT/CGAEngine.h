#ifndef CGAENGINE_H
#define CGAENGINE_H
#include<iostream>
#include<vector>
#include "GaDeclaration.h"
#include "GaUtilityFunction.h"
#include "../CAIEngineBase.h"
#include "../IACS.h"

using IACSZJUT::VarRangeStruct;
class CAIEngineBase;

namespace CALCZJUT{
    class CExeFitnessInterface;
    class CParameter;
    class CStructPoolBase;
}

namespace GAZJUT{

class CGpopulation;
class CGaOperatorBase;
class CGaparameter;

class CGAEngine:public CAIEngineBase
{
	public:
		CGAEngine(CALCZJUT::CParameter*);
		~CGAEngine();

//        void SetVarRange(const std::vector<IACSZJUT::VarRangeStruct>&);
//        std::vector<IACSZJUT::VarRangeStruct>& VarRange()const;
		void init();
		void evolve();

	protected:
	    std::vector<CGaOperatorBase*>                 m_GeneticOperator;
        CGpopulation                                 *m_pCurrentPopulation;
        CGaparameter                                 *m_pGaparameter;


};


}
#endif
