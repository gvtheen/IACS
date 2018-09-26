#ifndef CEVALUATOR_H
#define CEVALUATOR_H
#include "CGpopulation.h"
#include "CGaoperatorBase.h"
#include "CCalcFitnessInterface.h"
#include "CCalcModeStruct.h"
/*
The main function of class is that FitnessValuator conducts calculation to
 get originalValue and RawScore.

*/

namespace GAZJUT{

class CEvaluator:public CGaoperatorBase
{
	public:
		CEvaluator();
		CEvaluator(CCalcFitnessInterface*);
		~CEvaluator();

		void clone();
		void init();
		virtual void run( CGpopulation* );

        void setCalcFitnessInterface(CALCZJUT::CCalcFitnessInterface* CalcModeStruct);
        CALCZJUT::CCalcFitnessInterface* CalcFitnessInterface();
	protected:
	    CALCZJUT::CCalcFitnessInterface      *m_pEvaluator;
};

}
#endif
