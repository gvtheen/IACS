#ifndef CEVALUATOR_H
#define CEVALUATOR_H
#include "CGpopulation.h"
#include "CGaOperatorBase.h"
#include "../Util/Bitset.h"
/*
The main function of class is that FitnessValuator conducts calculation to
 get originalValue and RawScore.

*/

using util::Bitset;

namespace CALCZJUT{
  class CCalcFitnessInterface;
  class CCalcModeStruct;
}
namespace GAZJUT{

class CEvaluator:public CGaOperatorBase
{
	public:
		CEvaluator();
		CEvaluator(CALCZJUT::CCalcFitnessInterface*);
		~CEvaluator();

		void clone();
		void init();
		void run( CGpopulation* );

        void setCalcFitnessInterface(CALCZJUT::CCalcFitnessInterface* CalcModeStruct);
        CALCZJUT::CCalcFitnessInterface* CalcFitnessInterface();
	protected:
	    CALCZJUT::CCalcFitnessInterface      *m_pEvaluator;
	                             Bitset       pop_run_state;
};

}
#endif
