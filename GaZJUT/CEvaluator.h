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
  class CCalcStructureBasePool;
  class CParameter;
}
namespace GAZJUT{

class CEvaluator:public CGaOperatorBase
{
	public:
		CEvaluator();
		CEvaluator(CALCZJUT::CCalcFitnessInterface*,
                   CALCZJUT::CCalcStructureBasePool*);
		~CEvaluator();

		void clone();
		void init();
		void run( CGpopulation* );

        void setCalcFitnessInterface(CALCZJUT::CCalcFitnessInterface* CalcModeStruct);
        CALCZJUT::CCalcFitnessInterface* CalcFitnessInterface();

	protected:
	    CALCZJUT::CCalcFitnessInterface      *m_pEvaluator;
	    CALCZJUT::CCalcStructureBasePool     *m_pStructurePool;
	                //CALCZJUT::CParameter     *m_pParameter;
	                             Bitset       pop_run_state;
};

}
#endif
