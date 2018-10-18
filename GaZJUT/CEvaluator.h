#ifndef CEVALUATOR_H
#define CEVALUATOR_H
#include "CGpopulation.h"
#include "CGaOperatorBase.h"
#include "../Util/Bitset.h"
#include "GaDeclaration.h"
/*
The main function of class is that FitnessValuator conducts calculation to
 get originalValue and RawScore.

*/

using util::Bitset;

namespace CALCZJUT{
  class CExeFitnessInterface;
  class CModelBase;
  class CStructPoolBase;
  class CParameter;
}
namespace GAZJUT{

class CEvaluator:public CGaOperatorBase
{
	public:
		CEvaluator();
		CEvaluator(CALCZJUT::CExeFitnessInterface*,
                   CALCZJUT::CStructPoolBase*);
		~CEvaluator();

		void clone();
		void init();
		void run( CGpopulation* );

        void setCalcFitnessInterface(CALCZJUT::CExeFitnessInterface* CalcModeStruct);
        CALCZJUT::CExeFitnessInterface* CalcFitnessInterface();

	protected:
	    CALCZJUT::CExeFitnessInterface      *m_pEvaluator;
	          CALCZJUT::CStructPoolBase     *m_pStructurePool;
	                //CALCZJUT::CParameter     *m_pParameter;
	                             Bitset       pop_run_state;
};

}
#endif
