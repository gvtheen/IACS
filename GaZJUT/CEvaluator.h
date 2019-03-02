#ifndef CEVALUATOR_H
#define CEVALUATOR_H
#include <vector>
#include <map>
#include <iostream>

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
	    typedef enum{OPTIMAL_ENERGY = 0x9231A,
                        OPTIMAL_GAP = 0x9232A
	                 }OPT_TYPE;
		CEvaluator();
		CEvaluator(std::vector<CALCZJUT::CExeFitnessInterface*>,
                   CALCZJUT::CStructPoolBase*);
		~CEvaluator();

		void clone();
		void init();
		void run( CGpopulation* );


//        void setCalcFitnessInterface(CALCZJUT::CExeFitnessInterface* CalcModeStruct);
//        CALCZJUT::CExeFitnessInterface* CalcFitnessInterface();

	protected:
        void getTargetPopWithCondition(std::vector<size_t>& res,std::map <size_t, double>& mapIndexValue,
                                                       size_t num=1,double conditionValue = 0);
        void standardOutput(std::map <size_t, double>&);
        void standardOutput(CGpopulation*);
    private:
	    std::vector<CALCZJUT::CExeFitnessInterface*>      m_pEvaluatorPool;
	                  CALCZJUT::CExeFitnessInterface      *m_currentEvaluator;
	                        CALCZJUT::CStructPoolBase     *m_pStructurePool;
	                //CALCZJUT::CParameter     *m_pParameter;
	                                          Bitset       pop_run_state;
};

}
#endif
