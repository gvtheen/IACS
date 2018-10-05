#ifndef CSELECTOR_H
#define CSELECTOR_H
#include "CGpopulation.h"
#include "CGenome.h"
#include "CGaOperatorBase.h"

namespace GAZJUT{
// declarization

class CSelector:public CGaOperatorBase
{
	public:
		CSelector();
		CSelector(const CSelector&);
		~CSelector();

		CSelector* clone();
		void init();
		void run(CGpopulation*);
		// three different selectors
		void roulette_Wheel_Select(CGpopulation*);
		void rondom_Select(CGpopulation*);
		void tournament_Select(CGpopulation*); //TOURNAMENT

	protected:


};


}
#endif
