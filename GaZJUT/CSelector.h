#ifndef CSELECTOR_H
#define CSELECTOR_H
#include "CGpopulation.h"
#include "CGenome.h"
#include "CGene.h"
#include "CRandomgenerator.h"
namespace GAZJUT{

class CSelector:public CGaoperatorBase
{
	public:
		CSelector();
		CSelector(const CSelector&);
		~CSelector();

		void clone();
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
