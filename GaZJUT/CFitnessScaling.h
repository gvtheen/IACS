#ifndef CFITNESSSCALING_H
#define CFITNESSSCALING_H
#include "CGaoperatorBase.h"
#include "CGpopulation.h"
namespace GAZJUT{

class CFitnessScaling:public CGaoperatorBase
{
	public:
		CFitnessScaling();
		~CFitnessScaling();

		void clone();
		void init();
		void run(CGpopulation*);

		// scaling rawscore to fitness
		void linearScaling(CGpopulation*);
		void sigmaTruncScaling(CGpopulation*);
		void powerLawScaling(CGpopulation*);
	protected:
};


}
#endif
