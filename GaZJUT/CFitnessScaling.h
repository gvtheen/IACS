#ifndef CFITNESSSCALING_H
#define CFITNESSSCALING_H
#include "CGaOperatorBase.h"
#include "CGpopulation.h"
namespace GAZJUT{

class CFitnessScaling:public CGaOperatorBase
{
	public:
		CFitnessScaling();
		~CFitnessScaling();

		void clone();
		void init();
		void run(CGpopulation*);

    private:
		void linearScaling(CGpopulation*);
		void sigmaTruncScaling(CGpopulation*);
		void powerLawScaling(CGpopulation*);
		void noScaling(CGpopulation*);
	protected:
};


}
#endif
