#ifndef CROSS_H
#define CROSS_H
#include "CGpopulation.h"
#include "CGaOperatorBase.h"
namespace GAZJUT{

class CCross:public CGaOperatorBase
{
	public:
		CCross();
		CCross(double);
		~CCross();

		virtual void run(CGpopulation*);
		void X_Crossover(CGpopulation*,int,int);
	protected:
		double m_Crossprob;
};


}
#endif
