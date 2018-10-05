#ifndef CMUTATOR_H
#define CMUTATOR_H
#include "CGpopulation.h"
#include "CGaOperatorBase.h"

namespace GAZJUT{

class CMutator:public CGaOperatorBase
{
	public:
		CMutator();
		~CMutator();
		virtual void run(CGpopulation*);
	protected:

};

}
#endif
