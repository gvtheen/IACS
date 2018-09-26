#ifndef CMUTATOR_H
#define CMUTATOR_H
#include "CGpopulation.h"

namespace GAZJUT{

class CMutator:public CGaoperatorBase
{
	public:
		CMutator();
		~CMutator();
		virtual void run(CGpopulation*);
	protected:

};

}
#endif
