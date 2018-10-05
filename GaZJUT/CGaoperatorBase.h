#ifndef CGaOperatorBase_H
#define CGaOperatorBase_H
#include "CGpopulation.h"

namespace GAZJUT{

class CGaOperatorBase
{
	public:
		CGaOperatorBase();
		virtual ~CGaOperatorBase();
		virtual void run(CGpopulation*)=0;
	protected:
};

}
#endif
