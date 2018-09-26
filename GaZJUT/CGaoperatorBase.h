#ifndef CGaoperatorBase_H
#define CGaoperatorBase_H
#include "CGpopulation.h"

namespace GAZJUT{

class CGaoperatorBase
{
	public:
		CGaoperatorBase();
		virtual ~CGaoperatorBase();
		virtual void run(CGpopulation*)=0;
	protected:
};

}
#endif
