#ifndef CELIST_H
#define CELIST_H
#include "CGenome.h"
#include "CGpopulation.h"
#include "CGaoperatorBase.h"

namespace GAZJUT{

class CElist:public CGaoperatorBase
{
	public:
		CElist();
		~CElist();

        virtual void run(CGpopulation*);
	protected:
	    CGpopulation  *m_OldPopulation;
        CGenome       *m_pBestGenome;
	    double         m_MinOrignalValue;
        //double        m_BestOrigValue;    // For DFT, is is energy of relaxed configuration.
};


}
#endif
