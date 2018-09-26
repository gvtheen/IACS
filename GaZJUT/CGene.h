#ifndef CGENE_H
#define CGENE_H
#include<iostream>

#include "GaDeclaration.h"
#include "GaUtilityFunction.h"


namespace GAZJUT{

class CGene
{
	public:
		//constructor
		CGene();
		CGene(GENEVAR);
		CGene(const CGene&);
		~CGene();

		//operator function
		int    encode(double value);
		double decode();
		CGene* clone();
		void   randomInitial();
		void   init();
		void   init(GENEVAR);
		void   updatebit(std::vector <unsigned int>*);

		//input output function
		int    binaryBitNum();
		double  value();
	protected:
		std::vector <unsigned int> *m_gene;
		int                    m_bitNum;


};

}
#endif
