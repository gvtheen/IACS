#ifndef CCALCGAUSSIAN_H
#define CCALCGAUSSIAN_H
#include "CCalcFitnessInterface.h"
using namespace GAZJUT;


namespace CALCZJUT{

class CCalcGaussian:public CCalcFitnessInterface
{
	public:
		CCalcGaussian(CParameter* mpara);
		~CCalcGaussian();

		 void init();
		 double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
		 void ConvOrigToRawScore(std::vector<double>&);

		 //normal function
		 void   CheckInputFile();
         bool   IsNormalComplete();
         double readFinalEnergy();
         void   getRelaxedGeometryCoord();
	protected:
	    std::string *m_pInputFile;

};



}
#endif
