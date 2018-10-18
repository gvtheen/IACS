#ifndef CExeGaussian_H
#define CExeGaussian_H
#include "CExeFitnessInterface.h"
using namespace GAZJUT;


namespace CALCZJUT{

class CExeGaussian:public CExeFitnessInterface
{
	public:
		CExeGaussian(CParameter* mpara);
		~CExeGaussian();

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
