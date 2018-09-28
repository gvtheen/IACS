#ifndef CCALCDMOL_H
#define CCALCDMOL_H
#include "CCalcFitnessInterface.h"
using namespace GAZJUT;

namespace CALCZJUT{

class CCalcDMol:public CCalcFitnessInterface
{
	public:
		CCalcDMol(CParameter*);
		~CCalcDMol();

         void init();
		 double CalcuRawFit(std::vector<double>* RealValueOfGenome);
         void   ConvOrigToRawScore(std::vector<double>* OrigRawScore);

         void   CheckInputFile();
         bool   IsNormalComplete();
         double readFinalEnergy();
         void getRelaxedGeometryCoord();
	protected:
	     std::vector<std::string*>   m_pInputFile;
	     std::string *dmol_inputfile;
};



}
#endif
