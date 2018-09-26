#ifndef CCALCVASP_H
#define CCALCVASP_H
#include<string>
#include<vector>
#include "CCalcFitnessInterface.h"

namespace CALCZJUT{

class CCalcVASP:public CCalcFitnessInterface
{
	public:
		 CCalcVASP(CParameter*);
		~CCalcVASP();

		 //virtual function
		 void init();
		 double CalcuRawFit(std::vector<double>* RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
		 void ConvOrigToRawScore(std::vector<double>*);

		 //normal function
		 void   CheckInputFile();
         bool   IsNormalComplete();
         double readFinalEnergy();
         void getRelaxedGeometryCoord();

	protected:
		 std::vector<std::string*>   m_pInputFile;
};



}//NAMESPACE
#endif
