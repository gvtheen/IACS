#ifndef CExeVASP_H
#define CExeVASP_H
#include<string>
#include<vector>
#include "CExeFitnessInterface.h"

namespace CALCZJUT{

class CExeVASP:public CExeFitnessInterface
{
	public:
		 CExeVASP(CParameter*);
		~CExeVASP();

		 //virtual function
		 void init();
		 double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
		 void ConvOrigToRawScore(std::vector<double>&);

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
