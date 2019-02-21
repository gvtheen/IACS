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
		 CExeFitnessInterface* clone();
		 void init();
		 double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
		 void ConvOrigToRawScore(std::vector<double>&);
         std::string ExeName();
		 //normal function
		 void   CheckInputFile();
         bool   IsNormalComplete();
         double readFinalEnergy();
         void   getRelaxedGeometryCoord();

	protected:
		 std::vector<std::string*>   m_pInputFile;
         std::string *m_pParaFileAbsPath;
};



}//NAMESPACE
#endif
