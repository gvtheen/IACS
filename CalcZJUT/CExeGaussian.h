#ifndef CExeGaussian_H
#define CExeGaussian_H

#include<vector>
#include "CExeFitnessInterface.h"
using namespace GAZJUT;


namespace CALCZJUT{

class CExeGaussian:public CExeFitnessInterface
{
	public:
		CExeGaussian(CParameter* mpara);
		~CExeGaussian();

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

         std::string& inputFile()const;
         void setInputFile(const std::string&);
	protected:
	     std::string *m_pInputFile;
	     std::string *m_pParaFileAbsPath;

};



}
#endif
