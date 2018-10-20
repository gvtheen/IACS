#ifndef CExeDMol_H
#define CExeDMol_H
#include "CExeFitnessInterface.h"
using namespace GAZJUT;

namespace CALCZJUT{

class CExeDMol:public CExeFitnessInterface
{
	public:
		CExeDMol(CParameter*);
		~CExeDMol();

		 CExeFitnessInterface* clone();

         void init();
		 double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
		 void ConvOrigToRawScore(std::vector<double>&);

         void   CheckInputFile();
         bool   IsNormalComplete();
         double readFinalEnergy();
         void   getRelaxedGeometryCoord();

         std::string& inputFile()const;
         void setInputFile(const std::string&);
	protected:
	     std::vector<std::string*>   m_pInputFile;
	     std::string *dmol_inputfile;
};



}
#endif
