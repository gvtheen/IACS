#ifndef CEXEDFTB_H
#define CEXEDFTB_H

#include<vector>
#include<iostream>
#include<string>
#include "CExeFitnessInterface.h"

namespace CALCZJUT{

class CExeDFTB:public CExeFitnessInterface
{
    public:
        CExeDFTB(CParameter*);
        virtual ~CExeDFTB();

         CExeFitnessInterface* clone();

         void init();
		 double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
		 void ConvOrigToRawScore(std::vector<double>&);
         std::string ExeName();

         void   CheckInputFile();
         bool   IsNormalComplete();
         double readFinalEnergy();
         void   getRelaxedGeometryCoord();

         std::string& inputFile()const;
         void setInputFile(const std::string&);
    protected:

    private:
};


}
#endif // CEXEDFTB_H
