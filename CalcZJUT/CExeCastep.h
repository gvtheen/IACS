#ifndef CEXECASTEP_H
#define CEXECASTEP_H

namespace CALCZJUT{

class CExeCastep:public CExeFitnessInterface
{
    public:
        CExeCastep(CParameter*);
        virtual ~CExeCastep();

         CExeFitnessInterface* clone();

         void init();
		 double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
		 void ConvOrigToRawScore(std::vector<double>&);
         char* ExeName();

         void   CheckInputFile();
         bool   IsNormalComplete();
         double readFinalEnergy();
         void   getRelaxedGeometryCoord();

         std::string& inputFile()const;
         void setInputFile(const std::string&);
    protected:

    private:
        std::vector<std::string>   m_pInputFile;
};


}
#endif // CEXECASTEP_H
