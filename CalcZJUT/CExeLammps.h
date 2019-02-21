#ifndef CExeLammps_H
#define CExeLammps_H
#include<vector>
#include "CExeFitnessInterface.h"

namespace CALCZJUT{

class CExeLammps:public CExeFitnessInterface
{
    public:
        CExeLammps(CParameter* mpara);
        virtual ~CExeLammps();

        double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
        void ConvOrigToRawScore(std::vector<double>&);
        std::string ExeName();
        CExeFitnessInterface* clone();

    protected:

    private:

};


}
#endif // CExeLammps_H
