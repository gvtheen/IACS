#ifndef CCALCLAMMPS_H
#define CCALCLAMMPS_H

namespace CALCZJUT{

class CCalcLammps:public CCalcFitnessInterface
{
    public:
        CCalcLammps();
        virtual ~CCalcLammps();

        double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
		 void ConvOrigToRawScore(std::vector<double>&);

    protected:

    private:

};


}
#endif // CCALCLAMMPS_H
