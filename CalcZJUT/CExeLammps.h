#ifndef CExeLammps_H
#define CExeLammps_H

namespace CALCZJUT{

class CExeLammps:public CExeFitnessInterface
{
    public:
        CExeLammps();
        virtual ~CExeLammps();

        double CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist);
		 void ConvOrigToRawScore(std::vector<double>&);

    protected:

    private:

};


}
#endif // CExeLammps_H
