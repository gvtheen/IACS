#ifndef CCALCLAMMPS_H
#define CCALCLAMMPS_H

namespace CALCZJUT{

class CCalcLammps:public CCalcFitnessInterface
{
    public:
        CCalcLammps();
        virtual ~CCalcLammps();

    protected:

    private:
        GAZJUT::CEvaluator   *m_pEvaluator;
};


}
#endif // CCALCLAMMPS_H
