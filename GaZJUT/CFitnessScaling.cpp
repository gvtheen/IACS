//#include<math.h>
#include <iostream>
#include <math.h>
#include "CFitnessScaling.h"
#include "GaDeclaration.h"

namespace GAZJUT{

CFitnessScaling::CFitnessScaling()
{
}

CFitnessScaling::~CFitnessScaling()
{
}
void CFitnessScaling::run(CGpopulation* currentPopulation)
{
     E_SCALING_TYPE myStype = currentPopulation->m_pObjGaparameter->ScalingMode();
     switch(myStype )
     {
        case LINEAR:
            this->linearScaling(currentPopulation);
            break;
        case SIGMA:
            this->powerLawScaling(currentPopulation);
            break;
        case POWER:
            this->sigmaTruncScaling(currentPopulation);
            break;
        default:
            this->linearScaling(currentPopulation);
            break;
     }
}
void CFitnessScaling::linearScaling(CGpopulation* currentPopulation)
{
    double raw_min, raw_max, raw_avg;
    currentPopulation->raw_statistic();

    raw_min = (*currentPopulation)["minRaw"];
    raw_max = (*currentPopulation)["maxRaw"];
    raw_avg = (*currentPopulation)["avgRaw"];
    int a=0,b=0,delta=0;
    int c = currentPopulation->m_pObjGaparameter->ScaleLinearMultiplier;
    if(raw_avg==raw_max)
    {
       a = 1.0;
       b = 0.0;
    }else if( raw_min > ( c*raw_avg - raw_max/c - 1.0)){
       delta=raw_max-raw_avg;
       a = (c - 1.0) * raw_avg / delta;
       b = raw_avg*(raw_max - c*raw_avg)/delta;
    }else{
       delta = raw_avg - raw_min;
       a = raw_avg / delta;
       b = -1*raw_min * raw_avg / delta;
    }
    double temp_raw;
    int popnum = currentPopulation->popNum();
    //std::vector <CGenome*> *p_genome=currentPopulation->m_pGpopulation;
    for(int i=0;i<popnum;i++)
    {
       temp_raw = ((*currentPopulation)[i])->rawscore();
       temp_raw = temp_raw*a + b;
       if(temp_raw<0)
          temp_raw=0;
       ((*currentPopulation)[i])->setFitness(temp_raw);
    }

}
void CFitnessScaling::sigmaTruncScaling(CGpopulation* currentPopulation)
{
    currentPopulation->raw_statistic();

    double c = currentPopulation->m_pObjGaparameter->ScaleSigmaTruncMultiplier;
    double raw_avg = (*currentPopulation)["avgRaw"];
    double raw_dev = (*currentPopulation)["rawDev"];

    int popnum = currentPopulation->popNum();
    //std::vector <CGenome*> *p_genome=currentPopulation->m_pGpopulation;
    double temp_raw;
    for(int i=0;i<popnum;i++)
    {
       temp_raw = ((*currentPopulation)[i])->rawscore() - raw_avg;
       temp_raw = temp_raw + raw_dev*c;
       ((*currentPopulation)[i])->setFitness(temp_raw);
    }
}
void CFitnessScaling::powerLawScaling(CGpopulation* currentPopulation)
{
    currentPopulation->raw_statistic();

    double k = currentPopulation->m_pObjGaparameter->ScalePowerLawFactor;
    int popnum = currentPopulation->popNum();

    //std::vector <CGenome*> *p_genome=currentPopulation->m_pGpopulation;
    double temp_raw;
    for(int i=0;i<popnum;i++)
    {
       temp_raw = ((*currentPopulation)[i])->rawscore();
       temp_raw = std::pow(temp_raw,k);
       ((*currentPopulation)[i])->setFitness(temp_raw);
    }
}


}
