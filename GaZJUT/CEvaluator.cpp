#include "CEvaluator.h"
#include "../CalcZJUT/CCalcFitnessInterface.h"
#include "../CalcZJUT/CCalcModeStruct.h"
namespace GAZJUT{

CEvaluator::CEvaluator()
{
}
CEvaluator::CEvaluator(CCalcFitnessInterface* myEvaluator)
{
    this->m_pEvaluator = myEvaluator->clone();
}
CEvaluator::~CEvaluator()
{
}
void CEvaluator::run(CGpopulation* CurrentPopulation)
{
   bool runstate=false;
   int pop_num = CurrentPopulation->popNum();

   std::vector<double> *Temp_OrigScore = new (std::vector<double>)(pop_num);
   // calculate the fitness of each genome in population
   for(int i=0;i<pop_num;i++)
   {
      Temp_OrigScore->at(i) = m_pEvaluator->CalcuRawFit(((*CurrentPopulation)[i])->getDecValue(),runstate);
      ((*CurrentPopulation)[i])->setOrigValue( Temp_OrigScore->at(i) );
      ((*CurrentPopulation)[i])->setFinishState(runstate);
   }
   // call the function of m_pEvaluator for converting the original value to Raw score.
   m_pEvaluator->ConvOrigToRawScore(Temp_OrigScore);

   for(int i=0;i<pop_num;i++)
      ((*CurrentPopulation)[i])->setRawScore(Temp_OrigScore->at(i));
}
void CEvaluator::setCalcFitnessInterface(CALCZJUT::CCalcFitnessInterface* CalcFitness)
{
   this->m_pEvaluator=CalcFitness;
}
CALCZJUT::CCalcFitnessInterface* CEvaluator::CalcFitnessInterface()
{
   return this->m_pEvaluator;
}
}
