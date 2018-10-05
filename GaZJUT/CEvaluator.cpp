#include "CEvaluator.h"
#include "../CalcZJUT/CCalcFitnessInterface.h"
#include "../CalcZJUT/CCalcModeStruct.h"
namespace GAZJUT{

CEvaluator::CEvaluator()
{
}
CEvaluator::CEvaluator(CALCZJUT::CCalcFitnessInterface* myEvaluator)
{
    this->m_pEvaluator = myEvaluator->clone();
}
CEvaluator::~CEvaluator()
{
}
void CEvaluator::run(CGpopulation* CurrentPopulation)
{
   bool runstate=false;
   double tempOrigValue;
   size_t pop_num = CurrentPopulation->popNum();

   std::vector<double> OrigScore;
   std::vector<double> DecValueOfGenome;
   // calculate the fitness of each genome in population
   for(size_t i=0;i<pop_num;i++)
   {
      //clear all content of DecValueOfGenome
      DecValueOfGenome.clear();
      // get dec. Value of ith Genome.
      ((*CurrentPopulation)[i])->getDecValue(DecValueOfGenome);
      // Run evaluator, obtained raw value.
      tempOrigValue = m_pEvaluator->CalcuRawFit(DecValueOfGenome,i,runstate);

      OrigScore.push_back(tempOrigValue);
      ((*CurrentPopulation)[i])->setOrigValue( OrigScore[i] );
      ((*CurrentPopulation)[i])->setFinishState(runstate);
   }
      // call the function of m_pEvaluator for converting the original value to Raw score.
   m_pEvaluator->ConvOrigToRawScore(OrigScore);

   for(size_t i=0;i<pop_num;i++)
      ((*CurrentPopulation)[i])->setRawScore(OrigScore[i]);
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
