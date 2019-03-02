/******************************************************************************
**
** Copyright (C) 2019-2031 Dr.Gui-lin Zhuang <glzhuang@zjut.edu.cn>
** All rights reserved.
**
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/
#include <float.h>
#include "CEvaluator.h"
#include "CGenome.h"
#include "../CalcZJUT/CExeFitnessInterface.h"
#include "../CalcZJUT/CStructPoolBase.h"
#include "../CalcZJUT/CModelBase.h"
#include "../CalcZJUT/CParameter.h"
#include "../Util/Bitset.h"
#include "../Util/log.hpp"


using util::Bitset;
using util::Log;

namespace CALCZJUT{
   class CStructPoolBase;
   class CParameter;
}

namespace GAZJUT{

CEvaluator::CEvaluator()
{

}
CEvaluator::CEvaluator(std::vector<CALCZJUT::CExeFitnessInterface*>  myEvaluatorPool,
                       CALCZJUT::CStructPoolBase *myStructPool)
{
    // sample evaluator
    this->m_pEvaluatorPool.assign(myEvaluatorPool.begin(),myEvaluatorPool.end());
    this->m_currentEvaluator=nullptr;
    // sample structural pool
    this->m_pStructurePool = myStructPool;
    // sample evaluator monitor!
    pop_run_state.resize(m_pStructurePool->m_pParameter->popNum(),false);
}
CEvaluator::~CEvaluator()
{
}
void CEvaluator::run(CGpopulation* CurrentPopulation)
{
   assert(CurrentPopulation);

   #ifdef DEBUG
      Log::Debug<<"*********** CEvaluator::run***********"<< std::endl;
   #endif

   bool runstate=false;
   size_t pop_num = CurrentPopulation->popNum();
    #ifdef DEBUG
         Log::Debug<<"1-*********** CEvaluator::run***********"<< std::endl;
    #endif
   std::vector<double> OrigScoreVect;
   std::vector<double> DecValueOfGenome;

   std::map <size_t, double> OrigScoreMap;
   // calculate the fitness of each genome in population

   //put the 1st evaluator into current evaluator.
   assert(m_pEvaluatorPool[0]);
   m_currentEvaluator= m_pEvaluatorPool[0];
   #ifdef DEBUG
         Log::Debug<<"2-*********** CEvaluator::run***********"<< std::endl;
    #endif
   for(size_t i=0;i<pop_num;i++)
   {
      //clear all content of DecValueOfGenome
      DecValueOfGenome.clear();
      // get dec. Value of ith Genome.
      #ifdef DEBUG
         Log::Debug<<"3-*********** CEvaluator::run***********"<< std::endl;
      #endif
      ((*CurrentPopulation)[i])->getDecValue(DecValueOfGenome);

      //sett i th structure to evaluator
      m_currentEvaluator->setCalcModeStruct((*m_pStructurePool)[i]);

      #ifdef DEBUG
         (*m_pStructurePool)[i]->outputStructureToFile();
      #endif // DEBUG

      // Run evaluator, obtained raw value.
      OrigScoreMap[i] = m_currentEvaluator->CalcuRawFit(DecValueOfGenome,i,runstate);
      /** \brief
       *  \OrigScoreMap:    build the map relation between pop index and orig. Value;
       */

      // Get dec value from relaxed structure after calculation
      // Re-set this value to individual Genome.
      DecValueOfGenome.clear();
      (*m_pStructurePool)[i]->getGeneValuefromStruct(DecValueOfGenome);

      this->pop_run_state.set(i,true);
      //if some updating gene was gotten, set them to population
      if(DecValueOfGenome.size()!=0)
         ((*CurrentPopulation)[i])->updateDecValueGene(DecValueOfGenome);

      ((*CurrentPopulation)[i])->setOrigValue( OrigScoreMap[i] );
      ((*CurrentPopulation)[i])->setFinishState(runstate);
   }
   // convert original value to raw score by using specific method.
   // call the function of currentEvaluator for converting the original value to Raw score.
   // Then set the set value to the whole population.
   OrigScoreVect.resize(pop_num);
   for(size_t i=0;i<pop_num;i++)
      OrigScoreVect[i] = OrigScoreMap[i];
   m_currentEvaluator->ConvOrigToRawScore(OrigScoreVect);
   for(size_t i=0;i<pop_num;i++)
      ((*CurrentPopulation)[i])->setRawScore(OrigScoreVect[i]);
   // After one cycle calculation, output some structure and computational value.
   // output the result
   this->standardOutput(CurrentPopulation);
   // output the relaxed structure
   this->standardOutput(OrigScoreMap);

}
void CEvaluator::standardOutput(std::map <size_t, double>& mapIndexScore)
{
   std::vector<size_t> res;

   switch((int)(m_pStructurePool->m_pParameter->evaluatorCriterion)){
          case CALCZJUT::CParameter::ENERGY:
          case CALCZJUT::CParameter::FORCE:
               getTargetPopWithCondition(res, mapIndexScore, 1, 0);
               break;
          case CALCZJUT::CParameter::BAND_GAP:
               getTargetPopWithCondition(res, mapIndexScore, 1,
                                         m_pStructurePool->m_pParameter->optimal_gap_value);
               break;
          default:
               break;
   }
   for(size_t i=0;i<res.size();i++)
       (*m_pStructurePool)[i]->standardOutput(1);
}
void CEvaluator::standardOutput( CGpopulation* CurrentPopulation )
{
   assert(m_currentEvaluator);
   assert(m_pStructurePool);

   Log::Output<<"Computational results in "<<m_pStructurePool->m_pParameter->currentGenerationNum();
   Log::Output<<" th generation by using " << this->m_currentEvaluator->ExeName()<<std::endl;
   Log::Output<<"-------------------------------------------------------------------"<<std::endl;
   for(size_t i=0;i<CurrentPopulation->popNum();i++){
       switch((int)(m_pStructurePool->m_pParameter->evaluatorCriterion)){
          case CALCZJUT::CParameter::ENERGY:
               Log::Output<<"  E( Pop["<<i+1<<"] ) = ";
               break;
          case CALCZJUT::CParameter::FORCE:
               Log::Output<<"  F( Pop["<<i+1<<"] ) = ";
               break;
          case CALCZJUT::CParameter::BAND_GAP:
               Log::Output<<"Gap( Pop["<<i+1<<"] ) = ";
          default:
               break;
       }
       Log::Output<<(*(*CurrentPopulation)[i])["origvalue"]<<std::endl;
   }
   Log::Output<<"-------------------------------------------------------------------"<<std::endl;
}
//void CEvaluator::setCalcFitnessInterface(CALCZJUT::CExeFitnessInterface* CalcFitness)
//{
//    this->m_pEvaluator=CalcFitness;
//}
//CALCZJUT::CExeFitnessInterface* CEvaluator::CalcFitnessInterface()
//{
//    return this->m_pEvaluator;
//}
void CEvaluator::getTargetPopWithCondition(std::vector<size_t>& res, std::map <size_t, double>& mapIndexValue,
                                           size_t num, double conditionValue)
{
    Bitset pos(mapIndexValue.size());
    pos.set(true);

    if(res.size()!=0)
        res.clear();

    double minValue=DBL_MAX;
    size_t pos_index_min=0;
    while(num>0)
    {
        for(size_t i=0; i<mapIndexValue.size(); i++){
          switch((int)(m_pStructurePool->m_pParameter->evaluatorCriterion)){
             case CALCZJUT::CParameter::ENERGY:
             case CALCZJUT::CParameter::FORCE:
                 if(pos[i]==1 && mapIndexValue[i] < minValue ){
                    minValue = mapIndexValue[i];
                    pos_index_min = i;
                 }
                 break;
             case CALCZJUT::CParameter::BAND_GAP:
                 if(pos[i]==1 && std::fabs(mapIndexValue[i]-conditionValue) < minValue ){
                    minValue = std::fabs(mapIndexValue[i]-conditionValue);
                    pos_index_min = i;
                 }
                 break;
             default:
                break;
          }
        }
        res.push_back(pos_index_min);
        pos.set(pos_index_min,false);
        num--;
    }
}



}
