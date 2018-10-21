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
#include "CEvaluator.h"
#include "../CalcZJUT/CExeFitnessInterface.h"
#include "../CalcZJUT/CStructPoolBase.h"
#include "../CalcZJUT/CModelBase.h"
#include "../CalcZJUT/CParameter.h"

namespace CALCZJUT{
   class CStructPoolBase;
}

namespace GAZJUT{

CEvaluator::CEvaluator()
{

}
CEvaluator::CEvaluator(CALCZJUT::CExeFitnessInterface  *myEvaluator,
                       CALCZJUT::CStructPoolBase *myStructPool)
{
    // sample evaluator
    this->m_pEvaluator     = myEvaluator;
    //sample structural pool
    this->m_pStructurePool = myStructPool;
    // sample evaluator monitor!
    pop_run_state.resize(m_pStructurePool->m_pParameter->popNum(),false);
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

      //sett i th structure to evaluator
      m_pEvaluator->setCalcModeStruct((*m_pStructurePool)[i]);

      #ifdef DEBUG
         (*m_pStructurePool)[i]->outputStructureToFile();
      #endif // DEBUG

      // Run evaluator, obtained raw value.
      tempOrigValue = m_pEvaluator->CalcuRawFit(DecValueOfGenome,i,runstate);

      // Get dec value from relaxed structure after calculation
      // Re-set this value to individual Genome.
      DecValueOfGenome.clear();
      (*m_pStructurePool)[i]->getGeneValuefromStruct(DecValueOfGenome);

      this->pop_run_state.set(i,true);
      //if some updating gene was gotten, set them to population
      if(DecValueOfGenome.size()!=0)
         ((*CurrentPopulation)[i])->updateDecValueGene(DecValueOfGenome);

      OrigScore.push_back(tempOrigValue);
      ((*CurrentPopulation)[i])->setOrigValue( OrigScore[i] );
      ((*CurrentPopulation)[i])->setFinishState(runstate);
   }
      // call the function of m_pEvaluator for converting the original value to Raw score.
   m_pEvaluator->ConvOrigToRawScore(OrigScore);

   for(size_t i=0;i<pop_num;i++)
      ((*CurrentPopulation)[i])->setRawScore(OrigScore[i]);
}
void CEvaluator::setCalcFitnessInterface(CALCZJUT::CExeFitnessInterface* CalcFitness)
{
      this->m_pEvaluator=CalcFitness;
}
CALCZJUT::CExeFitnessInterface* CEvaluator::CalcFitnessInterface()
{
      return this->m_pEvaluator;
}




}
