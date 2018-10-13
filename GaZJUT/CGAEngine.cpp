#include "CGAEngine.h"
#include "CSelector.h"
#include "CEvaluator.h"
#include "CCross.h"
#include "CElist.h"
#include "CMutator.h"
#include "CEvaluator.h"
#include "CGaOperatorBase.h"
#include "CCalcFitnessInterface.h"
#include "CGaparameter.h"
#include "../CalcZJUT/CParameter.h"
#include "CGpopulation.h"
//
#include "../CalcZJUT/CCalcVASP.h"
#include "../CalcZJUT/CCalcGaussian.h"
#include "../CalcZJUT/CCalcDMol.h"
#include "../CalcZJUT/CCalcLammps.h"
#include "../CalcZJUT/CCalcStructBasePool.h"
#include "../CalcZJUT/CCalcClusterStructPool.h"
#include "../CalcZJUT/CCalcSupportStructPool.h"

//
namespace GAZJUT{

CGAEngine(CALCZJUT::CParameter* para)
:m_pParameter(para)
{
    m_pGaparameter = m_pParameter->GaParameter();
}

CGAEngine::~CGAEngine()
{
    delete m_pFitnessCalculator;
    delete m_pCurrentPopulation;
    delete this->m_pStructurePool;
    for(size_t i=0;i<m_GeneticOperator.size();i++)
          delete m_GeneticOperator[i];
    m_GeneticOperator.clear();
}
/*
1. According the inputting setting, construct Fitness calculator (such as VASP, GAUSSSIAN, DMOL, LAMMPS )
2. According the inputting setting, construct Calculation modes(such as pure cluster, adsorbent/cluster, adsorbent/2DSupport )
3. Construct population,initialize genome and gene.
4. Construct GAoperators: selector, crosser, mutator, elitor,
*/
void CGAEngine::init()
{
   if((*m_pGaparameter)["Evaluator"]=="VASP")
       m_pFitnessCalculator = new CALCZJUT::CCalcVASP(this->m_pParameter);
   else if((*m_pGaparameter)["Evaluator"]=="DMOL")
       m_pFitnessCalculator = new CALCZJUT::CCalcDMol(this->m_pParameter);
   else if((*m_pGaparameter)["Evaluator"]=="GAUSSIAN")
       m_pFitnessCalculator = new CALCZJUT::CCalcGaussian(this->m_pParameter);
   else if((*m_pGaparameter)["Evaluator"]=="LAMMPS")
       m_pFitnessCalculator = new CALCZJUT::CCalcLammps(this->m_pParameter);

   switch ((int)m_pParameter->simulationMode)
   {
       case CParameter::CLUSTER:
           m_pStructurePool = new CCalcClusterStructPool(this->m_pParameter);
       case CParameter::PERIODIC:
           break;
       case CParameter::MOL_2DMATERIAL:
       case CParameter::MOL_CLUSTER:
           m_pStructurePool = new CCalcSupportStructPool(this->m_pParameter);
           break;
       default:
           break;
   }
   /*
     Initialize structural pool( read structural file or random identify it)
     Obtain the gene-variable range for the construction of Population object
   */
   this->m_pStructurePool->init();
   //
   this->m_pFitnessCalculator->init();
   // until now, all parameters in object of Gaparameter were set.
   //
   m_pCurrentPopulation = new CGpopulation(m_pGaparameter);

   // sequence of operators!
   m_GeneticOperator.push_back(new CEvaluator());
   m_GeneticOperator.push_back(new CFitnessScaling());
   m_GeneticOperator.push_back(new CElist());
   m_GeneticOperator.push_back(new CSelector());
   m_GeneticOperator.push_back(new CCross());
   m_GeneticOperator.push_back(new CMutator());
}
void CGAEngine::evolve()
{
   // read the setting generation number
   size_t total_pop_num=m_pGaparameter->GenerationNum();

   while(true){

      for(size_t i=0;i<m_GeneticOperator.size();i++)
          m_GeneticOperator[i]->run(m_pCurrentPopulation);
      m_pGaparameter->add_Curr_Generation();
      if(m_pGaparameter->Curr_Generation>=total_pop_num)
         break;
   }
}

}
