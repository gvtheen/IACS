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
#include "CGAEngine.h"
#include "CSelector.h"
#include "CEvaluator.h"
#include "CCross.h"
#include "CElist.h"
#include "CMutator.h"
#include "CEvaluator.h"
#include "CFitnessScaling.h"
#include "CGaOperatorBase.h"
#include "CExeFitnessInterface.h"
#include "CGaparameter.h"
#include "../CalcZJUT/CParameter.h"
#include "CGpopulation.h"
#include "../CalcZJUT/CExeVASP.h"
#include "../CalcZJUT/CExeGaussian.h"
#include "../CalcZJUT/CExeDMol.h"
#include "../CalcZJUT/CExeLammps.h"
#include "../CalcZJUT/CStructPoolBase.h"
#include "../CalcZJUT/CStructPoolCluster.h"
#include "../CalcZJUT/CStructPoolSupported.h"

using CALCZJUT::CParameter;
namespace GAZJUT{

CGAEngine::CGAEngine(CALCZJUT::CParameter* para)
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
       m_pFitnessCalculator = new CALCZJUT::CExeVASP(this->m_pParameter);
   else if((*m_pGaparameter)["Evaluator"]=="DMOL")
       m_pFitnessCalculator = new CALCZJUT::CExeDMol(this->m_pParameter);
   else if((*m_pGaparameter)["Evaluator"]=="GAUSSIAN")
       m_pFitnessCalculator = new CALCZJUT::CExeGaussian(this->m_pParameter);
   else if((*m_pGaparameter)["Evaluator"]=="LAMMPS")
       m_pFitnessCalculator = new CALCZJUT::CExeLammps(this->m_pParameter);

   switch ((int)m_pParameter->simulationMode)
   {
       case CParameter::CLUSTER:
           m_pStructurePool = new CALCZJUT::CStructPoolCluster(this->m_pParameter);
           break;
       case CParameter::PERIODIC:
           break;
       case CParameter::MOL_2DMATERIAL:                 // It is just the same as that of MOL_CLUSTER.
       case CParameter::MOL_CLUSTER2DMATERIAL:
       case CParameter::MOL_CLUSTER:
           m_pStructurePool = new CALCZJUT::CStructPoolSupported(this->m_pParameter);
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
