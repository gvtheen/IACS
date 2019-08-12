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
#include "CCataCalcAdapter.h"

#include "../Util/log.hpp"

using util::Log;

namespace CALCZJUT{

CCataCalcAdapter::CCataCalcAdapter(CALCZJUT::CParameter* _Cparameter)
:CAdapterBase(),m_pParameter( _Cparameter )
{
   Log::Info<<"Initialize Catalyst-simulation Adapter...."<<std::endl;

   std::vector<std::string> res;

   m_pParameter->GetEvaluateEXE(res);

   for(size_t i=0;i<res.size();i++){
       if(res[i]=="VASP")
           m_FitnessCalculatorPool.push_back(new CALCZJUT::CExeVASP(this->m_pParameter));
       else if(res[i]=="DMOL")
           m_FitnessCalculatorPool.push_back(new CALCZJUT::CExeDMol(this->m_pParameter));
       else if(res[i]=="GAUSSIAN")
           m_FitnessCalculatorPool.push_back(new CALCZJUT::CExeGaussian(this->m_pParameter));
       else if(res[i]=="LAMMPS")
           m_FitnessCalculatorPool.push_back(new CALCZJUT::CExeLammps(this->m_pParameter));
       else if(res[i]=="DFTB")
           m_FitnessCalculatorPool.push_back(new CALCZJUT::CExeDFTB(this->m_pParameter));
       else if(res[i]=="CASTEP")
           m_FitnessCalculatorPool.push_back(new CALCZJUT::CExeCastep(this->m_pParameter));
   }
   if(this->m_FitnessCalculatorPool.size()==0){
      Log::Error<< "No fitness calculator!! CCataCalcAdapter::CCataCalcAdapter!\n";
      boost::throw_exception(std::runtime_error("No fitness calculator!!! CCataCalcAdapter::CCataCalcAdapter!\n"));
   }

   Log::Info<<"Initialize Structural Pool....."<< std::endl;

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

   this->m_pStructurePool->init();

   Log::Info<<"Initialize Fitness calculator Pool....."<< std::endl;

   for(size_t i=0;i<m_FitnessCalculatorPool.size();i++)
        m_FitnessCalculatorPool[i]->init();
}
void CCataCalcAdapter::computationToEngine(std::vector<IACSZJUT::VarRangeStruct>& _mvalue)
{
    assert(this->m_pStructurePool);

    this->m_pStructurePool->VarRangeStructRange( _mvalue);
}
void CCataCalcAdapter::computationToEngine(std::vector<double>& _mvalue)
{
    this->m_currentEvaluator->
}
void CCataCalcAdapter::engineToComputation(std::vector<double>& _input,
                                           size_t _index,
                                           bool& _state,
                                           double& _output)
{
   //CALCZJUT::CExeFitnessInterface *currentEvaluator = this->m_FitnessCalculator[0];
   Log::Info<<"Run "<<m_currentEvaluator->ExeName()<<" calculation of No "<<_index+1;

   double _output = m_currentEvaluator->CalcuRawFit(_mvalue,_index,_state);
}
double CCataCalcAdapter::convertOrigValueToRawvalue(const double& _mOriValue)
{
   swithch((int)(m_pParameter->searchMode)){
       case CParameter::MIN:
           return
           break;
       case CParameter::MAX:
           return
           break;
       case CParameter::SPE:
           return std::fabs(m_pParameter->optimal_specific_value - _mOriValue);
           break;
       default:
   }

}

}

CCataCalcAdapter::~CCataCalcAdapter()
{
    //dtor
}


}
