#include "CStructPoolPeriodic.h"
#include "CModelPeriodicStruct.h"
namespace CALCZJUT{

CStructPoolPeriodic::CStructPoolPeriodic(CParameter* othr)
:CStructPoolBase(othr)
{
    for(size_t i=0;i<this->m_pParameter->GaParameter()->PopNum();i++){
        this->m_CalcStructPool.push_back(new CModelPeriodicStruct(this->m_pParameter,i));
        m_CalcStructPool[m_CalcStructPool.size()-1]->periodicFramework()->setExcludeBond(m_pParameter->excludeBond);
        m_CalcStructPool[m_CalcStructPool.size()-1]->periodicFramework()->setTolerancefactor(m_pParameter->bondToleranceFactor);
    }
}

CStructPoolPeriodic::~CStructPoolPeriodic()
{

}


}
