#include "CAIEngineBase.h"

CAIEngineBase::CAIEngineBase(CALCZJUT::CParameter* mbf)
:m_pParameter(mbf)
{
    this->m_pStructurePool=nullptr;
}

CAIEngineBase::~CAIEngineBase()
{
    //dtor
}
