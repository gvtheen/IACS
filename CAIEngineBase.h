#ifndef CAIENGINEBASE_H
#define CAIENGINEBASE_H
#include<vector>
#include "IACS.h"

using IACSZJUT::VarRangeStruct;

namespace CALCZJUT{
    class CExeFitnessInterface;
    class CParameter;
    class CStructPoolBase;
}

class CAIEngineBase
{
    public:
        CAIEngineBase(CALCZJUT::CParameter*);
        virtual ~CAIEngineBase();

//        virtual void SetVarRange(const std::vector<VarRangeStruct>&)=0;
//        virtual std::vector<VarRangeStruct>& VarRange()const=0;
		virtual void init()=0;
		virtual void evolve()=0;

    public:
        std::vector<CALCZJUT::CExeFitnessInterface*>  m_FitnessCalculator;
        CALCZJUT::CStructPoolBase                    *m_pStructurePool;
	    CALCZJUT::CParameter                         *m_pParameter;
	    std::vector<IACSZJUT::VarRangeStruct>         m_VarRangeStruct;

    private:

};

#endif // CAIENGINEBASE_H
