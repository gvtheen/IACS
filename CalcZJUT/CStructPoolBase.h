#ifndef CCALCSTRUCTUREBASEPOOL_H
#define CCALCSTRUCTUREBASEPOOL_H
#include <iostream>
#include <vector>
#include "../GaZJUT/GaDeclaration.h"
#include "../IACS.h"

using IACSZJUT::VarRangeStruct;

namespace CALCZJUT{

class CModelBase;
class CParameter;
class CIOBase;

class CStructPoolBase
{
    public:
        CStructPoolBase(CParameter*);
        virtual ~CStructPoolBase();

        virtual void init();
        virtual void VarRangeStructRange(std::vector<VarRangeStruct>&);
        CModelBase* operator[](size_t index);


    public:
        std::vector<CModelBase*>      m_CalcStructPool;
                     CParameter      *m_pParameter;
                         CIOBase     *m_IO;
            std::vector<VarRangeStruct>     *m_pVarRangeStruct;


};




}//NAMESPACE
#endif // CCALCSTRUCTUREBASEPOOL_H
