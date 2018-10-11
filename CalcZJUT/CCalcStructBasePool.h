#ifndef CCALCSTRUCTUREBASEPOOL_H
#define CCALCSTRUCTUREBASEPOOL_H
#include <vector>

namespace CALCZJUT{

class CCalcModeStruct;
class CParameter;
class Cios;

class CCalcStructureBasePool
{
    public:
        CCalcStructureBasePool(CParameter*);
        virtual ~CCalcStructureBasePool();

        virtual void init();
        CCalcModeStruct* operator[](size index);


    public:
        std::vector<CCalcModeStruct*>      m_CalcStructPool;
                           CParameter     *m_pParameter;
                                 Cios     *m_IO;


};




}//NAMESPACE
#endif // CCALCSTRUCTUREBASEPOOL_H
