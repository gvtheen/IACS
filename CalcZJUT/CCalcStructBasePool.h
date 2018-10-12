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
        virtual void GeneVARRange(std::vector<GeneVAR>&);
        CCalcModeStruct* operator[](size index);


    public:
        std::vector<CCalcModeStruct*>      m_CalcStructPool;
                           CParameter     *m_pParameter;
                                 Cios     *m_IO;
                 std::vector<GeneVAR>     *m_pGeneVAR;


};




}//NAMESPACE
#endif // CCALCSTRUCTUREBASEPOOL_H
