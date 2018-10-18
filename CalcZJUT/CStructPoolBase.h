#ifndef CCALCSTRUCTUREBASEPOOL_H
#define CCALCSTRUCTUREBASEPOOL_H
#include <vector>

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
        virtual void GeneVARRange(std::vector<GeneVAR>&);
        CModelBase* operator[](size index);


    public:
        std::vector<CModelBase*>      m_CalcStructPool;
                     CParameter      *m_pParameter;
                         CIOBase     *m_IO;
            std::vector<GeneVAR>     *m_pGeneVAR;


};




}//NAMESPACE
#endif // CCALCSTRUCTUREBASEPOOL_H
