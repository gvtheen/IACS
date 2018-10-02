#ifndef CCALCMODESTRUCT_H
#define CCALCMODESTRUCT_H
#include <vector>
#include "../Util/Bitset.h"
#include "../GaZJUT/GaDeclaration.h"
//
using util::Bitset;
using GAZJUT::GENEVAR;

namespace CATAZJUT{
  class CPeriodicFramework;
}
namespace CALCZJUT{
    //declare this class in this namespace
class CParameter;

class CCalcModeStruct
{
    public:
        CCalcModeStruct(CParameter*);
        virtual ~CCalcModeStruct();

        virtual void setGeneValueToStruct(const std::vector<double>& realValueOfgene)=0;
        virtual std::vector<double>*  getGeneValuefromStruct()const =0;
        virtual std::vector<GENEVAR>* GeneVarRange()=0;
        // only effective for supported catalyst
        virtual void createSupport( const Bitset &);
        virtual void createMoleAdsorb( const Bitset &);

        virtual void init();
        virtual void createStructureAtGene();
        virtual void removeStructureOfGene();
        CATAZJUT::CPeriodicFramework* periodicFramework();
        CATAZJUT::CPeriodicFramework* periodicFramework(size_t index_int);
        void setPeriodicFramekwork(CATAZJUT::CPeriodicFramework*);

    public:
                          CParameter*          m_pParameter;
                std::vector<GENEVAR>*          m_pGeneVar;

        CATAZJUT::CPeriodicFramework*          m_pPeriodicFramework;
std::vector<CATAZJUT::CPeriodicFramework*>     m_PopuPeriodicFramework;
        CATAZJUT::CPeriodicFramework*          m_backupPeriodicFramework;

};


}
#endif // CCALCMODESTRUCT_H
