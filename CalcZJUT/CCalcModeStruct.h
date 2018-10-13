#ifndef CCALCMODESTRUCT_H
#define CCALCMODESTRUCT_H
#include <vector>
#include "../Util/Bitset.h"
#include "../GaZJUT/GaDeclaration.h"
//
using util::Bitset;
using GAZJUT::GeneVAR;

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

        virtual CCalcModeStruct* clone()=0;

        virtual void setGeneValueToStruct(const std::vector<double>& realValueOfgene)=0;
        virtual void getGeneValuefromStruct(std::vector<double>&) =0;
        virtual void GeneVARRange(std::vector<GeneVAR>&)=0;
        // only effective for supported catalyst
        virtual void createSupport( const Bitset &);
        virtual void createMoleAdsorb( const Bitset &);

        virtual void init();

        CATAZJUT::CPeriodicFramework* periodicFramework();
        void setPeriodicFramekwork(CATAZJUT::CPeriodicFramework*);

        void setRandomInitState(const bool&);
        bool RandomInitState();

        virtual std::vector<std::pair<std::string,size_t>>& chemicalFormula();
        virtual void setChemicalFormula(const std::vector<std::pair<std::string,size_t>>&);

    public:
                          CParameter*          m_pParameter;
                std::vector<GeneVAR>*          m_pGeneVAR;

        CATAZJUT::CPeriodicFramework*          m_pPeriodicFramework;
        CATAZJUT::CPeriodicFramework*          m_backupPeriodicFramework;

 std::vector<std::pair<std::string,size_t>>    m_chemicalFormula;
                                 bool          m_IsNeedRandomInit;

};


}
#endif // CCALCMODESTRUCT_H
