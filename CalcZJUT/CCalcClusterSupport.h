#ifndef CCALCCLUSTERSUPPORT_H
#define CCALCCLUSTERSUPPORT_H
#include <Eigen/Dense>
#include "../Util/Point-Vector.h"
#include "CCalcModeStruct.h"
#include "../GaZJUT/GaDeclaration.h"
#include "../Util/Bitset.h"

using GAZJUT::GENEVAR;
using util::Bitset;
namespace CATAZJUT{
   class CCrystalPlanes;
}
namespace CALCZJUT{

class CCalcSupportBase;
class CMoleculeAdsorbent;
class CCalcMoleculeAdsorbent;
class CParameter;

class CCalcClusterSupport:public CCalcModeStruct
{
    public:
        CCalcClusterSupport(CParameter*);
        virtual ~CCalcClusterSupport();

        //virtual function from CCalcModeStruct
        void setGeneValueToStruct(const std::vector<double>& realValueOfgene);
        std::vector<double>* getGeneValuefromStruct()const;
        std::vector<GENEVAR>* GeneVarRange();
        void backUpStructure();
        void fromBackupToCurrent();

        void createSupport(Bitset &);
        void createMoleAdsorb(Bitset &);
    protected:

    private:
          CCalcSupportBase*      m_pSupport;
        CMoleculeAdsorbent*      m_pAdsorbMolecule;
  CATAZJUT::CCrystalPlanes*      m_pCrystalPlanes;
                     Bitset      m_BitbackupSupport;
                     Bitset      m_BitbackupAdsorbMolecule;



};



}
#endif // CCalcClusterSupport_H
