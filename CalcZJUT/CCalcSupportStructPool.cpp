#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <strings.h>
#include <cstring>
#include <math.h>
#include "../CataZJUT/CPeriodicFramework.h"
#include "CCalcCluster.h"
#include "CParameter.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CCartesianCoordinates.h"
#include "../Util/log.hpp"
#include "../Util/foreach.h"
#include "../Util/Bitset.h"
#include "../Util/Point-Vector.h"
#include "../CataZJUT/Constant.h"
#include "../CataZJUT/CFragment.h"
#include "../Util/utilFunction.h"
#include "../Util/CRandomgenerator.h"
#include "../GaZJUT/CGaparameter.h"
#include "../CataZJUT/CElement.h"
#include "Cios.h"
#include "CIOMol.h"
#include "CIOCar.h"
#include "CIOGjf.h"
#include "CIOPoscar.h"
#include "CCalcSupportStructPool.h"


namespace CALCZJUT{

CCalcSupportStructPool::CCalcSupportStructPool(CParameter* oth)
:CCalcStructureBasePool(oth)
{
    for(size_t i=0;i<this->m_pParameter->GaParameter()->PopNum();i++){
        if(m_pParameter->simulationMode ==CParameter::MOL_2DMATERIAL)
            this->m_CalcStructPool.push_back(new CCalc2DSupport(this->m_pParameter));
         else
            this->m_CalcStructPool.push_back(new CCalcClusterSupport(this->m_pParameter));
         m_CalcStructPool[m_CalcStructPool.size()-1]->periodicFramework()->setExcludeBond(m_Parameter->excludeBond);
         m_CalcStructPool[m_CalcStructPool.size()-1]->periodicFramework()->setTolerancefactor(m_Parameter->bondToleranceFactor);
    }
}

CCalcSupportStructPool::~CCalcSupportStructPool()
{
    //dtor
}
void CCalcSupportStructPool::init()
{
       if(this->m_pParameter->adso_supp_Input_File.size()!= 0){
           //read coordinate of mixed adso_supp_Struct, add the pointer of m_pCalcModeStruct;
           CATAZJUT::CFragment *m1,*m2;
           for(size_t i=0;i<m_pParameter->adso_supp_Input_File.size();i++)
           {
               if(i >= m_CalcStructPool.size())
                  break;

               this->getIO(m_pParameter->adso_supp_Input_File[i],m_CalcStructPool[i]->periodicFramework())
               this->m_IO->input(m_pParameter->adso_supp_Input_File[i]);
                m_CalcStructPool[i]->periodicFramework()->perceiveBonds();
                m_CalcStructPool[i]->periodicFramework()->perceiveFragments();
                m1 = m_CalcStructPool[i]->periodicFramework()->fragment(0);
                m2 = m_CalcStructPool[i]->periodicFramework()->fragment(1);
               if(m1== nullptr || m2 ==nullptr){
                  Log::Error<<"Support and adsorbent setting of input file are error! CCalcSupportStructPool::init().\n";
                  boost::throw_exception(std::runtime_error("Support and adsorbent setting are error! Check the file: CCalcSupportStructPool::init()."));
               }
               if(m1->atomCount() > m2->atomCount())
               {
                  m_CalcStructPool[i]->createSupport(m1->bitSet());
                  m_CalcStructPool[i]->createMoleAdsorb(m2->bitSet());
               }else{
                  m_CalcStructPool[i]->createSupport(m2->bitSet());
                  m_CalcStructPool[i]->createMoleAdsorb(m1->bitSet());
               }
               m_CalcStructPool[i]->setRandomInitState(false);
           }
      }else{ // treat cluster-support
            //read coordinate of mixed adso_supp_Struct, add the pointer of m_pCalcModeStruct;
           //
           this->getIO(m_Parameter->supportStructFile,m_CalcStructPool[0]->periodicFramework());
           this->m_IO->input(m_Parameter->supportStructFile);


           this->getIO(m_Parameter->adsorbentStructFile,m_CalcStructPool[0]->periodicFramework());
           Bitset Bit_adsorbent = this->m_IO->input(m_Parameter->adsorbentStructFile,CParameter::MOL_CLUSTER);
           // get opposite bit for support
           m_CalcStructPool[0]->createSupport(~Bit_adsorbent);
           //set molecular adsorbent bits
           m_CalcStructPool[0]->createMoleAdsorb(Bit_adsorbent);
           m_CalcStructPool[0]->setRandomInitState(false);
     }
}
void CCalcSupportStructPool::getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework)
{
    std::vector<std::string> vectstr;
    boost::algorithm::split(vectstr,file_name,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
    boost::algorithm::trim(vectstr[1]);
    if(vectstr[1]=="mol")
        this->m_IO = new CIOMol(currentPeriodicFramework);
    else if(vectstr[1]=="car")
        this->m_IO = new CIOCar(currentPeriodicFramework);
    else if(vectstr[1]=="poscar")
        this->m_IO = new CIOPoscar(currentPeriodicFramework);
    else if(vectstr[1]=="gjf")
        this->m_IO = new CIOGjf(currentPeriodicFramework);
    else if(vectstr[1]=="cif")
        this->m_IO = new CIOCif(currentPeriodicFramework);
    else if(vectstr[1]=="cell")
        this->m_IO = new CIOCellFile(currentPeriodicFramework);
    else
        ;

    return nullptr;
}

}//namespace