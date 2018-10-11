#include <vector>
#include <boost/algorithm/string.hpp>
#include "../GACatalyst.h"
#include "CInitialization.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../CataZJUT/CFragment.h"
#include "CCalcCluster.h"
#include "CCalcClusterSupport.h"
#include "CCalc2DSupport.h"
#include "Cios.h"
#include "CIOCar.h"
#include "CIOMol.h"
#include "CIOPoscar.h"
#include "CIOGjf.h"
#include "CIOCif.h"
#include "CIOCellFile.h"
#include "../Util/Point-Vector.h"
#include "../GaZJUT/CGaparameter.h"
#include "../Util/log.hpp"
#include "../CataZJUT/CBondTolerance.h"

using util::Log;

namespace CALCZJUT{

CInitialization::CInitialization(CParameter* mpara)
:m_Parameter(mpara)
{
    //ctor
}

CInitialization::~CInitialization()
{
    if(this->m_IO!=nullptr)
        delete this->m_IO;
}
void CInitialization::InitClusterSupport()
{
         if(m_Parameter->simulationMode ==CParameter::MOL_2DMATERIAL){
             this->m_pCalcModeStruct=new CCalc2DSupport(this->m_Parameter);
         }else{
             this->m_pCalcModeStruct=new CCalcClusterSupport(this->m_Parameter);
         }

         if(m_Parameter->adso_supp_Struct != ""){
            // Only one initialized structure is required.
           //read coordinate of mixed adso_supp_Struct, add the pointer of m_pCalcModeStruct;
           //
            this->getIO(m_Parameter->adso_supp_Struct,m_pCalcModeStruct->periodicFramework());
            this->m_IO->input(m_Parameter->adso_supp_Struct);
           //construct chemical bond
            m_pCalcModeStruct->periodicFramework()->perceiveBonds();
           //construct chemical fragments to identify support and adsorbent
            m_pCalcModeStruct->periodicFramework()->perceiveFragments();
           //
           //if(m2dSupport->periodicFramework()->fragments().size())
            CATAZJUT::CFragment *m1,*m2;
            m1 = m_pCalcModeStruct->periodicFramework()->fragment(0);
            m2 = m_pCalcModeStruct->periodicFramework()->fragment(1);
            if(m1== nullptr || m2 ==nullptr){
               Log::Error<<"Support and adsorbent setting of input file are error! input_CCalcVASP.\n";
               boost::throw_exception(std::runtime_error("Support and adsorbent setting are error! Check the file: Error_information.txt."));
            }
            if(m1->atomCount() > m2->atomCount())
            {
                m_pCalcModeStruct->createSupport(m1->bitSet());
                m_pCalcModeStruct->createMoleAdsorb(m2->bitSet());
            }else{
                m_pCalcModeStruct->createSupport(m2->bitSet());
                m_pCalcModeStruct->createMoleAdsorb(m1->bitSet());
            }
         }else{ // treat cluster-support
            //read coordinate of mixed adso_supp_Struct, add the pointer of m_pCalcModeStruct;
           //
            this->getIO(m_Parameter->supportStructFile,m_pCalcModeStruct->periodicFramework());
            this->m_IO->input(m_Parameter->supportStructFile);

            this->getIO(m_Parameter->adsorbentStructFile,m_pCalcModeStruct->periodicFramework());
            Bitset Bit_adsorbent = this->m_IO->input(m_Parameter->adsorbentStructFile,CParameter::MOL_CLUSTER);
           // get opposite bit for support
            m_pCalcModeStruct->createSupport(~Bit_adsorbent);
           //set molecular adsorbent bits
            m_pCalcModeStruct->createMoleAdsorb(Bit_adsorbent);
         }

}
void CInitialization::getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework)
{
    std::vector<std::string> vectstr;
    boost::algorithm::split(vectstr,file_name,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
    boost::algorithm::trim(vectstr[1]);
    if (this->m_IO != nullptr)
        delete this->m_IO;

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
}
void CInitialization::Init2DSupport()
{
    this->InitClusterSupport();
}
void CInitialization::InitCluster()
{
   if(this->m_Parameter->cluster_Input_File.size()!=0){
        this->Initialization(m_Parameter->cluster_Input_File);
    }else if(this->m_Parameter->cluster_Formula!=""){
        this->Initialization(m_Parameter->cluster_Formula);
    }else{
       Log::Error<< " Chemical formula and structural files is required. Initialization in CCalcCluster!\n";
       boost::throw_exception(std::runtime_error("Chemical formula and structural files is required. init_CCalcCluster!!\n"));//ERROR TREATMENT;
    }
}

}//namespace
