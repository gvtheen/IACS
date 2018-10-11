#include "CCalcSupportStructPool.h"


namespace CALCZJUT{

CCalcSupportStructPool::CCalcSupportStructPool(CParameter* oth)
:CCalcStructureBasePool(oth)
{
    for(size_t i=0;i<this->m_pParameter->GaParameter()->PopNum();i++)
        if(m_pParameter->simulationMode ==CParameter::MOL_2DMATERIAL){
            this->m_CalcStructPool.push_back(new CCalc2DSupport(this->m_pParameter));
        }else{
            this->m_CalcStructPool.push_back(new CCalcClusterSupport(this->m_pParameter));
        }
}

CCalcSupportStructPool::~CCalcSupportStructPool()
{
    //dtor
}
void CCalcSupportStructPool::init()
{
       if(m_pParameter->adso_supp_Struct != ""){
           //read coordinate of mixed adso_supp_Struct, add the pointer of m_pCalcModeStruct;
           //
           Cios* tempIO = this->getIO(m_Parameter->adso_supp_Struct,m_pCalcModeStruct->periodicFramework());
           tempIO->input(m_Parameter->adso_supp_Struct);
           delete tempIO;
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
           Cios* tempIO = this->getIO(m_Parameter->supportStructFile,m_pCalcModeStruct->periodicFramework());
           tempIO->input(m_Parameter->supportStructFile);
           delete tempIO;

           Cios* tempIO_1 = this->getIO(m_Parameter->adsorbentStructFile,m_pCalcModeStruct->periodicFramework());
           Bitset Bit_adsorbent = tempIO_1->input(m_Parameter->adsorbentStructFile,CParameter::MOL_CLUSTER);
           delete tempIO_1;
           // get opposite bit for support
           m_pCalcModeStruct->createSupport(~Bit_adsorbent);
           //set molecular adsorbent bits
           m_pCalcModeStruct->createMoleAdsorb(Bit_adsorbent);
         }
}

}//namespace
