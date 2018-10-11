#include <boost/algorithm/string.hpp>
#include "CCalcFitnessInterface.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "CParameter.h"
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
using util::Bitset;

namespace CALCZJUT{

CCalcFitnessInterface::CCalcFitnessInterface(CParameter* mpara)
:m_Parameter(mpara)
{

}

CCalcFitnessInterface::~CCalcFitnessInterface()
{
}

void CCalcFitnessInterface::init()
{
     if(m_Parameter->simulationMode ==CParameter::MOL_2DMATERIAL || m_Parameter->simulationMode ==MOL_CLUSTER )
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
           Cios* tempIO = std::move(getIO(m_Parameter->adso_supp_Struct,m_pCalcModeStruct->periodicFramework()));
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
           Cios* tempIO = std::move(getIO(m_Parameter->supportStructFile,m_pCalcModeStruct->periodicFramework()));
           tempIO->input(m_Parameter->supportStructFile);
           delete tempIO;

           Cios* tempIO_1 = std::move(getIO(m_Parameter->adsorbentStructFile,m_pCalcModeStruct->periodicFramework()));
           Bitset Bit_adsorbent = tempIO_1->input(m_Parameter->adsorbentStructFile,CParameter::MOL_CLUSTER);
           delete tempIO_1;
           // get opposite bit for support
           m_pCalcModeStruct->createSupport(~Bit_adsorbent);
           //set molecular adsorbent bits
           m_pCalcModeStruct->createMoleAdsorb(Bit_adsorbent);
         }
     }else{    //pure cluster model
        CCalcCluster* cluster = new CCalcCluster(this->m_Parameter);
        //m_pIO = new CIOPoscar( cluster->periodicFramework());
        this->m_pCalcModeStruct = cluster;
        // construct initial structures
        this->m_pCalcModeStruct->init();
        //m_pIO->input("POSCAR");
     }

    //initialize gene varible
    //
    m_pCalcModeStruct->GeneVARRange(m_Parameter->GaParameter()->GeneVAR());

    m_pCalcModeStruct->periodicFramework()->m_pBondEvaluator->setExcludeBond(m_Parameter->excludeBond);

    m_pCalcModeStruct->periodicFramework()->m_pBondEvaluator->setTolerancefactor(m_Parameter->bondToleranceFactor);
}
void CCalcFitnessInterface::GetDecGeneAfterCalc(std::vector<double>& tmpVect)
{
    return this->m_pCalcModeStruct->getGeneValuefromStruct(tmpVect);
}
double CCalcFitnessInterface::CalcuRawFit(std::vector<double>& RealValueOfGenome,size_t& pop_index, bool& isNormalexist)
{
    return 0;
}

void CCalcFitnessInterface::ConvOrigToRawScore(std::vector<double>& othr)
{

}
void CCalcFitnessInterface::setCalcModeStruct(CCalcModeStruct* Temp_calcModeStruct)
{
    this->m_pCalcModeStruct = Temp_calcModeStruct;
}
CCalcModeStruct* CCalcFitnessInterface::calcModeStruct()
{
    return this->m_pCalcModeStruct;
}
void CCalcFitnessInterface::setIO(Cios* m_IO)
{
    this->m_pIO=m_IO;
}
Cios* CCalcFitnessInterface::IO()const
{
    return this->m_pIO;
}
Cios* getIO(std::string &file_name,CATAZJUT::CPeriodicFramework* currentPeriodicFramework)
{
    std::vector<std::string> vectstr;
    boost::algorithm::split(vectstr,file_name,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
    boost::algorithm::trim(vectstr[1]);
    if(vectstr[1]=="mol")
        return new CIOMol(currentPeriodicFramework);
    else if(vectstr[1]=="car")
        return new CIOCar(currentPeriodicFramework);
    else if(vectstr[1]=="poscar")
        return new CIOPoscar(currentPeriodicFramework);
    else if(vectstr[1]=="gjf")
        return new CIOGjf(currentPeriodicFramework);
    else if(vectstr[1]=="cif")
        return new CIOCif(currentPeriodicFramework);
    else if(vectstr[1]=="cell")
        return new CIOCellFile(currentPeriodicFramework);
    else
        ;

    return nullptr;
}

}
