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


using util::Log;
using util::Point3;
using util::Matrix;

namespace CALCZJUT{

CCalcCluster::CCalcCluster(CParameter* mPara)
:CCalcModeStruct(mPara)
{
    CATAZJUT::CConfigurationBase* temp = new CATAZJUT::CConfigurationBase(mPara);

    for(size_t i=0;i<mPara->GaParameter()->PopNum();i++)
        this->m_PopuPeriodicFramework.push_back(new CATAZJUT::CPeriodicFramework(*temp));

    this->m_pPeriodicFramework = new CATAZJUT::CPeriodicFramework(*temp);
     // In this situation, this->m_pPeriodicFramework->m_pUnitCell is empty pointer NULL;
    delete temp;

}

CCalcCluster::~CCalcCluster()
{
    //dtor
}
void CCalcCluster::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{

}
std::vector<double>*  CCalcCluster::getGeneValuefromStruct()const
{

}
std::vector<GENEVAR>* CCalcCluster::GeneVarRange()
{

}
void CCalcCluster::init()
{
    // perform all initialized work
    if(this->m_pParameter->cluster_Input_File.size()!=0){
        this->Initialization(m_pParameter->cluster_Input_File);
    }else if(this->m_pParameter->cluster_Formula!=""){
        this->Initialization(m_pParameter->cluster_Formula);
    }else{
       Log::Error<< " Chemical formula and structural files is required. init_CCalcCluster!\n";
       boost::throw_exception(std::runtime_error("Chemical formula and structural files is required. init_CCalcCluster!!\n"));//ERROR TREATMENT;
    }

}
void CCalcCluster::Initialization(const std::string& mth)
{
       // initialize from chemical
       // chemical formula = C10 H10 O12
/*
    e.g.:   mth = C10 N S10        // delimiter = blank " "
*/
    //dealwith chemical formula by blank
    std::string str=mth;
    std::vector<std::string> vecStr;
    boost::algorithm::trim(str);
    boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);

    std::pair<std::string,size_t> ChemicalStr;
    for(size_t i=0;i<vecStr.size();i++){
        ChemicalStr.first  = boost::algorithm::trim_copy_if(vecStr[i], boost::algorithm::is_digit());
        //check whether the label is OK.
        if( CATAZJUT::CElement::isValidSymbol(ChemicalStr.first) == false ){
           Log::Error<<vecStr[i]<< " is error. Initialization_CCalcCluster!\n";
           boost::throw_exception(std::runtime_error(vecStr[i] + " is error. Initialization_CCalcCluster!\n"));//ERROR TREATMENT;
        }
        str = boost::algorithm::trim_copy_if(vecStr[i], boost::algorithm::is_alpha());
        if ( str=="" )
            ChemicalStr.second = 1;
        else
            ChemicalStr.second = std::stoi(str);
        this->chemicalFormula.push_back(ChemicalStr);
    }
    for(size_t i=0;i<this->m_pParameter->GaParameter()->PopNum();i++)
        RandomBuildFromChemicalFormula(this->m_PopuPeriodicFramework[i]);

}
void CCalcCluster::Initialization(const char* mth)
{      // initialize from chemical formula
     std::string tmp(mth);
     this->Initialization(tmp);
}
void CCalcCluster::Initialization(const std::vector<std::string*>& inputfiles)  // initialize from exit
{
   std::string str;
   std::vector<std::string> vecStr;
   Cios* inputIO;
   size_t pos=0;
   for(size_t i=0;i<inputfiles.size();i++){
       str=*(inputfiles[i]);
       boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
       if(boost::iequals(vecStr[1],"gjf"))
          inputIO = new CIOGjf(this->m_PopuPeriodicFramework[pos++]);
       else if(boost::iequals(vecStr[1],"car"))
          inputIO = new CIOCar(this->m_PopuPeriodicFramework[pos++]);
       else if(boost::iequals(vecStr[1],"mol"))
          inputIO = new CIOMol(this->m_PopuPeriodicFramework[pos++]);
       else if(boost::iequals(vecStr[1],"poscar")==0)
          inputIO = new CIOPoscar(this->m_PopuPeriodicFramework[pos++]);
       else{
          Log::Error<<vecStr[i] << " file is not supported by cluster model! Initialization_CCalcCluster\n";
          boost::throw_exception(std::runtime_error(vecStr[i]+ " file is not supported by cluster model! Initialization_CCalcCluster!"));
       }
       // read structure and save m_PopuPeriodicFramework[pos++]
       inputIO->input(str);
   }
   this->chemicalFormula= *(m_PopuPeriodicFramework[0]->composition());

   for(size_t i=pos;i<this->m_pParameter->GaParameter()->PopNum();i++)
      RandomBuildFromChemicalFormula(this->m_PopuPeriodicFramework[i]);
}
void CCalcCluster::RandomBuildFromChemicalFormula(CATAZJUT::CPeriodicFramework* predict_struct)
{
     std::vector<CATAZJUT::CElement*> res;
     size_t atom_Sum=0;
     for(size_t i=0;i<chemicalFormula.size();i++){
        res.push_back(new CATAZJUT::CElement(chemicalFormula[i].first));
        atom_Sum = atom_Sum + chemicalFormula[i].second;
     }
     size_t cluster_type = ClusterType(res);
     if(cluster_type==1){        //pure metal
        spherePredict(predict_struct);
     }else if(cluster_type==2){  //pure nonmetal

     }else if(cluster_type==3){   //pure mixed nonmetal

     }
}
size_t CCalcCluster::ClusterType(std::vector<CATAZJUT::CElement*>& mht)
{
   // size_t M_Pointer=0,NoM_Pointer=0, H_Pointer=0;
    util::Bitset Element_Bit(2);
     //set all bit to be 1;
    Element_Bit.set();
    for(size_t i=0;i<mht.size();i++)
    {
       if( mht[i]->isMetal() )    Element_Bit.set(0,false);
       if( mht[i]->isNonmetal() ) Element_Bit.set(1,false);
    }
    if(Element_Bit.to_ulong()== BOOST_BINARY(10))
        return 1;
    if(Element_Bit.to_ulong()== BOOST_BINARY(01))
        return 2;
    if(Element_Bit.to_ulong()== BOOST_BINARY(00))
        return 3;
    return 0;
}
// for metal clusters
void CCalcCluster::spherePredict(CATAZJUT::CPeriodicFramework* predict_struct)
{
    util::Vector3 polar_coord, basic_coord;
    util::Point3 coordinate;
    double cosTheta,sinTheta,cosPhi,sinPhi;

    double covalent_Radius=0.0;

    size_t atom_Sum=0;
    for(size_t i=0;i<this->chemicalFormula.size();i++){
        covalent_Radius = covalent_Radius + (new CATAZJUT::CElement(this->chemicalFormula[i].first))->covalentRadius();
        atom_Sum = atom_Sum + this->chemicalFormula[i].second;
    }
    // clear all atoms and all bonds;
    predict_struct->clear();
    covalent_Radius = covalent_Radius/atom_Sum;
    double bondrange = 2*covalent_Radius*(m_pParameter->bondToleranceFactor.first + \
                                          m_pParameter->bondToleranceFactor.second)*0.5;
    // determine:  radius, theta, Phi
    basic_coord<< bondrange*std::pow(atom_Sum,1.0/3),180.0,360.0;
    util::CRandomgenerator* myRandomgenerator = new util::CRandomgenerator();
    for(size_t i=0;i<chemicalFormula.size();i++)
        for(size_t j=0;j<chemicalFormula[i].second;j++){
                    // using polar coordinate to predict it
            polar_coord =basic_coord.cwiseProduct(myRandomgenerator->randomVector01(i*j));
                   // transfer polar coordinate to cartesian coordinate.
            cosTheta= std::cos(polar_coord(1)*CATAZJUT::constants::DegreesToRadians);
            sinTheta= std::sin(polar_coord(1)*CATAZJUT::constants::DegreesToRadians);
            cosPhi = std::cos(polar_coord(2)*CATAZJUT::constants::DegreesToRadians);
            sinPhi = std::sin(polar_coord(2)*CATAZJUT::constants::DegreesToRadians);
            coordinate<<polar_coord(0)*cosTheta*cosPhi,polar_coord(0)*cosTheta*sinPhi,polar_coord(0)*sinTheta;
            predict_struct->addAtom(this->chemicalFormula[i].first,coordinate);
        }
    //predict_struct->perceiveBonds(); // it is not necessary to perform it.
    this->eliminateFragment(predict_struct);
    this->eliminateCloseContacts(predict_struct);
}
// nonmetal compounds
void CCalcCluster::planePredict(CATAZJUT::CPeriodicFramework* predict_struct)
{
    std::vector<std::pair<CATAZJUT::CElement*,size_t>> chemicalelement;
    for(size_t i=0;i<chemicalFormula.size();i++)
      chemicalelement.push_back(std::make_pair(new CATAZJUT::CElement(chemicalFormula[i].first),\
                                               chemicalFormula[i].second));
    size_t PBlockAtomNum=0;
    for(size_t i=0;i<chemicalelement.size();i++)
        if(chemicalelement[i].first->maxCoordinationNum()>1)
            PBlockAtomNum = PBlockAtomNum + chemicalelement[i].second;
    Point3 position;
    for(size_t i=0;i<chemicalelement.size();i++)
        if(chemicalelement[i].first->maxCoordinationNum()>1)
            for(size_t j=0;j<chemicalelement[i].second;j++ ){
                position = 1.20*std::pow(PBlockAtomNum,0.5)*(Point3::Random().normalized());
                predict_struct->addAtom(chemicalelement[i].first->symbol(),position);
            }

   predict_struct->perceiveBonds();
   this->eliminateCloseContacts(predict_struct);
   this->eliminateFragment(predict_struct);
   double bondlength;
   bool isAdd;
   for(size_t i=0;i<chemicalelement.size();i++)
      if(chemicalelement[i].first->maxCoordinationNum()<=1){
         for(size_t j=0;j<chemicalelement[i].second;j++){
            isAdd=false;
            foreach(CATAZJUT::CAtom* atom,predict_struct->atoms())
                if(atom->isBondSaturated() == false){
                    bondlength = chemicalelement[i].first->covalentRadius() + atom->element().covalentRadius();
                    predict_struct->addAtom(chemicalelement[i].first->symbol(),atom->position()+bondlength*atom->NewBondingVect());
                    isAdd=true;
                    break;
                }
            if(!isAdd){
                position = 1.20*std::pow(PBlockAtomNum,0.5)*(Point3::Random().normalized());
                predict_struct->addAtom(chemicalelement[i].first->symbol(),position);
              }
          }
      }
   predict_struct->perceiveBonds();
   this->eliminateCloseContacts(predict_struct);
   this->eliminateFragment(predict_struct);

   // clear heap space
   for(size_t i=0;i<chemicalelement.size();i++)
       delete chemicalelement[i].first;
   chemicalelement.erase();

}
void CCalcCluster::eliminateCloseContacts(CATAZJUT::CPeriodicFramework* curr_struct,double distanceCutOff)
{
    util::Vector3 vect;
    util::Point3  center_P;
    double eps=0.1;
    bool modifiedbol=true;
    while(modifiedbol)
    {
       modifiedbol=false;
       // get center pointer of structure
       center_P = curr_struct->center();

       foreach(CATAZJUT::CAtom* atom_s, curr_struct->atoms())
          foreach(CATAZJUT::CAtom* atom_m, curr_struct->atoms())
            if( curr_struct->distance(atom_m,atom_s) < distanceCutOff ){
                if( (center_P - atom_m->position()).norm() > (center_P - atom_s->position()).norm() ){
                    vect = atom_m->position() - atom_s->position();
                    vect = ( distanceCutOff - vect.norm()+eps )*(vect.normalized());
                    curr_struct->moveAtom(atom_m,vect);
                }else{
                    vect = atom_s->position() - atom_m->position();
                    vect = ( distanceCutOff - vect.norm()+eps )*(vect.normalized());
                    curr_struct->moveAtom(atom_s,vect);
                }
                modifiedbol=true;
            }
    }
}
void CCalcCluster::eliminateFragment(CATAZJUT::CPeriodicFramework* curr_struct)
{
    if(! curr_struct->fragmentsPerceived())    // analysize the fragments of the whole structure
         curr_struct->perceiveFragments();
    if( curr_struct->fragmentNum() > 1 ){
        size_t maxCount=0;
        CATAZJUT::CFragment* mainFragment;
        foreach(CATAZJUT::CFragment* fragment_s, curr_struct->fragments()){
            if(maxCount<fragment_s->atomCount()){
                maxCount=fragment_s->atomCount();
                mainFragment=fragment_s;
            }
        }
        //
        Matrix mainsphereEquation4(4,1), othersphereEquation4(4,1);
        Point3 maincenter,othercenter;
        double mainR, otherR, differ;
        Vector3 differVect;
        //main center, radius of the largest fragment.
        //all other fragments move toward it.
        mainsphereEquation4 = util::SphereEquationFromPoints(mainFragment->coordinates());
        maincenter<<mainsphereEquation4(0,0),mainsphereEquation4(1,0),mainsphereEquation4(2,0);
        mainR=mainsphereEquation4(3,0);
        foreach(CATAZJUT::CFragment* fragment_s, curr_struct->fragments())
            if(mainFragment!=fragment_s){
               othersphereEquation4 = util::SphereEquationFromPoints(fragment_s->coordinates());
               othercenter<<othersphereEquation4(0,0),othersphereEquation4(1,0),othersphereEquation4(2,0);
               otherR=othersphereEquation4(3,0);
               differ = (maincenter-othercenter).norm() - mainR - otherR;
               differVect= (differ-1.5)*(maincenter-othercenter).normalized();
               fragment_s->move(differVect);
            }
    }

}



}
