/******************************************************************************
**
** Copyright (C) 2019-2031 Dr.Gui-lin Zhuang <glzhuang@zjut.edu.cn>
** All rights reserved.
**
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <strings.h>
#include <cstring>
#include <math.h>
#include "../CataZJUT/CPeriodicFramework.h"
#include "CModelCluster.h"
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
#include "CIOBase.h"
#include "CIOMol.h"
#include "CIOCar.h"
#include "CIOGjf.h"
#include "CIOPoscar.h"
#include "CStructPoolCluster.h"
#include "../IACS.h"

using util::Log;
using util::Bitset;
using util::Vector4;

namespace CALCZJUT{

CStructPoolCluster::CStructPoolCluster(CParameter* othr)
:CStructPoolBase(othr)
{
    #ifdef DEBUG
      Log::Debug<<"*********** CStructPoolCluster::CStructPoolCluster***********"<< std::endl;
    #endif
    for(size_t i=0;i<this->m_pParameter->GaParameter()->PopNum();i++){
        this->m_CalcStructPool.push_back(new CModelCluster(this->m_pParameter,i));
       // m_CalcStructPool[m_CalcStructPool.size()-1]->periodicFramework()->setExcludeBond(m_pParameter->excludeBond);
       // m_CalcStructPool[m_CalcStructPool.size()-1]->periodicFramework()->setTolerancefactor(m_pParameter->bondToleranceFactor);
    }
    #ifdef DEBUG
      Log::Debug<<"End*********** CStructPoolCluster::CStructPoolCluster***********"<< std::endl;
    #endif
}

CStructPoolCluster::~CStructPoolCluster()
{

}

void CStructPoolCluster::init()
{
    if(this->m_pParameter->cluster_Input_File.size()!=0){
        this->Initialization(m_pParameter->cluster_Input_File);
    }else if(this->m_pParameter->cluster_Formula!=""){
        this->Initialization(m_pParameter->cluster_Formula);
    }else{
       Log::Error<< " Chemical formula or initially structural files is required. init_CModelCluster!\n";
       boost::throw_exception(std::runtime_error("Chemical formula and structural files is required. init_CModelCluster!!\n"));//ERROR TREATMENT;
    }
}
void CStructPoolCluster::GeneVARRange(std::vector<GeneVAR>&  mht)
{
   double maxRadius=0.0;
   std::vector<GeneVAR> resTemp;
   for(size_t i=0;i<this->m_CalcStructPool.size();i++){
        this->m_CalcStructPool[i]->GeneVARRange(resTemp);
        if(maxRadius > resTemp[0].max)
            maxRadius = resTemp[0].max;
   }
   maxRadius = maxRadius + 1.0;
   size_t num=this->m_CalcStructPool[0]->m_pPeriodicFramework->size();

   for(size_t i=0;i<3*num;i++)
      mht.push_back({-1*maxRadius,maxRadius,0.001});
}
void CStructPoolCluster::Initialization(const std::string& mth)
{
       // initialize from chemical
       // chemical formula = C10 H10 O12
/*
    e.g.:   mth = C10 N S10        // delimiter = blank " "
*/
    //dealwith chemical formula by blank
    #ifdef DEBUG
      Log::Debug<<"*********** CStructPoolCluster::Initialization***********"<< std::endl;
    #endif
    std::string str=mth;
    std::vector<std::string> vecStr;
    boost::algorithm::trim(str);
    boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);

    std::pair<std::string,size_t> ChemicalStr;
    std::vector<std::pair<std::string,size_t>> tempChemFormula;

    for(size_t i=0;i<vecStr.size();i++){
        ChemicalStr.first  = boost::algorithm::trim_copy_if(vecStr[i], boost::algorithm::is_digit());
        //check whether the label is OK.
        if( CATAZJUT::CElement::isValidSymbol(ChemicalStr.first) == false ){
           Log::Error<<vecStr[i]<< " is error. Initialization_CModelCluster!\n";
           boost::throw_exception(std::runtime_error(vecStr[i] + " is error. Initialization_CModelCluster!\n"));//ERROR TREATMENT;
        }
        str = boost::algorithm::trim_copy_if(vecStr[i], boost::algorithm::is_alpha());
        if ( str=="" )
            ChemicalStr.second = 1;
        else
            ChemicalStr.second = std::stoi(str);
        tempChemFormula.push_back(ChemicalStr);
    }
    #ifdef DEBUG
      Log::Debug<<"*********** 2- CStructPoolCluster::Initialization***********"<< std::endl;
    #endif
    for(size_t i=0;i<this->m_pParameter->GaParameter()->PopNum();i++)
        RandomBuildFromChemicalFormula(this->m_CalcStructPool[i]->periodicFramework(),tempChemFormula);
    tempChemFormula.clear();

    #ifdef DEBUG
      Log::Debug<<"*********** End CStructPoolCluster::Initialization***********"<< std::endl;
    #endif
}
void CStructPoolCluster::Initialization(const char* mth)
{      // initialize from chemical formula
     std::string tmp(mth);
     this->Initialization(tmp);
}
void CStructPoolCluster::Initialization(const std::vector<std::string*>& inputfiles)  // initialize from exit
{
   std::string str;
   std::vector<std::string> vecStr;
   CIOBase* inputIO;
   size_t pos=0;
   for(size_t i=0;i<inputfiles.size();i++){
       str=*(inputfiles[i]);
       boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
       if(boost::iequals(vecStr[1],"gjf"))
          inputIO = new CIOGjf(this->m_CalcStructPool[pos]->periodicFramework());
       else if(boost::iequals(vecStr[1],"car"))
          inputIO = new CIOCar(this->m_CalcStructPool[pos]->periodicFramework());
       else if(boost::iequals(vecStr[1],"mol"))
          inputIO = new CIOMol(this->m_CalcStructPool[pos]->periodicFramework());
       else if(boost::iequals(vecStr[1],"poscar")==0)
          inputIO = new CIOPoscar(this->m_CalcStructPool[pos]->periodicFramework());
       else{
          Log::Error<<vecStr[i] << " file is not supported by cluster model! Initialization_CModelCluster\n";
          boost::throw_exception(std::runtime_error(vecStr[i]+ " file is not supported by cluster model! Initialization_CModelCluster!"));
       }
       // read structure and save m_PopuPeriodicFramework[pos++]
       inputIO->input(str);

       m_CalcStructPool[pos]->setChemicalFormula( m_CalcStructPool[pos]->periodicFramework()->composition() );
       pos++;
   }

   for(size_t i=pos;i<this->m_pParameter->GaParameter()->PopNum();i++){
      this->m_CalcStructPool[i]->setChemicalFormula(this->m_CalcStructPool[0]->chemicalFormula());
      RandomBuildFromChemicalFormula(this->m_CalcStructPool[i]->periodicFramework(),
                                     this->m_CalcStructPool[i]->chemicalFormula());
      }
}
void CStructPoolCluster::RandomBuildFromChemicalFormula(CATAZJUT::CPeriodicFramework* predict_struct,
                                                            std::vector<std::pair<std::string,size_t>>& chemFormula)
{
     std::vector<CATAZJUT::CElement*> res;
     size_t atom_Sum=0;
     for(size_t i=0;i<chemFormula.size();i++){
        res.push_back(new CATAZJUT::CElement(chemFormula[i].first));
        atom_Sum = atom_Sum + chemFormula[i].second;
     }
     #ifdef DEBUG
      Log::Debug<<"*********** CStructPoolCluster::RandomBuildFromChemicalFormula***********"<< std::endl;
     #endif
     size_t cluster_type = ClusterType(res);
     if(cluster_type==1){        //pure metal
        metalClusterPredict(predict_struct,chemFormula);
     }else if(cluster_type==2){  //pure nonmetal
        nonMetalClusterPredict(predict_struct,chemFormula);
     }else if(cluster_type==3){   //pure mixed nonmetal
        mixedClusterPredict(predict_struct,chemFormula);
     }
}
size_t CStructPoolCluster::ClusterType(std::vector<CATAZJUT::CElement*>& mht)
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
void CStructPoolCluster::metalClusterPredict(CATAZJUT::CPeriodicFramework* predict_struct,
                                                 std::vector<std::pair<std::string,size_t>>& chemFormula)
{
    util::Vector3 polar_coord, basic_coord;
    util::Point3 coordinate;
    double cosTheta,sinTheta,cosPhi,sinPhi;

    double covalent_Radius=0.0;

    size_t atom_Sum=0;
    for(size_t i=0;i<chemFormula.size();i++){
        covalent_Radius = covalent_Radius + (new CATAZJUT::CElement(chemFormula[i].first))->covalentRadius();
        atom_Sum = atom_Sum + chemFormula[i].second;
    }
    // clear all atoms and all bonds;
    predict_struct->clear();
    covalent_Radius = covalent_Radius/chemFormula.size();
    double bondrange = 2.0*covalent_Radius*(m_pParameter->bondToleranceFactor.first + \
                                          m_pParameter->bondToleranceFactor.second)*0.5;
    // determine:  radius, theta, Phi
//    #ifdef DEBUG
//        std::cout<<"covalent_Radius:"<<covalent_Radius<<std::endl;
//        std::cout<<"bondrange:"<<bondrange<<std::endl;
//     #endif
    basic_coord<< bondrange*std::pow(atom_Sum,1.0/3),180.0,360.0;
    util::CRandomgenerator* myRandomgenerator = new util::CRandomgenerator();
    for(size_t i=0;i<chemFormula.size();i++)
        for(size_t j=0;j<chemFormula[i].second;j++){
                    // using polar coordinate to predict it
            polar_coord =basic_coord.cwiseProduct(myRandomgenerator->randomVector01((i+2)*(j+3)));
                   // transfer polar coordinate to cartesian coordinate.
            cosTheta= std::cos(polar_coord(1)*CATAZJUT::constants::DegreesToRadians);
            sinTheta= std::sin(polar_coord(1)*CATAZJUT::constants::DegreesToRadians);
            cosPhi = std::cos(polar_coord(2)*CATAZJUT::constants::DegreesToRadians);
            sinPhi = std::sin(polar_coord(2)*CATAZJUT::constants::DegreesToRadians);
            coordinate<<polar_coord(0)*cosTheta*cosPhi,polar_coord(0)*cosTheta*sinPhi,polar_coord(0)*sinTheta;
//            #ifdef DEBUG
//               std::cout<<chemFormula[i].first<<j+1<<": "<<coordinate.transpose()<<std::endl;
//            #endif
            predict_struct->addAtom(chemFormula[i].first,coordinate);
        }
     predict_struct->perceiveBonds(); // it is not necessary to perform it.
//     #ifdef DEBUG
//        Log::Debug<<"*********** output metalClusterPredict***********"<< std::endl;
//        this->outPutStructure(predict_struct);
//     #endif
    this->eliminateFragment(predict_struct);
    this->eliminateCloseContacts(predict_struct);
}
// nonmetal compounds
void CStructPoolCluster::nonMetalClusterPredict(CATAZJUT::CPeriodicFramework* predict_struct,
                                                    std::vector<std::pair<std::string,size_t>>& chemFormula)
{
    std::vector<std::pair<CATAZJUT::CElement*,size_t>> chemicalelement;
    for(size_t i=0;i<chemFormula.size();i++)
      chemicalelement.push_back(std::make_pair(new CATAZJUT::CElement(chemFormula[i].first),\
                                               chemFormula[i].second));
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
   chemicalelement.clear();

}
// nonmetal-metal cluster
void CStructPoolCluster::mixedClusterPredict(CATAZJUT::CPeriodicFramework* predict_struct,
                                             std::vector<std::pair<std::string,size_t>>& chemFormula)
{
    //GACatalysis(nousing)
}
void CStructPoolCluster::eliminateCloseContacts(CATAZJUT::CPeriodicFramework* curr_struct,double distanceCutOff)
{
    util::Vector3 vect;
    util::Point3  center_P;
    double eps=0.1;
    bool modifiedbol=true;
    size_t less_cycle=0;

    assert(curr_struct);

//    #ifdef DEBUG
//      Log::Debug<<"*********** eliminateCloseContacts***********"<< std::endl;
//    #endif
    while(modifiedbol)
    {
       modifiedbol=false;
       // get center pointer of structure
       center_P = curr_struct->center();

       foreach(CATAZJUT::CAtom* atom_s, curr_struct->atoms())
          foreach(CATAZJUT::CAtom* atom_m, curr_struct->atoms()){
            distanceCutOff = (atom_s->CovalentRadius() + atom_m->CovalentRadius())*0.6;
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
       less_cycle++;
       if(less_cycle > 256){
          Log::Warn<<"Cycle of function is more than 256 in CStructPoolCluster::eliminateCloseContacts!"<<std::endl;
          break;
       }
    }
}
void CStructPoolCluster::eliminateFragment(CATAZJUT::CPeriodicFramework* curr_struct)
{
    assert(curr_struct);

    size_t less_cycle;
    if(! curr_struct->fragmentsPerceived())    // analysize the fragments of the whole structure
         curr_struct->perceiveFragments();
    if( curr_struct->fragmentNum() > 1 ){
        size_t maxCount=0;
        CATAZJUT::CFragment* mainFragment;
        foreach(CATAZJUT::CFragment* fragment_s, curr_struct->fragments()){
            if(maxCount < fragment_s->atomCount()){
                maxCount = fragment_s->atomCount();
                mainFragment = fragment_s;
            }
            #ifdef DEBUG
                  Log::Debug<<"*2- CStructPoolCluster::eliminateFragment:"<<maxCount<< std::endl;
                  Log::Debug<<"*2- CStructPoolCluster::eliminateFragment:fragmentNum"<<curr_struct->fragmentNum()<< std::endl;
            #endif
        }
        //
        #ifdef DEBUG
           Log::Debug<<"***********3- CStructPoolCluster::eliminateFragment***********"<< std::endl;
        #endif
        Vector4 mainsphereEquation4, othersphereEquation4;
        Point3 maincenter,othercenter;
        double mainR, otherR, differ;
        Vector3 differVect;
        //main center, radius of the largest fragment.
        //all other fragments move toward it.
        mainsphereEquation4 = util::SphereEquationFromPoints(mainFragment->coordinates());
        maincenter<<mainsphereEquation4(0,0),mainsphereEquation4(1,0),mainsphereEquation4(2,0);
        mainR=mainsphereEquation4(3,0);
        #ifdef DEBUG
           Log::Debug<<"***********4- CStructPoolCluster::eliminateFragment***********"<< std::endl;
        #endif
        foreach(CATAZJUT::CFragment* fragment_s, curr_struct->fragments())
            if(mainFragment!=fragment_s){
               if(fragment_s->atomCount()>1){
                   othersphereEquation4 = util::SphereEquationFromPoints(fragment_s->coordinates());
                   othercenter<<othersphereEquation4(0,0),othersphereEquation4(1,0),othersphereEquation4(2,0);
                   otherR=othersphereEquation4(3,0);
                   differ = (maincenter-othercenter).norm() - mainR - otherR;
                   differVect= (differ-1.5)*(maincenter-othercenter).normalized();
                   fragment_s->move(differVect);
               }else{
                   othercenter=fragment_s->coordinates()[0];
                   differ = (maincenter-othercenter).norm() - mainR;
                   differVect= (differ-1.5)*(maincenter-othercenter).normalized();
                   fragment_s->move(differVect);
               }
               // further judge whether two fragments is bonded.
               less_cycle=0;
               #ifdef DEBUG
                  Log::Debug<<"***********3- CStructPoolCluster::eliminateFragment***********"<< std::endl;
               #endif
               while(mainFragment->isBondTo(fragment_s)!=true){
                   differVect= 0.2*(maincenter-othercenter).normalized();
                   fragment_s->move(differVect);
                   less_cycle++;
//                   #ifdef DEBUG
//                     std::cout<<less_cycle<<std::endl;
//                   #endif // DEBUG
                   if(less_cycle>256){
                      Log::Warn<<"Cycle of function is more than 65536 in CStructPoolCluster::eliminateFragment!"<<std::endl;
                      break;
                   }
               }
            }
    }
}
void CStructPoolCluster::outPutStructure(CATAZJUT::CPeriodicFramework* structure)
{
    assert(structure);
    CIOMol* outmol =new CIOMol(structure);
    outmol->output("Pd10_temp.mol");
    delete outmol;
}




}
