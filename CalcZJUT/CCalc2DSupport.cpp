#include <Eigen/Core>
#include <Eigen/Dense>
#include <float.h>
#include "CCalc2DSupport.h"
#include "CPeriodicFramework.h"
#include "CCalcMoleculeAdsorbent.h"
#include "CPeriodicFramework.h"
#include "GaDeclaration.h"
#include "CCartesianCoordinates.h"
#include "CFractionCoordinates.h"
#include "CAtom.h"
#include "CUnitCell.h"
#include "CCalcSupportBase.h"
#include "foreach.h"
#include "CPlane.h"
#include "CBondTolerance.h"
#include "GaUtilityFunction.h"

using GAZJUT::ERROR_OUTPUT;
namespace CALCZJUT{


CCalc2DSupport::CCalc2DSupport(CParameter* mPara)
:CCalcModeStruct(mPara)
{
    m_pPeriodicFramework = new CATAZJUT::CPeriodicFramework(mPara);
    m_backupPeriodicFramework = m_pPeriodicFramework;
    m_latticeDirection = CCalc2DSupport::NONE_DIR;
}

CCalc2DSupport::~CCalc2DSupport()
{
     if( m_pSupport !=nullptr) delete m_pSupport;
     if( m_pAdsorbMolecule !=nullptr) delete m_pAdsorbMolecule;
     if( m_support_surface !=nullptr) delete m_support_surface;
    //dtor
}

//attribute operation
void CCalc2DSupport::createSupport(Bitset& mht)
{
    assert(m_pPeriodicFramework);
    m_BitbackupSupport = mht;
    m_pSupport= new CCalcSupportBase(this->m_pPeriodicFramework,mht);
}
void CCalc2DSupport::createMoleAdsorb(Bitset& mht)
{
    assert(m_pPeriodicFramework);
    m_BitbackupAdsorbMolecule = mht;
    m_pAdsorbMolecule = new CCalcMoleculeAdsorbent(this->m_pPeriodicFramework,mht);
}
void CCalc2DSupport::setPeriodicFramekwork(CATAZJUT::CPeriodicFramework* mbf)
{
    this->m_pPeriodicFramework=mbf;
    this->m_pAdsorbMolecule->setConfiguration(mbf);
    this->m_pSupport->setConfiguration(mbf);
}
//void CCalc2DSupport::backUpStructure()
//{
//   assert(m_pPeriodicFramework);
//   assert(m_pSupport);
//   assert(m_pAdsorbMolecule);
//
//   m_backupPeriodicFramework = new (CATAZJUT::CPeriodicFramework)(*m_pPeriodicFramework);
//   m_BitbackupSupport = m_pSupport->bitSet();
//   m_BitbackupAdsorbMolecule = m_pAdsorbMolecule->bitSet();
//}
//void CCalc2DSupport::fromBackupToCurrent()
//{
//   if(m_pPeriodicFramework->coordinateType()==CATAZJUT::DEFINED::Cartesian){
//       m_pPeriodicFramework->coordinates()->clear();
//       foreach(const CATAZJUT::CAtom* atom,m_backupPeriodicFramework->atoms())
//            m_pPeriodicFramework->coordinates()->setPosition(atom->index(),atom->position());
//
//   }else if(m_pPeriodicFramework->coordinateType()==CATAZJUT::DEFINED::Fraction){
//       m_pPeriodicFramework->Fractioncoordinates()->clear();
//       foreach(const CATAZJUT::CAtom* atom,m_backupPeriodicFramework->atoms())
//           m_pPeriodicFramework->Fractioncoordinates()->setPosition(atom->index(),atom->position());
//   }
//}
void CCalc2DSupport::createStructureAtGene()
{
   CCalcModeStruct::createStructureAtGene();
   m_pSupport->setConfiguration(this->periodicFramework());
   m_pAdsorbMolecule->setConfiguration(this->periodicFramework());
}
void CCalc2DSupport::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{
    CATAZJUT::Point3 vect;
    CATAZJUT::Point3 target_p;
    double a,b;
    //identify the surface on the support
    if( m_support_surface ==nullptr)
        this->perceiveSupportSurface()

    if(m_latticeDirection == CCalc2DSupport::C_AXIS){

        a=realValueOfgene[1];     //in the coordinate system of lattice
        b=realValueOfgene[2];
        target_p<<a,b,0;
        //transform to
        target_p= m_pPeriodicFramework->unitcell()->NormilizedBravaisMatrix()*target_p;
        target_p(2)=target_p(2) + m_support_surface->Distance(target_p) + realValueOfgene[0];

        target_p=target_p - m_pAdsorbMolecule->gravityCentre();
        m_pAdsorbMolecule->moveBy( target_p );
        //rotate along axis and angle
        target_p<<realValueOfgene[3],realValueOfgene[4],realValueOfgene[5];
        m_pAdsorbMolecule->rotate(target_p,realValueOfgene[6]);
    }else if(m_latticeDirection == CCalc2DSupport::B_AXIS{
        a=realValueOfgene[1];     //in the coordinate system of lattice
        b=realValueOfgene[2];
        target_p<<a,0,b;
        //transform to
        target_p= m_pPeriodicFramework->unitcell()->NormilizedBravaisMatrix()*target_p;
        target_p(1)=target_p(1) + m_support_surface->Distance(target_p) + realValueOfgene[0];

        target_p=target_p - m_pAdsorbMolecule->gravityCentre();
        m_pAdsorbMolecule->moveBy( target_p );
        //rotate along axis and angle
        target_p<<realValueOfgene[3],realValueOfgene[4],realValueOfgene[5];
        m_pAdsorbMolecule->rotate(target_p,realValueOfgene[6]);
    }else if(m_latticeDirection == CCalc2DSupport::A_AXIS{
        a=realValueOfgene[1];     //in the coordinate system of lattice
        b=realValueOfgene[2];
        target_p<<0,a,b;
        //transform to
        target_p= m_pPeriodicFramework->unitcell()->NormilizedBravaisMatrix()*target_p;
        target_p(0)=target_p(0) + m_support_surface->Distance(target_p) + realValueOfgene[0];

        target_p=target_p - m_pAdsorbMolecule->gravityCentre();
        m_pAdsorbMolecule->moveBy( target_p );
        //rotate along axis and angle
        target_p<<realValueOfgene[3],realValueOfgene[4],realValueOfgene[5];
        m_pAdsorbMolecule->rotate(target_p,realValueOfgene[6]);

    }else{
       ERROR_OUTPUT("Lattice Direction is error!","setGeneValueToStruct","CCalc2DSupport");
       boost::throw_exception(std::runtime_error("Lattice Direction is error!! Check the file: Error_information.txt."));
    }

     // eliminate the closecontact in new structure
    this->eliminateCloseContacts();
}
void CCalc2DSupport::eliminateCloseContacts(double distanceCutOff)
{
    CATAZJUT::Vector3 vect;
    double eps=0.01;
    bool modifiedbol=true;
    while(modifiedbol)
    {
       modifiedbol=false;
       foreach(CATAZJUT::CAtom* atom_s, m_pSupport->atoms()){
          foreach(CATAZJUT::CAtom* atom_m, m_pAdsorbMolecule->atoms()){
            if( m_pPeriodicFramework->distance(atom_m,atom_s) <distanceCutOff ){
                vect = atom_m->position() - atom_s->position();
                vect = (distanceCutOff-vect.norm()+eps)*(vect.normalized());
                this->m_pAdsorbMolecule->moveBy(vect);
                modifiedbol=true;
            }
          }
       }
    }
}
std::vector<double>* CCalc2DSupport::getGeneValuefromStruct()const
{
   return nullptr;
}
std::vector<GENEVAR>* CCalc2DSupport::GeneVarRange()
{
    std::vector<GENEVAR>* res= new (std::vector<GENEVAR>);
    //distance of adsorbent on 2d support
    res->push_back({1.0,2.5,0.01});
    // a axes displacing
    if(m_latticeDirection==CCalc2DSupport::NONE_DIR)
        IdentifyvacuumLayerDirection();

    if(m_latticeDirection == CCalc2DSupport::C_AXIS){
        //c a,b
        res->push_back({0,m_pPeriodicFramework->unitcell()->a(),0.01});
        res->push_back({0,m_pPeriodicFramework->unitcell()->b(),0.01});
    }else if (m_latticeDirection == CCalc2DSupport::B_AXIS){
        //b c  a
        res->push_back({0,m_pPeriodicFramework->unitcell()->c(),0.01});
        res->push_back({0,m_pPeriodicFramework->unitcell()->a(),0.01});
        m_latticeDirection = CCalc2DSupport::B_AXIS;
    }else if(m_latticeDirection == CCalc2DSupport::A_AXIS){
        //a b  c
        res->push_back({0,m_pPeriodicFramework->unitcell()->b(),0.01});
        res->push_back({0,m_pPeriodicFramework->unitcell()->c(),0.01});
    }else{
         ERROR_OUTPUT("It has no vacuum layer!","GeneVarRange","CCalc2DSupport");
         boost::throw_exception(std::runtime_error("It has no vacuum layer!! Check the file: Error_information.txt."));
    }
    //rotation axis
    res->push_back({-1,1,0.01});
    res->push_back({-1,1,0.01});
    res->push_back({-1,1,0.01});

    //rotation angle
    res->push_back({0.0,360.0,2.0});

    return res;

}
void CCalc2DSupport::IdentifyvacuumLayerDirection()  //-1: none, 0:a,  1:b  2:c
{

    Eigen::Matrix<double, 3, 3> latticeVect(m_pPeriodicFramework->unitcell()->MatrixOfBravaisLattice());

    double dis;
    bool breakBol=false;
    CATAZJUT::Vector3 disVect;
    for(size_t i=0;i<3;i++){
       breakBol=false;
       foreach(const CATAZJUT::CAtom* atom_1,m_pSupport->atoms()){
         foreach(const CATAZJUT::CAtom* atom_2,m_pSupport->atoms()){
             disVect=latticeVect.col(i);
             disVect= latticeVect*(disVect.normalized());
              dis= (atom_1->position() - (atom_2->position()+disVect)).norm();
              if(m_pPeriodicFramework->m_pBondEvaluator->IsBond(atom_1->Symbol(), \
                 atom_2->Symbol(),dis)==true){
                    breakBol=true;
                    break;
              }
         }
         if(breakBol)
            break;
       }
       if(!breakBol){
          switch (i){
            case 0:
                m_latticeDirection = CCalc2DSupport::A_AXIS;
                break;
            case 1:
                m_latticeDirection = CCalc2DSupport::B_AXIS;
                break;
            case 2:
                m_latticeDirection = CCalc2DSupport::C_AXIS;
                break;
            default:
                m_latticeDirection = CCalc2DSupport::NONE_DIR;
                break;
          }
          break;
       }
    }
}
void CCalc2DSupport::perceiveSupportSurface()
{
    std::vector<size_t> coord_Num;
    int direction;
    if ( m_latticeDirection == CCalc2DSupport::C_AXIS)
        direction=2;
    else if(m_latticeDirection == CCalc2DSupport::B_AXIS)
        direction=1;
    else //if (m_latticeDirection = CCalc2DSupport::A_AXIS)
        direction=0;

    double max_latticedirection=DBL_MIN;
    foreach(const CATAZJUT::CAtom *atom,m_pSupport->atoms()){
        if(atom->position()[direction]>max_latticedirection)
            max_latticedirection=atom->position()[direction];
    }
    std::vector<CATAZJUT::CAtom*> matom;
    double PlaneCutoff = 0.1;
    while(true)
    {
        matom.clear();
        foreach(CATAZJUT::CAtom *atom,m_pSupport->atoms()){
           if(std::fabs(atom->position()[direction]-max_latticedirection)<PlaneCutoff)
              matom.push_back(atom);
        }
        if(matom.size()>4)
            break;
        else if(PlaneCutoff >1.0){
            ERROR_OUTPUT(" 2D support has serious error!","perceiveSupportSurface","CCalc2DSupport");
            boost::throw_exception(std::runtime_error("2D support has serious error! Check the file: Error_information.txt."));
        }else
            PlaneCutoff=PlaneCutoff + 0.1;
    }
    Eigen::MatrixXd* currentMat = new (Eigen::MatrixXd)(matom.size(),3);
    for(size_t i=0;i<matom.size();i++)
        (*currentMat).row(i) = matom[i]->position();

    this->m_support_surface = new CATAZJUT::CPlane(currentMat);
    delete currentMat;
}


}
