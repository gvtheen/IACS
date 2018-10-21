#include <Eigen/Core>
#include <Eigen/Dense>
#include <float.h>
#include "CModel2DSupport.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "CModelMoleculeAdsorbent.h"
#include "../GaZJUT/GaDeclaration.h"
#include "../CataZJUT/CCartesianCoordinates.h"
#include "../CataZJUT/CFractionCoordinates.h"
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CUnitCell.h"
#include "CModelSupport.h"
#include "../Util/foreach.h"
#include "../CataZJUT/CPlane.h"
#include "../CataZJUT/CBondTolerance.h"
#include "../GaZJUT/GaUtilityFunction.h"
#include  "../Util/Point-Vector.h"

using util::Log;
using util::Point3;
using util::Vector3;

namespace CALCZJUT{


CModel2DSupport::CModel2DSupport(CParameter* mPara,
                                 CATAZJUT::CPeriodicFramework** copy_ppPeriodicFramework,
                                 size_t index)
:CModelBase(mPara,index)
{
    m_pPeriodicFramework = new CATAZJUT::CPeriodicFramework(mPara);
    m_ppBackupPeriodicFramework = copy_ppPeriodicFramework;
    m_latticeDirection = CModel2DSupport::NONE_DIR;
}

CModel2DSupport::~CModel2DSupport()
{
     if( m_pSupport !=nullptr) delete m_pSupport;
     if( m_pAdsorbMolecule !=nullptr) delete m_pAdsorbMolecule;
     if( m_support_surface !=nullptr) delete m_support_surface;
    //dtor
}
CModelBase* CModel2DSupport::clone()
{
    CModel2DSupport* res= new CModel2DSupport(this->m_pParameter,         \
                                              this->m_ppBackupPeriodicFramework,    \
                                              0);
    res->setPeriodicFramekwork(this->periodicFramework());
    res->createMoleAdsorb(this->MoleAdsorbBit());
    res->createSupport(this->SupportBit());
    res->setLatticeDirection(this->latticeDirection());
    res->setSupportSurfacePlane(this->supportSurfacePlane());
    return res;

}
//attribute operation
void CModel2DSupport::createSupport(const Bitset& mht)
{
    assert(m_pPeriodicFramework);
    m_BitbackupSupport = mht;
    m_pSupport= new CModelSupport(this->m_pPeriodicFramework,m_BitbackupSupport);
}
void CModel2DSupport::createMoleAdsorb( const Bitset& mht)
{
    assert(m_pPeriodicFramework);
    m_BitbackupAdsorbMolecule = mht;
    m_pAdsorbMolecule = new CModelMoleculeAdsorbent(this->m_pPeriodicFramework,m_BitbackupAdsorbMolecule);
}
Bitset CModel2DSupport::SupportBit()
{
    return this->m_BitbackupSupport;
}
Bitset CModel2DSupport::MoleAdsorbBit()
{
   return this->m_BitbackupAdsorbMolecule;
}
void CModel2DSupport::setPeriodicFramekwork(CATAZJUT::CPeriodicFramework* mbf)
{
    if(this->m_pPeriodicFramework!=nullptr)
        delete this->m_pPeriodicFramework;
    this->m_pPeriodicFramework=new CATAZJUT::CPeriodicFramework(*mbf);
    this->m_pAdsorbMolecule->setConfiguration(this->m_pPeriodicFramework);
    this->m_pSupport->setConfiguration(this->m_pPeriodicFramework);
}

CATAZJUT::CPlane* CModel2DSupport::supportSurfacePlane()
{
    return this->m_support_surface;
}
void CModel2DSupport::setSupportSurfacePlane(CATAZJUT::CPlane* mth)
{
   if(this->m_support_surface!=nullptr)
       delete this->m_support_surface;
   this->m_support_surface = new CATAZJUT::CPlane(*mth);
}

CModel2DSupport::LATT_DIRECTION CModel2DSupport::latticeDirection()
{
   return this->m_latticeDirection;
}
void CModel2DSupport::setLatticeDirection(CModel2DSupport::LATT_DIRECTION mth)
{
   this->m_latticeDirection=mth;
}

void CModel2DSupport::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{
    Point3 vect,target_p;
    double a,b;
    //identify the surface on the support
    if(this->RandomInitState()==false){
        this->setRandomInitState(true);
        goto RETURN_Random_Label;
    }

    /*
        check backup configuration;
        import it to the currentPointer of configuration;
    */
    if( *m_ppBackupPeriodicFramework == nullptr ){
       Log::Error<<"Backup pointer of the configuration is null!  setGeneValueToStruct_CModel2DSupport";
       boost::throw_exception(std::runtime_error("Backup pointer of the configuration is null.setGeneValueToStruct_CCalc2DSuppor"));
    }else{
       delete this->m_pPeriodicFramework;
       this->m_pPeriodicFramework = (*m_ppBackupPeriodicFramework)->clone();
    }
    if( m_support_surface ==nullptr)
        this->perceiveSupportSurface();

    if(m_latticeDirection == CModel2DSupport::C_AXIS){

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
    }else if(m_latticeDirection == CModel2DSupport::B_AXIS){
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
    }else if(m_latticeDirection == CModel2DSupport::A_AXIS){
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
       Log::Error<<"Lattice Direction is error! setGeneValueToStruct_CModel2DSupport";
       boost::throw_exception(std::runtime_error("Lattice Direction is error!! Check the file: Error_information.txt."));
    }

     // eliminate the closecontact in new structure
    this->eliminateCloseContacts();

RETURN_Random_Label:
    ;
    //nothing;
}
void CModel2DSupport::eliminateCloseContacts(double distanceCutOff)
{
    Vector3 vect;
    double eps=0.01;
    bool modifiedbol=true;
    while(modifiedbol)
    {
       modifiedbol=false;
       foreach(CATAZJUT::CAtom* atom_s, m_pSupport->atoms())
          foreach(CATAZJUT::CAtom* atom_m, m_pAdsorbMolecule->atoms()){
            distanceCutOff = (atom_s->CovalentRadius() + atom_m->CovalentRadius())*0.6;
            if( m_pPeriodicFramework->distance(atom_m,atom_s) <distanceCutOff ){
                vect = atom_m->position() - atom_s->position();
                vect = (distanceCutOff-vect.norm()+eps)*(vect.normalized());
                this->m_pAdsorbMolecule->moveBy(vect);
                modifiedbol=true;
            }
        }
    }
}
void CModel2DSupport::getGeneValuefromStruct(std::vector<double>& geneRealValue)
{
   geneRealValue.clear();

}
void CModel2DSupport::GeneVARRange(std::vector<GeneVAR>& GeneVARibleVect)
{
    //clear vector!
    GeneVARibleVect.clear();
    //distance of adsorbent on 2d support
    GeneVARibleVect.push_back({1.5,3.0,0.01});
    // a axes displacing
    if(m_latticeDirection==CModel2DSupport::NONE_DIR)
        IdentifyvacuumLayerDirection();

    if(m_latticeDirection == CModel2DSupport::C_AXIS){
        //c a,b
        GeneVARibleVect.push_back({0,m_pPeriodicFramework->unitcell()->a(),0.01});
        GeneVARibleVect.push_back({0,m_pPeriodicFramework->unitcell()->b(),0.01});
    }else if (m_latticeDirection == CModel2DSupport::B_AXIS){
        //b c  a
        GeneVARibleVect.push_back({0,m_pPeriodicFramework->unitcell()->c(),0.01});
        GeneVARibleVect.push_back({0,m_pPeriodicFramework->unitcell()->a(),0.01});
        m_latticeDirection = CModel2DSupport::B_AXIS;
    }else if(m_latticeDirection == CModel2DSupport::A_AXIS){
        //a b  c
        GeneVARibleVect.push_back({0,m_pPeriodicFramework->unitcell()->b(),0.01});
        GeneVARibleVect.push_back({0,m_pPeriodicFramework->unitcell()->c(),0.01});
    }else{
         Log::Error<<"It has no vacuum layer! GeneVARRange_CModel2DSupport!\n";
         boost::throw_exception(std::runtime_error("It has no vacuum layer!! GeneVARRange_CModel2DSupport.\n"));
    }
    //rotation axis
    GeneVARibleVect.push_back({-1,1,0.01});
    GeneVARibleVect.push_back({-1,1,0.01});
    GeneVARibleVect.push_back({-1,1,0.01});

    //rotation angle
    GeneVARibleVect.push_back({0.0,360.0,2.0});


}
void CModel2DSupport::IdentifyvacuumLayerDirection()  //-1: none, 0:a,  1:b  2:c
{

    Eigen::Matrix<double, 3, 3> latticeVect(m_pPeriodicFramework->unitcell()->MatrixOfBravaisLattice());

    double dis;
    bool breakBol=false;
    Vector3 disVect;
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
                m_latticeDirection = CModel2DSupport::A_AXIS;
                break;
            case 1:
                m_latticeDirection = CModel2DSupport::B_AXIS;
                break;
            case 2:
                m_latticeDirection = CModel2DSupport::C_AXIS;
                break;
            default:
                m_latticeDirection = CModel2DSupport::NONE_DIR;
                break;
          }
          break;
       }
    }
}
void CModel2DSupport::perceiveSupportSurface()
{
   // std::vector<size_t> coord_Num;
    int direction;
    if ( m_latticeDirection == CModel2DSupport::C_AXIS)
        direction=2;
    else if(m_latticeDirection == CModel2DSupport::B_AXIS)
        direction=1;
    else //if (m_latticeDirection = CModel2DSupport::A_AXIS)
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
            Log::Error<<" 2D support has serious error! perceiveSupportSurface_CModel2DSupport!";
            boost::throw_exception(std::runtime_error("2D support has serious error! perceiveSupportSurface_CModel2DSupport!"));
        }else
            PlaneCutoff=PlaneCutoff + 0.1;
    }
    Eigen::MatrixXd* currentMat = new (Eigen::MatrixXd)(matom.size(),3);
    for(size_t i=0;i<matom.size();i++)
        (*currentMat).row(i) = matom[i]->position();

    this->m_support_surface = new CATAZJUT::CPlane(currentMat);
    // clear space
    delete currentMat;
}


}
