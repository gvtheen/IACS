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
#include <math.h>
#include "CModelClusterLoaded2DSupport.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "CModelMoleculeAdsorbent.h"
#include "../GaZJUT/GaDeclaration.h"
#include "../CataZJUT/CCartesianCoordinates.h"
#include "../CataZJUT/CFractionCoordinates.h"
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CUnitCell.h"
#include "CModelSupport.h"
#include "../Util/foreach.h"
#include "../CataZJUT/CPlane.h"
#include "../CataZJUT/CSphere.h"
#include "../CataZJUT/CFragment.h"
#include "../CataZJUT/CBondTolerance.h"
#include "../GaZJUT/GaUtilityFunction.h"
#include "../Util/Point-Vector.h"
#include "../CataZJUT/Constant.h"

using util::Log;
using util::Point3;
using util::Vector3;
using CATAZJUT::constants::Pi;

namespace CALCZJUT{

CModelClusterLoaded2DSupport::CModelClusterLoaded2DSupport(CParameter* mParameter,
                                                           CATAZJUT::CConfigurationBase** copy_ppPeriodicFramework,size_t index)
:CModelBase(mParameter,index)
{
    m_pPeriodicFramework = new CATAZJUT::CConfigurationBase(mParameter);
    m_ppBackupPeriodicFramework = copy_ppPeriodicFramework;
    m_latticeDirection = CModelClusterLoaded2DSupport::NONE_DIR;

    m_pClusterOnSupport = nullptr;
    m_pAdsorbMolecule = nullptr;
    m_support_surface = nullptr;
    m_pSphereCluster = nullptr;
}

CModelClusterLoaded2DSupport::~CModelClusterLoaded2DSupport()
{
     if( m_pClusterOnSupport!= nullptr)
        delete m_pClusterOnSupport;
     if( m_pAdsorbMolecule != nullptr)
        delete m_pAdsorbMolecule;
     if( m_support_surface != nullptr)
        delete m_support_surface;
     if( m_pSphereCluster!= nullptr)
        delete m_pSphereCluster;
}

CModelBase* CModelClusterLoaded2DSupport::clone()
{
    CModelClusterLoaded2DSupport* res= new CModelClusterLoaded2DSupport(this->m_pParameter,
                                                                        this->m_ppBackupPeriodicFramework,0);
    res->setPeriodicFramekwork(this->periodicFramework());
    // four types moieties
    res->createMoleAdsorb(this->MoleAdsorbBit());
    res->createSupport(this->SupportBit());
    res->createSupportSurface(this->SupportSurfaceBit());
    res->createSupportedCluster(this->ClusterBit());
    // two types elements
    res->setLatticeDirection(this->latticeDirection());
    res->createClusterSphere();
    res->createSupportPlane();

    return res;
}
void CModelClusterLoaded2DSupport::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{
    Point3 resCoordinate;
    resCoordinate<<realValueOfgene[0],realValueOfgene[1],realValueOfgene[2];
    resCoordinate = m_pSphereCluster->CartesianCoordAtGeneOf(resCoordinate);

     // let adsorbing molecule move to new position
    m_pAdsorbMolecule->moveBy(resCoordinate);

     // rotation some angel along specific axis
    resCoordinate<<realValueOfgene[3],realValueOfgene[4],realValueOfgene[5];
    m_pAdsorbMolecule->rotate(resCoordinate,realValueOfgene[6]);

     // check whether close contacts is avaiable.
    this->eliminateCloseContacts();
}
void CModelClusterLoaded2DSupport::getGeneValuefromStruct(std::vector<double>& currentGeneRealValue)
{

}
void CModelClusterLoaded2DSupport::VarRangeStructRange(std::vector<VarRangeStruct>& currentVarRangeStructible)
{
        double distance = m_support_surface->Distance(m_pSphereCluster->SphereCenter());

        double R = m_pSphereCluster->Radius();

        currentVarRangeStructible.push_back({1.5,3.0,0.01});
        // 2nd gene:    phi: [0,2*PI]
        currentVarRangeStructible.push_back({0.0,2*Pi,0.01});
        // 3rd gene thea: [0,PI]
        currentVarRangeStructible.push_back({0.0,Pi-std::acos(distance/R),0.01});

        // following four genes:  rotation axis and angle of adsorbing molecule over the sphere.

        //rotation axis
        currentVarRangeStructible.push_back({-1.0,1.0,0.01});
        currentVarRangeStructible.push_back({-1.0,1.0,0.01});
        currentVarRangeStructible.push_back({-1.0,1.0,0.01});
        //rotation angle
        currentVarRangeStructible.push_back({0.0,2*Pi,0.01});
}
void CModelClusterLoaded2DSupport::IdentifyvacuumLayerDirection()
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
                m_latticeDirection = CModelClusterLoaded2DSupport::A_AXIS;
                break;
            case 1:
                m_latticeDirection = CModelClusterLoaded2DSupport::B_AXIS;
                break;
            case 2:
                m_latticeDirection = CModelClusterLoaded2DSupport::C_AXIS;
                break;
            default:
                m_latticeDirection = CModelClusterLoaded2DSupport::NONE_DIR;
                break;
          }
          break;
       }
    }
}
void CModelClusterLoaded2DSupport::perceiveSupportSurface()
{
    size_t direction;

    if(m_latticeDirection==CModelClusterLoaded2DSupport::NONE_DIR)
        IdentifyvacuumLayerDirection();

    // convert the direction to coordinate index;
    if ( m_latticeDirection == CModelClusterLoaded2DSupport::C_AXIS)
        direction=2;
    else if(m_latticeDirection == CModelClusterLoaded2DSupport::B_AXIS)
        direction=1;
    else //if (m_latticeDirection = CModel2DSupport::A_AXIS)
        direction=0;

    // search the most_highest atom along the vacuumLayer direction

    double max_latticedirection=DBL_MIN;
    CATAZJUT::CAtom *most_highest_atom;

    Bitset  tempBit(m_pPeriodicFramework->atomCount(),false);

    foreach(CATAZJUT::CAtom *atom,m_pSupport->atoms())
        if(atom->position()[direction]>max_latticedirection){
            max_latticedirection=atom->position()[direction];
            most_highest_atom = atom;
        }
    tempBit.set(most_highest_atom->index(),true);
    //std::string atomStr =  most_highest_atom->Symbol();

    //std::vector<CATAZJUT::CAtom*> clusterAtom;
    double monitor_height=max_latticedirection;
    while(true){
        monitor_height = monitor_height - 0.2;
        foreach(const CATAZJUT::CAtom *atom,m_pSupport->atoms())
             if(atom->position()[direction]>monitor_height)
                tempBit.set(atom->index(),true);
        if(this->m_pPeriodicFramework->FromMoietyGetFragmentsAtom(tempBit,most_highest_atom)==true){
            this->m_BitbackupCluster = tempBit;
            break;
        }
    }
    this->createSupportedCluster(this->m_BitbackupCluster);

    m_BitBackupSupportNoCluster = ~(m_BitbackupCluster & m_BitbackupSupport);

    this->createSupportSurface(m_BitBackupSupportNoCluster);
    //construct the mode sphere of the cluster on the support.
    this->createClusterSphere();

    // identify the pure surface on the support
    this->createSupportPlane();

}
void CModelClusterLoaded2DSupport::eliminateCloseContacts(double distanceCutOff)
{
    Vector3 vect;
    double eps=0.2;
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
void CModelClusterLoaded2DSupport::createSupport(const Bitset& mht )
{
    assert(m_pPeriodicFramework);
    m_BitbackupSupport = mht;
    m_pSupport= new CModelSupport(this->m_pPeriodicFramework,m_BitbackupSupport);
}
void CModelClusterLoaded2DSupport::createMoleAdsorb(const Bitset& mht)
{
    assert(m_pPeriodicFramework);
    m_BitbackupAdsorbMolecule = mht;
    m_pAdsorbMolecule = new CModelMoleculeAdsorbent(this->m_pPeriodicFramework,m_BitbackupAdsorbMolecule);
}
void CModelClusterLoaded2DSupport::createSupportedCluster(const Bitset& mht)
{
    assert(this->m_pPeriodicFramework);
    m_BitbackupCluster = mht;
    m_pClusterOnSupport = new CModelMoleculeAdsorbent(m_pPeriodicFramework,m_BitbackupCluster);

}
void CModelClusterLoaded2DSupport::createSupportSurface(const Bitset& mht)
{
    assert(this->m_pPeriodicFramework);
    m_BitBackupSupportNoCluster = mht;
    m_pPureSupportSurface = new CModelSupport( m_pPeriodicFramework,m_BitBackupSupportNoCluster );
}
void CModelClusterLoaded2DSupport::createClusterSphere()
{
    Eigen::MatrixXd*  clusterCoordinate=new (Eigen::MatrixXd)(m_pPureSupportSurface->atomCount(),3);
    size_t num=0;
    foreach(CATAZJUT::CAtom* atom,m_pPureSupportSurface->atoms())
        clusterCoordinate->row(num++)=atom->position();
    this->m_pSphereCluster = new CATAZJUT::CSphere(clusterCoordinate);
    delete clusterCoordinate;
}
void CModelClusterLoaded2DSupport::createSupportPlane()
{
    size_t direction;
    double max_latticedirection=DBL_MIN;
    if ( m_latticeDirection == CModelClusterLoaded2DSupport::C_AXIS)
        direction=2;
    else if(m_latticeDirection == CModelClusterLoaded2DSupport::B_AXIS)
        direction=1;
    else //if (m_latticeDirection = CModel2DSupport::A_AXIS)
        direction=0;
    foreach(const CATAZJUT::CAtom *atom,m_pPureSupportSurface->atoms()){
        if(atom->position()[direction]>max_latticedirection)
            max_latticedirection=atom->position()[direction];
    }
    std::vector<CATAZJUT::CAtom*> matom;
    double PlaneCutoff = 0.1;
    while(true){
        matom.clear();
        foreach(CATAZJUT::CAtom *atom,m_pPureSupportSurface->atoms()){
           if(std::fabs(atom->position()[direction]-max_latticedirection)<PlaneCutoff)
              matom.push_back(atom);
        }
        if(matom.size()>4)
            break;
        else if(PlaneCutoff >1.0){
            Log::Error<<" 2D support has serious error! perceiveSupportSurface_CModelClusterLoaded2DSupport::perceiveSupportSurface!";
            boost::throw_exception(std::runtime_error("2D support has serious error! perceiveSupportSurface_CModelClusterLoaded2DSupport::perceiveSupportSurface!"));
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
Bitset CModelClusterLoaded2DSupport::SupportSurfaceBit()
{
    return this->m_BitBackupSupportNoCluster;
}
Bitset CModelClusterLoaded2DSupport::ClusterBit()
{
    return this->m_BitbackupCluster;
}
Bitset CModelClusterLoaded2DSupport::SupportBit()
{
      return this->m_BitbackupSupport;
}
Bitset CModelClusterLoaded2DSupport::MoleAdsorbBit()
{
      return this->m_BitbackupAdsorbMolecule;
}
        //overload the father class's function
void CModelClusterLoaded2DSupport::setPeriodicFramekwork(CATAZJUT::CConfigurationBase* mbf)
{
     if(this->m_pPeriodicFramework!=nullptr)
        delete this->m_pPeriodicFramework;
     this->m_pPeriodicFramework=new CATAZJUT::CConfigurationBase(*mbf);
     this->m_pAdsorbMolecule->setConfiguration(this->m_pPeriodicFramework);
     this->m_pSupport->setConfiguration(this->m_pPeriodicFramework);
}
CModelClusterLoaded2DSupport::LATT_DIRECTION CModelClusterLoaded2DSupport::latticeDirection()
{
    return this->m_latticeDirection;
}
void CModelClusterLoaded2DSupport::setLatticeDirection(CModelClusterLoaded2DSupport::LATT_DIRECTION mht)
{
    this->m_latticeDirection = mht;
}



}   //namespace
