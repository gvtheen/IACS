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
#include "CModelClusterSupport.h"
#include "CModelSupport.h"
#include "CModelMoleculeAdsorbent.h"
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CCrystalPlanes.h"
#include "../CataZJUT/CSphere.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../CataZJUT/Constant.h"
#include "../Util/Bitset.h"
#include "../Util/foreach.h"
#include "../Util/Point-Vector.h"

using util::Bitset;
using util::Log;
using util::Matrix;
using util::Vector3;
using CATAZJUT::constants::Pi;

namespace CALCZJUT{

CModelClusterSupport::CModelClusterSupport(CParameter* mPara,
                                           CATAZJUT::CConfigurationBase** copy_ppPeriodicFramework,
                                           size_t index)
:CModelBase(mPara,index)
{
    m_pPeriodicFramework = new CATAZJUT::CConfigurationBase(mPara);
    m_ppBackupPeriodicFramework = copy_ppPeriodicFramework;
    m_pSupport = nullptr;
    m_pAdsorbMolecule = nullptr;
    m_pCrystalPlanes = nullptr;
    m_ClusterModelPerceived = false;
    //m_ClusterModelType=CModelClusterSupport::POLYHEDRON;
}
CModelClusterSupport::CModelClusterSupport(CModelClusterSupport& obj)
:CModelBase(obj)
{
//   #ifdef DEBUG
//        Log::Debug<<"CModelClusterSupport::CModelClusterSupport()" << std::endl;
//   #endif // DEBU
   this->createSupport(obj.SupportBit());
   this->createMoleAdsorb(obj.MoleAdsorbBit());
//   #ifdef DEBUG
//        Log::Debug<<"1-CModelClusterSupport::CModelClusterSupport()" << std::endl;
//   #endif // DEBU
   if(obj.crystalPlanes()!=nullptr)
     this->setCrystalPlanes(obj.crystalPlanes());
   this->setRandomInitState(obj.RandomInitState());

//   #ifdef DEBUG
//        Log::Debug<<"2-CModelClusterSupport::CModelClusterSupport()" << std::endl;
//   #endif // DEBU
}
CModelBase* CModelClusterSupport::clone()
{
      return new CModelClusterSupport(*this);
}
CModelClusterSupport::~CModelClusterSupport()
{
    //dtor
}

void CModelClusterSupport::createSupport(const Bitset& mht)
{
   assert(m_pPeriodicFramework);

   m_BitbackupSupport = mht;
   m_pSupport= new CModelSupport(this->m_pPeriodicFramework,m_BitbackupSupport);
}
void CModelClusterSupport::createMoleAdsorb(const Bitset& mht)
{
   assert(m_pPeriodicFramework);

   m_BitbackupAdsorbMolecule = mht;
   m_pAdsorbMolecule = new CModelMoleculeAdsorbent(this->m_pPeriodicFramework,m_BitbackupAdsorbMolecule);
}
// crystal plane is set and get.
CATAZJUT::CCrystalPlanes* CModelClusterSupport::crystalPlanes()
{
   return this->m_pCrystalPlanes;
}
void CModelClusterSupport::setCrystalPlanes(CATAZJUT::CCrystalPlanes* mth)
{
   if(this->m_pCrystalPlanes!=nullptr)
       delete this->m_pCrystalPlanes;
   this->m_pCrystalPlanes = mth->Clone();
}
void CModelClusterSupport::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{
   Point3 resCoordinate;

   if(this->RandomInitState()==false){
        this->setRandomInitState(true);
        goto RETURN_Label;
    }
    if(realValueOfgene.size()==0)
        goto RETURN_Label;

   if(m_ClusterModelType==CModelClusterSupport::POLYHEDRON){
     //get coordinate of gravity center of adsorbing molecule.
     resCoordinate = this->m_pCrystalPlanes->CartesianCoordinateAtGene((int)realValueOfgene[0],realValueOfgene[1],
                                                                       realValueOfgene[2],realValueOfgene[3]);
     // let adsorbing molecule move to new position
     this->m_pAdsorbMolecule->moveBy(resCoordinate);
     // rotation some angel along specific axis
     resCoordinate<<realValueOfgene[4],realValueOfgene[5],realValueOfgene[6];
     this->m_pAdsorbMolecule->rotate(resCoordinate,realValueOfgene[7]);
     // check whether close contacts is avaiable.
     this->eliminateCloseContacts();
   }else{  //CModelClusterSupport::SPHERE
     //get coordinate of gravity center of adsorbing molecule.
     resCoordinate<<realValueOfgene[0],realValueOfgene[1],realValueOfgene[2];
     resCoordinate = this->m_pSphere->CartesianCoordAtGeneOf(resCoordinate);

     // let adsorbing molecule move to new position
     this->m_pAdsorbMolecule->moveBy(resCoordinate);

     // rotation some angel along specific axis
     resCoordinate<<realValueOfgene[3],realValueOfgene[4],realValueOfgene[5];
     this->m_pAdsorbMolecule->rotate(resCoordinate,realValueOfgene[6]);

     // check whether close contacts is avaiable.
     this->eliminateCloseContacts();
   }
RETURN_Label:
    ;// Nothing is done.
}
void CModelClusterSupport::getGeneValuefromStruct(std::vector<double>& currentGeneRealValue)
{

}
void CModelClusterSupport::VarRangeStructRange(std::vector<VarRangeStruct>& currentVarRangeStructible)
{
   if(currentVarRangeStructible.size()!=0)
       currentVarRangeStructible.clear();  // clear all data of in currentVarRangeStructible

   if( !this->m_ClusterModelPerceived )
        perceiveClusterModel();

   if(m_ClusterModelType==CModelClusterSupport::POLYHEDRON){
    // index of crystal planes in the cluster.                    1st gene descriptor.
    // max is m_pCrystalPlanes->crystalPlaneNum() - 1; substracting 1 is necessary!
        currentVarRangeStructible.push_back({0,m_pCrystalPlanes->crystalPlaneNum()-1.0,0.5});

        // Height the adsorbing molecule on the index crystal plane.  2nd gene descriptor.
        currentVarRangeStructible.push_back({1.5,3.0,0.01});
        // ratio
        currentVarRangeStructible.push_back({0.0,1.0,0.01});
        // angle
        currentVarRangeStructible.push_back({0.0,2*Pi,0.01});

        //rotation axis
        currentVarRangeStructible.push_back({-1.0,1.0,0.01});
        currentVarRangeStructible.push_back({-1.0,1.0,0.01});
        currentVarRangeStructible.push_back({-1.0,1.0,0.01});
        //rotation angle
        currentVarRangeStructible.push_back({0.0,2*Pi,0.01});
   }else{  //CModelClusterSupport::SPHERE
        // 1st gene:  height: the adsorbing molecule on the index crystal plane.
        currentVarRangeStructible.push_back({1.5,3.0,0.01});
        // 2nd gene:    phi: [0,2*PI]
        currentVarRangeStructible.push_back({0.0,2*Pi,0.01});
        // 3rd gene thea: [0,PI]
        currentVarRangeStructible.push_back({0.0,Pi,0.01});

        // following four genes:  rotation axis and angle of adsorbing molecule over the sphere.
        //rotation axis
        currentVarRangeStructible.push_back({-1.0,1.0,0.01});
        currentVarRangeStructible.push_back({-1.0,1.0,0.01});
        currentVarRangeStructible.push_back({-1.0,1.0,0.01});
        //rotation angle
        currentVarRangeStructible.push_back({0.0,2*Pi,0.01});
   }

}
void CModelClusterSupport::perceiveClusterModel()
{
//    #ifdef DEBUG
//                     Log::Debug<<"CModelClusterSupport::perceiveClusterModel()" << std::endl;
//    #endif // DEBU
    if(m_ClusterModelType==CModelClusterSupport::POLYHEDRON){
        std::map<std::string,size_t> maxCoordNum;
        this->m_pPeriodicFramework->maxCoordinationNum(maxCoordNum);

        //Eigen::MatrixXd *atom_Coordinate =  new (Eigen::MatrixXd)(m_pSupport->atomCount(),3);
        std::vector<size_t> atom_Coordinate;

        Point3 tempP;
        #ifdef DEBUG
            Log::Debug<<"21-CModelClusterSupport::perceiveClusterModel()" << std::endl;
        #endif // DEBU
        foreach(CATAZJUT::CAtom* atom_s,this->m_pSupport->atoms())
           if(atom_s->bondCount()<maxCoordNum[atom_s->Symbol()])
              atom_Coordinate.push_back(atom_s->index());

        if(atom_Coordinate.size()<4){
           foreach(CATAZJUT::CAtom* atom_s,this->m_pSupport->atoms())
                if(atom_s->bondCount()==maxCoordNum[atom_s->Symbol()])
                   atom_Coordinate.push_back(atom_s->index());
        }
        #ifdef DEBUG
            Log::Debug<<"22-CModelClusterSupport::perceiveClusterModel():" <<atom_Coordinate.size()<< std::endl;
        #endif // DEBU
        this->m_pCrystalPlanes =  new CATAZJUT::CCrystalPlanes(atom_Coordinate,this->m_pPeriodicFramework);
        this->m_pCrystalPlanes->CreateCrystalPlane();
        this->m_pCrystalPlanes->outputCrystalPlane();

    }else{  //CModelClusterSupport::SPHERE
        Eigen::MatrixXd* coordinatePoints = new Eigen::MatrixXd(m_pSupport->atomCount(),3);
        #ifdef DEBUG
                     Log::Debug<<"2-CModelClusterSupport::perceiveClusterModel()" << std::endl;
        #endif // DEBU
        size_t num=0;
        foreach(CATAZJUT::CAtom* atom_s,this->m_pSupport->atoms())
           coordinatePoints->row(num++)=atom_s->position();

        #ifdef DEBUG
                     Log::Debug<<"3-CModelClusterSupport::perceiveClusterModel()" << std::endl;
                     Log::Debug<<coordinatePoints->rows()<<std::endl;
        #endif // DEBU
        this->m_pSphere = new CATAZJUT::CSphere(coordinatePoints);
        #ifdef DEBUG
                     Log::Debug<<"4-CModelClusterSupport::perceiveClusterModel()" << std::endl;
        #endif // DEBU
        this->m_pSphere->CreateSphere();
        this->m_pSphere->output();
        delete coordinatePoints;
    }
    #ifdef DEBUG
                     Log::Debug<<"4-CModelClusterSupport::perceiveClusterModel()" << std::endl;
    #endif // DEBU

}
void CModelClusterSupport::eliminateCloseContacts(double distanceCutOff)
{
    Vector3 vect;
    double eps=0.01;
    bool modifiedbol=true;
    size_t modify_num=0;
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
                modify_num++;
            }
        }
       if(modify_num>256){
         Log::Warn<<"Cycle of function is more than 65536 in CModelClusterSupport::eliminateCloseContacts!"<<std::endl;
         break;
       }
    }
}
Bitset CModelClusterSupport::SupportBit()
{
    return this->m_BitbackupSupport;
}
Bitset CModelClusterSupport::MoleAdsorbBit()
{
    return this->m_BitbackupAdsorbMolecule;
}



}
