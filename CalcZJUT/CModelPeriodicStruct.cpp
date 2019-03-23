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
#include "CModelPeriodicStruct.h"
#include "CParameter.h"
#include "../Util/log.hpp"
#include "../Util/foreach.h"
#include "../Util/Bitset.h"
#include "../Util/Point-Vector.h"
#include "../Util/utilFunction.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CFragment.h"

using util::Log;
using util::Point3;
using util::Vector4;

namespace CALCZJUT{

CModelPeriodicStruct::CModelPeriodicStruct(CParameter* mPara,
                                           size_t index)
:CModelBase(mPara,index)
{
     this->m_pPeriodicFramework = new CATAZJUT::CConfigurationBase(mPara);
}

CModelPeriodicStruct::~CModelPeriodicStruct()
{
    delete this->m_pPeriodicFramework;
}

CModelBase* CModelPeriodicStruct::clone()
{
    CModelPeriodicStruct* res =  new CModelPeriodicStruct (this->m_pParameter,0);

    res->m_VarRangeStruct.assign(this->m_VarRangeStruct.begin(),this->m_VarRangeStruct.end());
    res->m_IsNeedRandomInit = this->RandomInitState();
    res->setChemicalFormula(this->chemicalFormula());

    return res;
}
void CModelPeriodicStruct::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{
     size_t index=0;
     Point3 tempPoint;

     if( realValueOfgene.size()%3 != 0 ){
        Log::Error<<"The size() of realValueofGene is error! CModelPeriodicStruct::setGeneValueToStruct!\n";
        boost::throw_exception(std::runtime_error("The size() of realValueofGene is error! CModelPeriodicStruct::setGeneValueToStruct!"));
     }
     if(realValueOfgene.size()==0)
        goto RETURN_Label;


      if(this->m_pParameter->currentGenerationNum()>0){
          try{
             foreach(CATAZJUT::CAtom* atom_s, this->m_pPeriodicFramework->atoms()){
                tempPoint<<realValueOfgene[index],realValueOfgene[index+1],realValueOfgene[index+2];
                atom_s->SetPosition(tempPoint);
                index += 3;
                if( index >= realValueOfgene.size() ) break;
            }
          }catch(std::exception const &e){
             Log::Error<<e.what();
             boost::throw_exception(std::runtime_error(e.what()));
          }
      }

RETURN_Label:
    ;
}
void CModelPeriodicStruct::getGeneValuefromStruct(std::vector<double>& currentGeneRealValue)
{
    if(currentGeneRealValue.size()!=0)
       currentGeneRealValue.clear();
    Point3 tempPoint;
    foreach(CATAZJUT::CAtom* atom_s, this->m_pPeriodicFramework->atoms()){
       tempPoint = atom_s->position();
       currentGeneRealValue.push_back(tempPoint(0,0));
       currentGeneRealValue.push_back(tempPoint(1,0));
       currentGeneRealValue.push_back(tempPoint(2,0));
    }
}
void CModelPeriodicStruct::VarRangeStructRange(std::vector<VarRangeStruct>& currentVarRangeStructible)
{
   //
}
void CModelPeriodicStruct::eliminateCloseContacts(CATAZJUT::CConfigurationBase* curr_struct,
                                                  double distanceCutOff)
{
    util::Vector3 vect;
    util::Point3  center_P;
    double eps=0.1;
    bool modifiedbol=true;
    size_t less_cycle=0;

    while(modifiedbol)
    {
       modifiedbol=false;
       // get center pointer of structure
       center_P = curr_struct->center();

       foreach(CATAZJUT::CAtom* atom_s, curr_struct->atoms())
          foreach(CATAZJUT::CAtom* atom_m, curr_struct->atoms()){
            distanceCutOff = (atom_s->CovalentRadius() + atom_m->CovalentRadius())*0.6;
            if( curr_struct->distance(atom_m,atom_s) < distanceCutOff ){
                /*
                  Let this atom move far away from center atom.
                */
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
       if(less_cycle > 65536){
          Log::Warn<<"Cycle of function is more than 65536 in CModelCluster::eliminateCloseContacts!"<<std::endl;
          break;
       }
    }

}
void CModelPeriodicStruct::eliminateFragment(CATAZJUT::CConfigurationBase* curr_struct)
{
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
        }
        //
        Vector4 mainsphereEquation4, othersphereEquation4;
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
               // further judge whether two fragments is bonded.
               less_cycle=0;
               while(mainFragment->isBondTo(fragment_s)!=true){
                   differVect= 0.2*(maincenter-othercenter).normalized();
                   fragment_s->move(differVect);
                   less_cycle++;
                   if(less_cycle>65536){
                      Log::Warn<<"Cycle of function is more than 65536 in CModelCluster::eliminateFragment!"<<std::endl;
                      break;
                   }
               }
            }
    }
}
std::vector<std::pair<std::string,size_t>>& CModelPeriodicStruct::chemicalFormula()
{
    if(this->m_chemicalFormula.size()==0)
        this->m_chemicalFormula = this->m_pPeriodicFramework->composition();

    return this->m_chemicalFormula;
}
void CModelPeriodicStruct::setChemicalFormula(const std::vector<std::pair<std::string,size_t>>& mth)
{
   this->m_chemicalFormula.assign(mth.begin(),mth.end());
}


}
