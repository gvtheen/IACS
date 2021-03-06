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
#include <algorithm>
#include <boost/scoped_ptr.hpp>
#include <boost/make_shared.hpp>
#include <sstream>
#include "CConfigurationBase.h"
#include "CConfigurationBase_sub.h"
#include "../Util/foreach.h"
#include "CAtom.h"
#include "CBond.h"
#include "CElement.h"
#include "CFragment.h"
#include "CUnitCell.h"
#include "CMolecularSymmetry.h"
#include "CSpaceGroup.h"
#include "CInternalCoordinates.h"
#include "CCartesianCoordinates.h"
#include "CFractionCoordinates.h"
#include "CBondTolerance.h"
#include "CConfigurationPrivateData.h"
#include "../Util/log.hpp"


/*
The target of this Configurationbase class:
   (1) keep the date of atoms, bonds and et.al except for coordinate.
       Note: the coordinate date is saved in the object of individual coordinate system.
   (2) the operation(such as read,write,remove)of atom
   (3) adjust coordinate by using coordinate pointer.
   (4) analyze fragments of whole configuration
   (5) control the format of output and input
   (3) build the bridge with work-modes(such as cluster-GA, cluster-molecule-GA and 2D-structure-molecule-GA)
*/
using util::Log;

namespace CATAZJUT{
/*
definition of class CConfigurationPrivateData

*/
CConfigurationPrivateData::CConfigurationPrivateData()
{
    fragmentsPerceived=false;
}
CConfigurationPrivateData::~CConfigurationPrivateData()
{
    for(size_t i=0;i<this->bonds.size();i++)
        delete bonds[i];
        bonds.clear();

    for(size_t i=0;i<this->fragments.size();i++)
        delete fragments[i];
        fragments.clear();

    for(size_t i=0;i<this->atomBonds.size();i++)
    {
        for(size_t j=0;j<atomBonds[i].size();j++)
            delete atomBonds[i][j];
        atomBonds[i].clear();
    }

    for(size_t i=0;i<this->bondAtoms.size();i++)
    {
        delete bondAtoms[i].first;
        delete bondAtoms[i].second;
    }
    bondAtoms.clear();
}
/*
definition of class CConfigurationBase

*/
CConfigurationBase::CConfigurationBase(CALCZJUT::CParameter* mPara)
:m_pData(new CConfigurationPrivateData()),m_pParameter(mPara)
{
     this->m_pCartesian=nullptr;
     m_pBondEvaluator = new CBondTolerance(this);

     m_pBondEvaluator->setTolerancefactor(m_pParameter->bondToleranceFactor);
     m_pBondEvaluator->setExcludeBond(m_pParameter->excludeBond);
     this->m_CoordinateType  = CATAZJUT::DEFINED::Cartesian;
     this->m_DimensionalType = CATAZJUT::DEFINED::Periodic;
}

CConfigurationBase::CConfigurationBase(CConfigurationBase& other)
:m_pData(new CConfigurationPrivateData())
{
     /*
        std::vector<CElement>      m_Element;
        std::vector<CAtom*>        m_Atom;
        CBondTolerance            *m_pBondEvaluator;

        Three variables were set during the following process addatom, addbond;
     */
//     #ifdef DEBUG
//      Log::Debug<<"***********CConfigurationBase::CConfigurationBase(CConfigurationBase& other)***********"<< std::endl;
//     #endif
     this->m_pCartesian=nullptr;
     m_pBondEvaluator = new CBondTolerance(this);
     this->m_pParameter=other.sysParameter();
     m_pBondEvaluator->setTolerancefactor(m_pParameter->bondToleranceFactor);
     m_pBondEvaluator->setExcludeBond(m_pParameter->excludeBond);
//     #ifdef DEBUG
//      Log::Debug<<"***********2-CConfigurationBase::CConfigurationBase(CConfigurationBase& other)***********"<< std::endl;
//     #endif
     if(other.atomCount()!=0){
         this->m_CoordinateType = other.coordinateType();
         this->m_DimensionalType= other.dimensionalType();

         m_pData->name = other.m_pData->name;
         std::map<const CAtom*,CAtom*> oldToNew;
         foreach(CAtom* atom, other.atoms())
              oldToNew[atom]=this->addAtomCopy(atom);
//         #ifdef DEBUG
//           Log::Debug<<"***********21-CConfigurationBase::CConfigurationBase(CConfigurationBase& other)***********"<< std::endl;
//         #endif
         foreach(CBond* bond,other.bonds())
         {
             //CBone* newBond = addBond(oldToNew[bond->atom1()],oldToNew[bond->atom2()]);
             this->addBond(oldToNew[bond->atom1()],oldToNew[bond->atom2()]);
         }
//         #ifdef DEBUG
//            Log::Debug<<"***********22-CConfigurationBase::CConfigurationBase(CConfigurationBase& other)***********"<< std::endl;
//         #endif
         this->constraintBits = other.constraintBit();
         oldToNew.clear();
     }

//     #ifdef DEBUG
//      Log::Debug<<"***********3-CConfigurationBase::CConfigurationBase(CConfigurationBase& other)***********"<< std::endl;
//     #endif
}
CConfigurationBase* CConfigurationBase::clone()
{
    return new CConfigurationBase(*this);
}


CConfigurationBase::~CConfigurationBase()
{
    foreach(CAtom* atom, m_Atom)
         delete atom;

    foreach(CBond* bond,m_pData->bonds)
         delete bond;

    foreach(CFragment* fragment,m_pData->fragments)
         delete fragment;
//    foreach( boost::shared_ptr<CCoordinateSet> coordinateSet,m_pData->coordinateSets)
//    {
//          delete coordinateSet;
//    }
    delete m_pBondEvaluator;
    delete m_pCartesian;
    delete m_pData;
}

std::vector<std::pair<std::string,size_t>>& CConfigurationBase::composition()
{
    if(this->m_Composition.size()!=0)
        return this->m_Composition;

    std::map<unsigned char,size_t>* comp = new (std::map<unsigned char,size_t>);
    foreach(const CAtom* atom,atoms())
        (*comp)[atom->AtomNumber()]++;
    //default, map is sorted by the key ;

    for(std::map<unsigned char,size_t>::iterator it=comp->begin();it!=comp->end();it++)
    {
        std::pair<std::string,size_t>  mPair;
        mPair.first  = (new CElement(it->first))->symbol();
        mPair.second = it->second;
        this->m_Composition.push_back(mPair);
    }
    delete comp;

    return this->m_Composition;

}
std::string CConfigurationBase::formula()
{
    std::vector<std::pair<std::string,size_t>> composit=composition();
    std::stringstream formula;

    std::vector<std::pair<std::string,size_t>>::iterator it;
    for(it=composit.begin();it!=composit.end();it++)
    {
        formula << it->first;
        if(it->second > 1)
           formula << it->second;
    }
    return formula.str();
}
//return cartesian coordinate of geometry
CCartesianCoordinates* CConfigurationBase::coordinates()
{
    boost::shared_ptr<CCoordinateSet> pSP;
    if(m_pData->coordinateSets.empty()|| m_pData->coordinateSets.front()->type()\
                                               == CATAZJUT::DEFINED::None){
            // if no, construct new coordinate
        m_pCartesian= new CCartesianCoordinates(this);
        //std::cout<<"size of atomnum:"<<this->atomCount()<<std::endl;
        m_pData->coordinateSets.push_back(boost::make_shared<CCoordinateSet>(m_pCartesian));
    }else{
        const boost::shared_ptr<CCoordinateSet> &coordinateSet = m_pData->coordinateSets.front();
        switch(coordinateSet->type()){
            case CATAZJUT::DEFINED::Cartesian:
                m_pCartesian = coordinateSet->cartesianCoordinates();
                break;
            case CATAZJUT::DEFINED::Internal:
                pSP=this->coordinateSet(CATAZJUT::DEFINED::Cartesian);
                if(pSP != nullptr)
                    this->removeCoordinateSet(pSP);
//                    #ifdef DEBUG
//                       Log ::Debug<<"2-CConfigurationBase::coordinates() "<<std::endl;
//                    #endif
                m_pCartesian = coordinateSet->internalCoordinates()->toCartesianCoordinates();
                m_pData->coordinateSets.push_back(boost::make_shared<CCoordinateSet>(m_pCartesian));
                break;
            case CATAZJUT::DEFINED::Fraction:
                pSP=this->coordinateSet(CATAZJUT::DEFINED::Cartesian);
                if(pSP != nullptr)
                    this->removeCoordinateSet(pSP);
//                    #ifdef DEBUG
//                       Log ::Debug<<"3-CConfigurationBase::coordinates() "<<std::endl;
//                    #endif
                m_pCartesian = coordinateSet->fractionCoordinates()->toCartesianCoordinates();
                m_pData->coordinateSets.push_back(boost::make_shared<CCoordinateSet>(m_pCartesian));

                break;
            case CATAZJUT::DEFINED::None:
                break;
        }
//           #ifdef DEBUG
//              Log ::Debug<<"4-CConfigurationBase::coordinates() "<<std::endl;
//           #endif
       this->stickCoordinateSet(CATAZJUT::DEFINED::Cartesian);
    }
    return this->m_pCartesian;
}
CInternalCoordinates* CConfigurationBase::Internalcoordinates()
{
    //int returnindex=-1;
    if(m_pData->coordinateSets.empty()|| m_pData->coordinateSets.front()->type()\
                                               == CATAZJUT::DEFINED::None){
            // if no, construct new coordinate
            CInternalCoordinates* temp_Internal= new CInternalCoordinates(this);
            m_pData->coordinateSets.push_back(boost::make_shared<CCoordinateSet>(temp_Internal));
//            returnindex = m_pData->coordinateSets.size() - 1;
//            return m_pData->coordinateSets[returnindex]->internalCoordinates();
    }else{
        boost::shared_ptr<CCoordinateSet> pSP = this->coordinateSet(CATAZJUT::DEFINED::Internal);

        if(pSP == nullptr){
            CInternalCoordinates* temp_Internal= new CInternalCoordinates(this,atomCount());
            m_pData->coordinateSets.push_back(boost::make_shared<CCoordinateSet>(temp_Internal));
//            returnindex = m_pData->coordinateSets.size() - 1;
//            return m_pData->coordinateSets[returnindex]->internalCoordinates();
        }//else
           // return pSP->internalCoordinates();
    }
    // push Internal coordinateset TOP
    this->stickCoordinateSet(CATAZJUT::DEFINED::Internal);

    return m_pData->coordinateSets.front()->internalCoordinates();
}
CFractionCoordinates*  CConfigurationBase::Fractioncoordinates()
{
    //int returnindex=-1;
    CFractionCoordinates* temp_Internal=nullptr;
    if(m_pData->coordinateSets.empty()|| m_pData->coordinateSets.front()->type()\
                                               == CATAZJUT::DEFINED::None){
            // if no, construct new coordinate
        temp_Internal = new CFractionCoordinates(this);
        m_pData->coordinateSets.push_back(boost::make_shared<CCoordinateSet>(temp_Internal));
//        returnindex = m_pData->coordinateSets.size() - 1;
//        return m_pData->coordinateSets[returnindex]->fractionCoordinates();
    }else{
        boost::shared_ptr<CCoordinateSet> pSP = this->coordinateSet(CATAZJUT::DEFINED::Fraction);
        if(pSP == nullptr)
        {
            temp_Internal= new CFractionCoordinates(this);
            m_pData->coordinateSets.push_back(boost::make_shared<CCoordinateSet>(temp_Internal));
           // returnindex = m_pData->coordinateSets.size() - 1;
           // return m_pData->coordinateSets[returnindex]->fractionCoordinates();
        }//else
           // return pSP->fractionCoordinates();
    }
    this->stickCoordinateSet(CATAZJUT::DEFINED::Fraction);

    return m_pData->coordinateSets.front()->fractionCoordinates();
//    #ifdef DEBUG
//        Log ::Debug<<"CConfigurationBase::Fractioncoordinates() "<<std::endl;
//    #endif
}
//coordinate property: Cartesian,Internal,Fraction
/*

*/
void CConfigurationBase::addCoordinateSet(const boost::shared_ptr<CCoordinateSet> &coordinates)
{
     m_pData->coordinateSets.push_back(coordinates);
}
void CConfigurationBase::addCoordinateSet(CCartesianCoordinates *coordinates)
{
     addCoordinateSet(boost::make_shared<CCoordinateSet>(coordinates));
}
void CConfigurationBase::addCoordinateSet(CInternalCoordinates *coordinates)
{
     addCoordinateSet(boost::make_shared<CCoordinateSet>(coordinates));
}
void CConfigurationBase::addCoordinateSet(CFractionCoordinates *coordinates)
{
     addCoordinateSet(boost::make_shared<CCoordinateSet>(coordinates));
}
bool CConfigurationBase::removeCoordinateSet(const boost::shared_ptr<CCoordinateSet> &coordinates)
{
    std::vector<boost::shared_ptr<CCoordinateSet>>::iterator it = \
                         std::find(m_pData->coordinateSets.begin(),m_pData->coordinateSets.end(),coordinates);
    if(it != m_pData->coordinateSets.end()){
        m_pData->coordinateSets.erase(it);
        return true;
    }
    return false;
}
boost::shared_ptr<CCoordinateSet> CConfigurationBase::coordinateSet(size_t index) const
{
    assert(m_pData->coordinateSets.size()>index);
    return m_pData->coordinateSets[index];
}
boost::shared_ptr<CCoordinateSet> CConfigurationBase::coordinateSet(CATAZJUT::DEFINED::CoordinateType type)
{

    foreach(const boost::shared_ptr<CCoordinateSet> &coordinates,m_pData->coordinateSets)
        if(coordinates->type() == type)
            return coordinates;
    return boost::shared_ptr<CCoordinateSet>();
}
int CConfigurationBase::coordinateSetIndex(CATAZJUT::DEFINED::CoordinateType type)const
{
    int res=-1;
    for(size_t i=0;i<this->m_pData->coordinateSets.size();i++)
        if(m_pData->coordinateSets[i]->type()==type){
            res=i;
            break;
        }
    return res;
}
CConfigurationBase::CoordinateSetRange CConfigurationBase::coordinateSets() const
{
    return boost::make_iterator_range(m_pData->coordinateSets.begin(),
                                      m_pData->coordinateSets.end());
}
// push the coordinateset of type the top
void CConfigurationBase::stickCoordinateSet(CATAZJUT::DEFINED::CoordinateType type)
{
    int index=this->coordinateSetIndex(type);

    if(index==-1){
         Log::Error<<"None is coordinate type! CConfigurationBase::stickCoordinateSet!"<<std::endl;
         boost::throw_exception(std::runtime_error("None is coordinate type! CConfigurationBase::stickCoordinateSet!\n"));
    }
    boost::shared_ptr<CCoordinateSet>tempSet=m_pData->coordinateSets.front();
    if(tempSet->type()!=type){
       m_pData->coordinateSets[0]=this->coordinateSet(type);
       m_pData->coordinateSets[index]=tempSet;
    }
}
//atom property
/*

*/
void CConfigurationBase::addAtom(const std::string& elementName, const Point3& position)
{
     CElement* newElement =new CElement(elementName);
//     #ifdef DEBUG
//        Log ::Debug<<"CConfigurationBase::addAtom "<<std::endl;
//     #endif
     //Check whether new atom lable is valid
     //otherwise the program quit!!!
     this->checkElement(*newElement);

     //add atomic property into the whole atomic-list.
     //and the m_Atom vector only save this pointer of this m_Element
     this->m_Element.push_back(*newElement);

     //add element into bond evaluator for the identification of chemical bonds.
     m_pBondEvaluator->addElement(*newElement);

     CAtom* tempAtom = new CAtom(this,this->m_Atom.size());
     this->m_Atom.push_back(tempAtom);

     //add default charge into vector.
     this->m_pData->partialCharges.push_back(0);

     this->m_pData->atomBonds.push_back(std::vector<CBond*>());

     if(m_CoordinateType==CATAZJUT::DEFINED::Fraction)
        this->Fractioncoordinates()->append(position);
     else if(m_CoordinateType==CATAZJUT::DEFINED::Cartesian)
        this->coordinates()->append(position);

     #ifdef DEBUG
//         Log::Debug<<coordinates()->position(m_pCartesian->size()-1)(0)<<"  ";
//         Log::Debug<<coordinates()->position(m_pCartesian->size()-1)(1)<<"  ";
//         Log::Debug<<coordinates()->position(m_pCartesian->size()-1)(2)<< std::endl;
     #endif
}
 void CConfigurationBase::addAtom(const std::string& elementName, const double& x, const double& y, const double& z)
 {
     Point3 tempPoint;
     tempPoint<<x,y,z;
     addAtom(elementName,tempPoint);
 }
void CConfigurationBase::addAtom(const std::string& elementName, const Point3i& connection,const Point3& coordinate)
{
    CElement* newElement =new CElement(elementName);
    //Check whether new atom lable is valid
    //otherwise the program quit!!!
    this->checkElement(*newElement);
    this->m_Element.push_back(std::move(*newElement));
    //add element into bond evaluator for the identification of chemical bonds.
    this->m_pBondEvaluator->addElement(*newElement);

    CAtom* tempAtom = new CAtom(this,this->m_Atom.size());
    this->m_Atom.push_back(tempAtom);

    //add default charge into vector.
    this->m_pData->partialCharges.push_back(0);

    this->Internalcoordinates()->addInternalCoordinates(connection,coordinate);
}
void CConfigurationBase::addAtom(const std::string& elementName, int atom1, double dist, \
                                  int atom2, double angle, int  atom3, double twist)
{
     Point3 tempCoord;
     tempCoord<<dist,angle,twist;
     util::Point3i tempConnection;
     tempConnection<<atom1,atom2,atom3;
     addAtom(elementName,tempConnection,tempCoord);
}
CAtom* CConfigurationBase::addAtom(const CElement &element)
{
     //Check whether new atom lable is valid
    //otherwise the program quit!!!
     //this->checkElement(element);
     this->m_Element.push_back(element);

     CAtom* tempAtom = new CAtom(this,this->m_Atom.size());
     this->m_Atom.push_back(tempAtom);

     this->m_pData->atomBonds.push_back(std::vector<CBond*>());

     //add default charge into vector.
     this->m_pData->partialCharges.push_back(0);

     if(m_CoordinateType==CATAZJUT::DEFINED::Fraction)
        this->Fractioncoordinates()->append(0,0,0);
     else if(m_CoordinateType==CATAZJUT::DEFINED::Cartesian)
        this->coordinates()->append(0,0,0);
     return tempAtom;
}
CAtom* CConfigurationBase::addAtomCopy(const CAtom *atom)
{
    CAtom *tempAtom = addAtom(atom->element());

    tempAtom->setPartialCharge(atom->PartialCharge());

    tempAtom->SetPosition(atom->position());
    return tempAtom;
}
void CConfigurationBase::removeAtom(CAtom *atom)
{
    if(contains(atom)==false){
        return;
    }
    this->removeBonds(atom->bonds());
    // coupling operation by erase and move.
    m_Atom.erase(std::remove(m_Atom.begin(),m_Atom.end(),atom),m_Atom.end());

    m_Element.erase(m_Element.begin() + atom->index());
    m_pData->atomBonds.erase(m_pData->atomBonds.begin() + atom->index());
    m_pData->partialCharges.erase(m_pData->partialCharges.begin() + atom->index());

    if(atom->index()<m_pData->atomTypes.size())
        m_pData->atomTypes.erase(m_pData->atomTypes.begin() + atom->index());
    if(m_pCartesian)
        m_pCartesian->remove(atom->index());

    for(size_t i=atom->index();i<m_Atom.size();i++)
        m_Atom[i]->m_index--;
    atom->m_pConfiguration=nullptr;
    delete atom;
}
void CConfigurationBase::removeAtoms(const std::vector<CAtom *> &atoms)
{
    BOOST_REVERSE_FOREACH(CAtom *atom, atoms)
    {
        removeAtom(atom);
    }
}
bool CConfigurationBase::contains(const CAtom *atom) const
{
    if(atom->Configuration()==this)
        return true;
    return false;
}
bool CConfigurationBase::contains(const CElement &element) const
{
    return std::find(m_Element.begin(),m_Element.end(),element) != m_Element.end();
}
// sort atom sequence by atomic number.
void CConfigurationBase::sortAtomsViaElements()
{
    std::sort(m_Atom.begin(),m_Atom.end(),[](CAtom* atom_a, CAtom* atom_b)   \
              { return atom_a->AtomNumber() < atom_b->AtomNumber();});

    //reset coordinate library!
    CCartesianCoordinates* tempCoordinate = new CCartesianCoordinates(*m_pCartesian);
    size_t tmp_index=0;
    foreach(CAtom* atom_s, this->atoms()){
       m_pCartesian->setPosition(tmp_index,tempCoordinate->position(atom_s->index()));
       atom_s->setIndex(tmp_index);
       tmp_index++;
    }
    delete tempCoordinate;
    // re-construct bond connectivity of whole configuration
    this->removeBonds(m_pData->bonds);
    this->perceiveBonds();
}


//bond property
/*

*/
void CConfigurationBase::perceiveBonds()
{
    assert(this->m_pBondEvaluator);

    //clear all bonds;
    if(this->bondCount()!=0)
      removeBonds(m_pData->bonds);

    for(size_t i=0;i<m_Atom.size()-1;i++)
        for(size_t j=i+1;j<m_Atom.size();j++)
        {
            if(m_pBondEvaluator->IsBond(m_Atom[i],m_Atom[j]) == true){
               // first generate bond and identify bond order
//               #ifdef DEBUG
//                  Log::Debug<<"*********** 1-CConfigurationBase::perceiveBonds()***********"<< std::endl;
//               #endif // DEBU
                this->setBondOrderType(this->addBond(m_Atom[i],m_Atom[j]));
            }
        }
//    #ifdef DEBUG
//       Log::Debug<<"***********2- CConfigurationBase::perceiveBonds()***********"<< std::endl;
//    #endif // DEBUG
}

CBond* CConfigurationBase::addBond(CAtom *a, CAtom *b)
{
    if(a == b)
        return nullptr;
    if(!contains(a)||!contains(b))
        return nullptr;
    if(a->isBondedTo(b))
        return bond(a,b);
    //create new bond;
    CBond *newBond =  new CBond(this,m_pData->bonds.size());
    m_pData->atomBonds[a->index()].push_back(newBond);
    m_pData->atomBonds[b->index()].push_back(newBond);
    m_pData->bonds.push_back(newBond);
    m_pData->bondAtoms.push_back(std::make_pair(a, b));

    return newBond;
}
CBond* CConfigurationBase::addBond(size_t a, size_t b)
{
    return addBond(m_Atom[a],m_Atom[b]);
}
void CConfigurationBase::removeBond(CBond *bond)
{
    assert(bond->Configuration()==this);

    m_pData->bonds.erase(m_pData->bonds.begin() + bond->index());

    std::vector<CBond*> &bondA = m_pData->atomBonds[bond->atom1()->index()];
    std::vector<CBond*> &bondB = m_pData->atomBonds[bond->atom2()->index()];
    bondA.erase(std::find(bondA.begin(),bondA.end(),bond));
    bondB.erase(std::find(bondB.begin(),bondB.end(),bond));

    m_pData->atomBonds.erase(m_pData->atomBonds.begin()+bond->index());

    m_pData->bondOrders.erase(m_pData->bondOrders.begin()+bond->index());

    for(size_t i=bond->index();m_pData->bonds.size();i++)
        m_pData->bonds[i]->m_index--;

    delete bond;
}
void CConfigurationBase::removeBond(CAtom *a, CAtom *b)
{
    CBond *myBond = this->bond(a,b);
    if( myBond != nullptr )
       removeBond(myBond);
}
void CConfigurationBase::removeBond(size_t a, size_t b)
{
    removeBond(bond(a,b));
}
void CConfigurationBase::removeBonds(const std::vector<CBond *> &bonds)
{
    BOOST_REVERSE_FOREACH(CBond* mybond,bonds){
          removeBond(mybond);
    }
}
CBond* CConfigurationBase::bond(size_t index) const
{
    return m_pData->bonds[index];
}
CBond* CConfigurationBase::bond(const CAtom *a, const CAtom *b) const
{
    return a->bondTo(b);
}
CBond* CConfigurationBase::bond(size_t a, size_t b) const
{
    return bond(m_Atom[a],m_Atom[b]);
}
CConfigurationBase::BondRange CConfigurationBase::bonds() const
{
    return boost::make_iterator_range(m_pData->bonds.begin(),m_pData->bonds.end());
}
size_t CConfigurationBase::bondCount() const
{
    return m_pData->bonds.size();
}
bool CConfigurationBase::contains(const CBond *bond) const
{
    foreach(CBond* bond,m_pData->bonds)
    {
        if(bond->containsBoth(bond->atom1(),bond->atom2())== true )
            return true;
    }
    return false;
}
bool CConfigurationBase::isBondBetween(const CAtom* atom1, const CAtom* atom2)
{
   return m_pBondEvaluator->IsBond(atom1,atom2);
}
void CConfigurationBase::clear()
{
     removeBonds(m_pData->bonds);
     removeAtoms(m_Atom);
}
void CConfigurationBase::setBondOrderType(CBond* mthBond)
{
    if(mthBond->atom1()->maxCoordinationNum()==1 || mthBond->atom1()->maxCoordinationNum()==1)
        this->m_pData->bondOrders.push_back(CBond::Single);
    else if(mthBond->atom1()->element().isTransitionMetal() || mthBond->atom2()->element().isTransitionMetal() )
        this->m_pData->bondOrders.push_back(CBond::Single);
    else if(mthBond->atom1()->element().isPblockElement()&& mthBond->atom2()->element().isPblockElement()){
        double coventL1,coventL2;
        coventL1=mthBond->atom1()->CovalentRadius();
        coventL2=mthBond->atom2()->CovalentRadius();
        if(mthBond->length() > (coventL1+coventL2-0.2))
            this->m_pData->bondOrders.push_back(CBond::Single);
        else if(mthBond->length()<=(coventL1+coventL2-0.2) && (mthBond->length()>(coventL1+coventL2-0.3)))
            this->m_pData->bondOrders.push_back(CBond::Double);
        else
            this->m_pData->bondOrders.push_back(CBond::Triple);
    }else
        this->m_pData->bondOrders.push_back(CBond::Single);
//    #ifdef DEBUG
//       Log::Debug<<"*********** 1-CConfigurationBase::setBondOrderType(CBond* mthBond)***********"<< std::endl;
//    #endif // DEBU

}
void CConfigurationBase::maxCoordinationNum(std::map<std::string,size_t>& max_Coord_Num_Map)
{
      std::map<std::string,size_t>::iterator it;

      foreach(const CAtom* atoms,this->m_Atom){
          it=max_Coord_Num_Map.find(atoms->Symbol());
          if(it==max_Coord_Num_Map.end())
            max_Coord_Num_Map.insert(std::pair<std::string,size_t>
                                     (atoms->Symbol(),atoms->bondCount()));
          else{
            if(it->second < atoms->bondCount())
                it->second=atoms->bondCount();
          }
      }
}


//*****Fragments****************//

void CConfigurationBase::perceiveFragments()
{
     if(this->isEmpty()==true)
        return;
     Bitset totalBit(m_Atom.size());
     //set all bit to be 1;
     totalBit.set();
     //
     size_t pos=0;
     while(true){
        //default, the initial value of all bits is set to 0
        Bitset partBit(m_Atom.size());
        //
        std::vector<const CAtom*> row;
        // this pos  is controlled by totalBit vector.
        // it is the position with bit value of 1
        row.push_back(m_Atom[pos]);
        //
        while(!row.empty()){
           std::vector<const CAtom*> nextRow;
           foreach(const CAtom* atom,row){
               //set the bit of atom-index to be 0
               totalBit.set(atom->index(),false);
               ////set the bit of atom-index to be 1
               partBit.set(atom->index());
               foreach(const CAtom* neighborAtom,atom->neighbors()){
                    if(totalBit[neighborAtom->index()])
                        nextRow.push_back(neighborAtom);
               }
           }
           row = nextRow;
        }
        CFragment* newFragment =  new CFragment(this,partBit);
        m_pData->fragments.push_back(newFragment);

        pos = totalBit.find_next(pos);
        if(pos == Bitset::npos)
            break;
    }
    setFragmentsPerceived(true);
}
bool CConfigurationBase::FromMoietyGetFragmentsAtom(Bitset& mBit,CAtom* atom_label)
{
   std::vector<size_t> indexAtom;
   for(size_t i=0;i<mBit.size();i++)
      if( mBit[i] == 1 )
         indexAtom.push_back(i);

   std::vector<CFragment*> temp_Fragment;

     Bitset totalBit(indexAtom.size());
     totalBit.set();    //set elements of totalBit 1
     size_t pos=0;      // 1st totalBit[0] corresponding with indexAtom[0]->Atom

     while(true){
        Bitset partBit(m_Atom.size());   // all value is set to 0

        std::vector<const CAtom*> row;
        row.push_back(m_Atom[indexAtom[pos]]);

        while(!row.empty()){
            std::vector<const CAtom*> nextRow;
            // analyze bonding atoms of current atoms ( row )
            foreach(const CAtom* atom,row){
                int atomIndex=std::distance(indexAtom.begin(),\
                                  std::find(indexAtom.begin(),indexAtom.end(),atom->index()));
                //set the value to 1 at atom->index() in parBit, tag this atoms to be selected
                partBit.set(atom->index());

                //set the value to 1 at atom->index() in totalBit, tag this atoms to be selected
                //It will make it be not selected by pos pointer.
                totalBit.set(atomIndex,false);
                std::vector<size_t>::iterator it;
                foreach(const CAtom* neighborAtom,atom->neighbors()){
                    it=std::find(indexAtom.begin(),indexAtom.end(),neighborAtom->index());
                    if(it!=indexAtom.end() && totalBit[std::distance(indexAtom.begin(),it)])
                       nextRow.push_back(neighborAtom);
                }
            }
            row = nextRow;
        }

        CFragment* newFragment =  new CFragment(this,partBit);
        temp_Fragment.push_back(newFragment);
        // search next position with 1 value and
        // assign it to pos.
        pos = totalBit.find_next(pos);
        if(pos==Bitset::npos)  // Jude whether this bitset is over.
            break;
     }
      if(temp_Fragment.size()>1){
         for(size_t i=0;i<temp_Fragment.size();i++)
             if(temp_Fragment[i]->contains(atom_label)){
                 mBit=temp_Fragment[i]->bitSet();
                 break;
             }
         return true;
      }else
         return false;

}
bool CConfigurationBase::isFragments(Bitset& mBit)
{
   std::vector<size_t> res;
   for(size_t i=0;i<mBit.size();i++)
      if( mBit[i] == 1 )
         res.push_back(i);
   return isFragments(res);
}
bool CConfigurationBase::isFragments( std::vector<size_t>& indexAtom)
{
     if(this->isEmpty()==true)
        return false;
     Bitset totalBit(indexAtom.size());
     totalBit.set();    //set elements of totalBit 1
     size_t pos=0;      // 1st totalBit[0] corresponding with indexAtom[0]->Atom
     size_t fragmentNum=0;
     while(true){
        Bitset partBit(m_Atom.size());   // all value is set to 0

        std::vector<const CAtom*> row;
        row.push_back(m_Atom[indexAtom[pos]]);

        while(!row.empty()){
            std::vector<const CAtom*> nextRow;
            // analyze bonding atoms of current atoms ( row )
            foreach(const CAtom* atom,row){
                int atomIndex=std::distance(indexAtom.begin(),\
                                  std::find(indexAtom.begin(),indexAtom.end(),atom->index()));
                //set the value to 1 at atom->index() in parBit, tag this atoms to be selected
                partBit.set(atom->index());

                //set the value to 1 at atom->index() in totalBit, tag this atoms to be selected
                //It will make it be not selected by pos pointer.
                totalBit.set(atomIndex,false);
                std::vector<size_t>::iterator it;
                foreach(CAtom* neighborAtom,atom->neighbors()){
                    it=std::find(indexAtom.begin(),indexAtom.end(),neighborAtom->index());
                    if(it!=indexAtom.end() && totalBit[std::distance(indexAtom.begin(),it)])
                       nextRow.push_back(neighborAtom);
                }
            }
            row = nextRow;
        }
        fragmentNum = fragmentNum+1;

        // search next position with 1 value and
        // assign it to pos.
        pos = totalBit.find_next(pos);
        if(pos==Bitset::npos)  // Jude whether this bitset is over.
            break;
     }
      if(fragmentNum>1)
        return true;
      else
        return false;
}
CFragment* CConfigurationBase::fragmentForAtom(const CAtom *atom)
{
    foreach(CFragment* fragment,m_pData->fragments)
    {
         if(fragment->contains(atom)==true)
            return fragment;
    }
    return nullptr;
}
CConfigurationBase::FragmentRange CConfigurationBase::fragments()
{
    if(this->fragmentsPerceived()==false){
       perceiveFragments();
       this->setFragmentsPerceived(true);
    }
    return boost::make_iterator_range(m_pData->fragments.begin(),m_pData->fragments.end());
}
CFragment* CConfigurationBase::fragment(size_t index)
{
    if(index>=m_pData->fragments.size())
        return nullptr;

    return m_pData->fragments[index];
}
void CConfigurationBase::setFragmentsPerceived(bool perceived)
{
//    if( perceived==m_pData->fragmentsPerceived && m_pData->fragmentsPerceived==true)
//        return;
    if(!perceived){
        foreach(const CFragment* fragment, m_pData->fragments){
            delete fragment;
        }
        m_pData->fragments.clear();
    }
    m_pData->fragmentsPerceived = perceived;
}
size_t CConfigurationBase::fragmentNum()
{
    if(fragmentsPerceived()==false)
        perceiveFragments();
    return m_pData->fragments.size();
}

bool CConfigurationBase::fragmentsPerceived()const
{
    return m_pData->fragmentsPerceived;
}
double CConfigurationBase::distance(const CAtom *a, const CAtom *b)
{
    return coordinates()->distance(a->index(),b->index());
}
double CConfigurationBase::bondAngle(const CAtom *a, const CAtom *b, const CAtom *c)
{
   return this->coordinates()->angle(a->index(),b->index(),c->index());
}
double CConfigurationBase::torsionAngle(const CAtom *a, const CAtom *b, const CAtom *c, const CAtom *d)
{
   return this->coordinates()->torsionAngle(a->index(),b->index(),c->index(),d->index());
}
double CConfigurationBase::wilsonAngle(const CAtom *a, const CAtom *b, const CAtom *c, const CAtom *d)
{
   return this->coordinates()->wilsonAngle(a->index(),b->index(),c->index(),d->index());
}
void CConfigurationBase::setCenter(const Point3 &position)
{
     const Vector3 &vector1 = position - this->center();

    foreach(CAtom *atom, m_Atom){
        atom->SetPosition(atom->position() + vector1);
    }
}
Point3 CConfigurationBase::center()
{
    if(this->m_pCartesian==nullptr){
        return Point3(0,0,0);
    }
    return this->m_pCartesian->center();
}
void CConfigurationBase::setCenter(double x, double y, double z)
{
    setCenter(Point3(x,y,z));
}

void CConfigurationBase::setCoordinateType(const CATAZJUT::DEFINED::CoordinateType currentType)
{
    this->m_CoordinateType = currentType;
}


CATAZJUT::DEFINED::CoordinateType CConfigurationBase::coordinateType()const
{
    return this->m_CoordinateType;
}


void CConfigurationBase::setDimensionalType(const CATAZJUT::DEFINED::DimensionalType currentType)
{
    this->m_DimensionalType = currentType;
}


CATAZJUT::DEFINED::DimensionalType CConfigurationBase::dimensionalType()const
{
    return this->m_DimensionalType;
}

Bitset CConfigurationBase::constraintBit()const
{
    return this->constraintBits;
}
void CConfigurationBase::setConstraintBit(const Bitset& othrBit)
{
    this->constraintBits = othrBit;
}
void CConfigurationBase::moveAtom(CAtom* atom, Point3 vect)
{
    atom->SetPosition(atom->position()+vect);
}
void CConfigurationBase::setTolerancefactor(std::pair<double,double> &mht)
{
    this->m_pBondEvaluator->setTolerancefactor(mht);
}
void CConfigurationBase::setExcludeBond(std::vector<std::pair<std::string*,std::string*>>&mht)
{
    this->m_pBondEvaluator->setExcludeBond(mht);
}
CALCZJUT::CParameter *CConfigurationBase::sysParameter()
{
    return this->m_pParameter;
}
CUnitCell* CConfigurationBase::unitcell()
{
    if(this->m_pUnitCell==nullptr)
        this->m_pUnitCell=new CUnitCell();

    return this->m_pUnitCell;
}

std::string CConfigurationBase::SymmetrySymbol()
{
    if(m_pParameter->simulationMode==CALCZJUT::CParameter::CLUSTER ||
       m_pParameter->simulationMode == CALCZJUT::CParameter::MOL_CLUSTER ){
        CMolecularSymmetry* molsym = new CMolecularSymmetry(this);
        std::string res("Molecular symmetry: point group ");
        res=res + molsym->GetPointGroup();
        delete molsym;
        return res;
    }else{
        CSpaceGroup* spaceGroup = new CSpaceGroup(this);
        std::string res("Crystallography symmetry: space group ");
        res=res + spaceGroup->GetSpaceGroup();
        delete spaceGroup;
        return res;
    }
}

} //namespace CATAZJUT
