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
#include <string>
#include<iostream>
#include<fstream>
#include "unistd.h"
#include "stdlib.h"
#include "CBondTolerance.h"
#include "CConfigurationBase.h"
#include "CAtom.h"
#include "CElement.h"
#include "CCoordinateSet.h"
#include "CBondPrivate.h"
#include "CBond.h"
#include "../Util/log.hpp"
#include "../GaZJUT/GaUtilityFunction.h"

using util::Log;
namespace CATAZJUT{

CBondTolerance::~CBondTolerance()
{
    for(size_t i=0; i<this->m_pBondType.size();i++)
        delete m_pBondType[i];
}
/*
default values are derived from MS software
lower_tolerance_factor:   0.6
upper_tolerance_factor:   1.15
*/
CBondTolerance::CBondTolerance(CConfigurationBase* curConfig)
   :m_pConfiguration(curConfig),lower_tolerance_factor(0.6), upper_tolerance_factor(1.15)
{
    //ctor
}
CBondTolerance::CBondTolerance(CConfigurationBase* curConfig,double minF,double maxF)
   :m_pConfiguration(curConfig),lower_tolerance_factor(minF), upper_tolerance_factor(maxF)
{

}
CBondTolerance::CBondTolerance(CConfigurationBase* curConfig,std::vector<std::pair<std::string*,std::string*>> excludebond)
  :m_pConfiguration(curConfig),m_pExcludeBond(excludebond)
{
   lower_tolerance_factor = 0.6;
   upper_tolerance_factor = 1.15;
}

CBondTolerance::CBondTolerance(CBondTolerance& other)
{

}
void CBondTolerance::addElement(const CElement& othr)
{
    CElement* newElement = new CElement(othr);
    std::vector<CElement>::iterator it = std::find(m_pElement.begin(),m_pElement.end(),*newElement);
    if(it==m_pElement.end())
    {
        for(it=m_pElement.begin();it!=m_pElement.end();it++)
           AddBondType(*it,*newElement);
        m_pElement.push_back(*newElement);
    }else
        delete newElement;
}
void CBondTolerance::removeElement(CElement& othr)
{
    m_pElement.erase(std::find(m_pElement.begin(),m_pElement.end(),othr));

    std::vector<CBondPrivate*>::iterator it;
    for(it=m_pBondType.begin(); it!=m_pBondType.end();it++)
        if((*it)->m_1stAtom == othr.symbol() || (*it)->m_2ndAtom == othr.symbol())
            m_pBondType.erase(it);
}

int CBondTolerance::IsExistBond(const std::string& _1stAtom,const std::string& _2ndAtom)
{
    int res= -1;
//    #ifdef DEBUG
//         Log::Debug<<"CBondTolerance::IsExistBond"<<std::endl;
//    #endif
    for(size_t i =0;i<m_pBondType.size();i++)
        if((m_pBondType[i]->m_1stAtom == _1stAtom && m_pBondType[i]->m_2ndAtom == _2ndAtom)|| \
           (m_pBondType[i]->m_1stAtom == _2ndAtom && m_pBondType[i]->m_2ndAtom == _1stAtom) )
           {
               res=i;
               break;
           }
    return res;
}
void CBondTolerance::AddBondType(const std::string& _1stAtom,const std::string& _2ndAtom,const double& mind,const double& maxd)
{
    int index = IsExistBond(_1stAtom,_2ndAtom);
    bool exclude_bond = isExcludeBond(_1stAtom,_2ndAtom);
    if( index == -1 && !exclude_bond )
    {
        CBondPrivate *temp =new CBondPrivate();
        temp->m_1stAtom = _1stAtom;
        temp->m_2ndAtom = _2ndAtom;
        temp->m_MinBondLength = mind;
        temp->m_MaxBondLength = maxd;
        m_pBondType.push_back(temp);
    }else if (index >=0 && !exclude_bond ) {
        m_pBondType[index]->m_MinBondLength=mind;
        m_pBondType[index]->m_MaxBondLength=maxd;
    }
}
void CBondTolerance::AddBondType(const std::string& _1stAtom,const std::string& _2ndAtom)
{

//    #ifdef DEBUG
//         Log::Debug<<"CBondTolerance::AddBondType"<<std::endl;
//    #endif
    int index = IsExistBond(_1stAtom,_2ndAtom);
//    #ifdef DEBUG
//         Log::Debug<<"2-CBondTolerance::AddBondType"<<std::endl;
//    #endif
    if( index==-1 && isExcludeBond(_1stAtom,_2ndAtom)==false )
    {
        CBondPrivate *temp =new CBondPrivate();
        temp->m_1stAtom = _1stAtom;
        temp->m_2ndAtom = _2ndAtom;
        double vdw_12= (new CElement(_1stAtom))->vanDerWaalsRadius() + \
                       (new CElement(_2ndAtom))->vanDerWaalsRadius();
        temp->m_MinBondLength = vdw_12*this->lower_tolerance_factor;
        temp->m_MaxBondLength = vdw_12*this->upper_tolerance_factor;
        m_pBondType.push_back(temp);
    }
//    #ifdef DEBUG
//         Log::Debug<<"end-CBondTolerance::AddBondType"<<std::endl;
//    #endif
}
bool CBondTolerance::IsBond(const std::string& _1stAtom,const std::string& _2ndAtom,const double& dis)
{
//    #ifdef DEBUG
//         Log::Debug<<"CBondTolerance::IsBond"<<std::endl;
//    #endif
    int index = IsExistBond(_1stAtom,_2ndAtom);
    //std::cout<<"CBondTolerance::IsBond: index: "<<index<<std::endl;

    if(index>0 && index < (int)m_pBondType.size())
       return ( m_pBondType[index]->m_MinBondLength <= dis && \
                m_pBondType[index]->m_MaxBondLength >= dis );

    else if(isExcludeBond(_1stAtom,_2ndAtom)==false){
//        #ifdef DEBUG
//         Log::Debug<<"2-CBondTolerance::IsBond"<<std::endl;
//         std::cout<<"m_pBondType.size():"<<m_pBondType.size()<<std::endl;
//        #endif
        AddBondType(_1stAtom,_2ndAtom);
//        #ifdef DEBUG
//         Log::Debug<<"3-CBondTolerance::IsBond"<<std::endl;
//        #endif
        int index_2 = m_pBondType.size()-1;
        return ( m_pBondType.at(index_2 )->m_MinBondLength <= dis && \
                 m_pBondType.at(index_2 )->m_MaxBondLength >= dis );
    }else
        return false;
}
bool CBondTolerance::IsBond(const size_t& _1st,const size_t& _2nd,const double& dis)
{
    return this->IsBond(m_pConfiguration->m_Atom[_1st]->Symbol(),\
                        m_pConfiguration->m_Atom[_2nd]->Symbol(),dis);
}
bool CBondTolerance::IsBond(const CAtom* atom1,const CAtom* atom2)
{
   double dis = atom1->distance(atom2);
   return IsBond(atom1->element().symbol(),atom2->element().symbol(),dis);
}
void CBondTolerance::RemoveBondType(const std::string& _1stAtom,const std::string& _2ndAtom)
{
   size_t index = IsExistBond(_1stAtom,_2ndAtom);
   m_pBondType.erase(m_pBondType.begin()+index);
}
void CBondTolerance::AddBondType(const CElement& Elem_1,const CElement& Elem_2,const double& mind,const double& maxd)
{
     AddBondType(Elem_1.symbol(),Elem_2.symbol(),mind,maxd);
}
void CBondTolerance::AddBondType(const CElement& Elem_1,const CElement& Elem_2)
{
     AddBondType(Elem_1.symbol(),Elem_2.symbol());
}
void CBondTolerance::OutputTofile(std::string& outfile)
{
    std::ofstream out(outfile,std::ios::app);
    out.setf(std::ios::fixed, std::ios::floatfield);
    out.precision(10);

    if (out.is_open()){
        out<<"**Bond tolerance information**"<<std::endl;
        for(size_t i=0;i<m_pBondType.size();i++)
        {
            out<<"||"<<m_pBondType[i]->m_1stAtom<<"-"<<m_pBondType[i]->m_2ndAtom<<":  ";
            out<<m_pBondType[i]->m_MinBondLength<<"   "<<m_pBondType[i]->m_MaxBondLength;
            out<<std::endl;
        }
        out<<"**End bond tolerance information**"<<std::endl;
    }
    out.close();
}
bool CBondTolerance::isExcludeBond(const std::string& e1,const std::string& e2)
{
    if(m_pExcludeBond.size()==0)
        return false;

    for(size_t i=0;i<m_pExcludeBond.size();i++)
       if( (*(m_pExcludeBond[i].first)==e1 && *(m_pExcludeBond[i].second)==e2 ) || \
           (*(m_pExcludeBond[i].first)==e2 && *(m_pExcludeBond[i].second)==e1 ) )
              return true;
//    #ifdef DEBUG
//         Log::Debug<<"CBondTolerance::isExcludeBond"<<std::endl;
//    #endif
    return false;
}
bool CBondTolerance::isExcludeBond(CElement& e1,CElement& e2)
{
    return isExcludeBond(e1.symbol(),e2.symbol());
}
void CBondTolerance::setTolerancefactor(std::pair<double,double> &mht)
{
    if( mht.first != 0 && mht.second != 0 )
       setTolerancefactor(mht.first,mht.second);
}
void CBondTolerance::setTolerancefactor(const double minV,const double maxV)
{
//    #ifdef DEBUG
//      Log::Debug<<"*********** CBondTolerance::setTolerancefactor***********"<< std::endl;
//    #endif
    this->lower_tolerance_factor=minV;
    this->upper_tolerance_factor=maxV;
    double vdw_12;
    for(size_t i=0;i<this->m_pBondType.size();i++)
    {
        vdw_12=(new CElement(m_pBondType[i]->m_1stAtom))->vanDerWaalsRadius() +
               (new CElement(m_pBondType[i]->m_2ndAtom))->vanDerWaalsRadius();

        m_pBondType[i]->m_MinBondLength = vdw_12*this->lower_tolerance_factor;
        m_pBondType[i]->m_MaxBondLength = vdw_12*this->upper_tolerance_factor;
    }
//    #ifdef DEBUG
//      Log::Debug<<"2-*********** CBondTolerance::setTolerancefactor***********"<< std::endl;
//    #endif
}
void CBondTolerance::setExcludeBond(std::vector<std::pair<std::string*,std::string*>>& tmp)
{
    this->m_pExcludeBond.insert(m_pExcludeBond.end(),tmp.begin(),tmp.end());
//    #ifdef DEBUG
//      Log::Debug<<"2-*********** CBondTolerance::setExcludeBond***********"<< std::endl;
//    #endif
}





}
