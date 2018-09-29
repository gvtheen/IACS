#include "CClusterGeneVar.h"
#include "../Util/log.hpp"
#include "../GACatalyst.h"
using util::Log;

namespace CALCZJUT{


CClusterGeneVar::CClusterGeneVar()
{
    //ctor
}
CClusterGeneVar::CClusterGeneVar(VAR_TYPE mth, int i, Point3i m)
:m_type(mth)
{
    this->setIndex(i,m);
}
CClusterGeneVar::~CClusterGeneVar()
{
    this->m_index.clear();
}
CClusterGeneVar::VAR_TYPE CClusterGeneVar::type()
{
    return this->m_type;
}
void CClusterGeneVar::setType(const CClusterGeneVar::VAR_TYPE& mth)
{
    this->m_type = mth;
}
std::vector<size_t> CClusterGeneVar::index()
{
   return this->m_index;
}
void CClusterGeneVar::setIndex(std::vector<size_t>& mht)
{
    this->m_index.assign(mht.begin(),mht.end());
}
void CClusterGeneVar::setIndex(int i,int j,int k,int h)
{
   if( i<0 || j <0){
      Log::Error<<" i and j are error!! setIndex_CClusterGeneVar\n";
      boost::throw_exception(std::runtime_error("i and j are error!! setIndex_CClusterGeneVar\n"));
   }else{
      this->m_index.push_back(i);
      this->m_index.push_back(j);
   }
   if( k > -1 ){
      this->m_index.push_back(k);
      if ( h > -1 )
        this->m_index.push_back(h);
   }

}
void CClusterGeneVar::setIndex(size_t i,Point3i m)
{
    size_t j,k,h;
    j=m[0];
    k=m[1];
    h=m[2];

    this->setIndex(i,j,k,h);
}
bool CClusterGeneVar::operator==(CClusterGeneVar& mth)
{
    if( this->m_type != mth.type())
        return false;
    if( this->m_index.size() != mth.index().size())
        return false;
    if(this->m_index == mth.index())
        return true;
    else
        return false;
}
void CClusterGeneVar::operator = (CClusterGeneVar& mth)
{
    this->m_type = mth.type();
    this->m_index = mth.index();
}





}
