#include <boost/scoped_ptr.hpp>
#include "CCoordinateSet.h"
#include "CCartesianCoordinates.h"
#include "CInternalCoordinates.h"
#include "CFractionCoordinates.h"
namespace CATAZJUT{

CCoordinateSet::CCoordinateSet()
{
    m_type = None;
    m_cartesianCordinates = nullptr;
}
CCoordinateSet::CCoordinateSet(CCartesianCoordinates *coordinates)
{
    this->m_type=Cartesian;
    this->m_cartesianCordinates = coordinates;
}
CCoordinateSet::CCoordinateSet(CInternalCoordinates  *coordinates)
{
    this->m_type=Internal;
    this->m_internalCoordinates = coordinates;
}
CCoordinateSet::CCoordinateSet(CFractionCoordinates  *coordinates)
{
    this->m_type=Fraction;
    this->m_fractionCoordinates = coordinates;
}
CCoordinateSet::Type CCoordinateSet::type() const
{
    return this->m_type;
}
size_t CCoordinateSet::size() const
{
    switch (m_type)
    {
        case Cartesian:
             return m_cartesianCordinates->size();
             break;
        case Internal:
             return m_internalCoordinates->size();
             break;
        case Fraction:
             return m_fractionCoordinates->size();
             break;
        default:
             return 0;
             break;
    }
}
bool CCoordinateSet::isEmpty() const
{
    return this->size()==0;
}
void CCoordinateSet::setCoordinates(CCartesianCoordinates *coordinates)
{
    this->clear();
    this->m_cartesianCordinates=coordinates;
    this->m_type=Cartesian;
}
void CCoordinateSet::setCoordinates(CInternalCoordinates *coordinates)
{
    this->clear();
    this->m_internalCoordinates=coordinates;
    this->m_type=Internal;
}
void CCoordinateSet::setCoordinates(CFractionCoordinates *coordinates)
{
    this->clear();
    this->m_fractionCoordinates=coordinates;
    this->m_type=Fraction;
}
void CCoordinateSet::clear()
{
    switch (m_type)
    {
        case Cartesian:
             delete m_cartesianCordinates;
             m_cartesianCordinates=nullptr;
             break;
        case Internal:
             delete m_internalCoordinates;
             m_internalCoordinates=nullptr;
             break;
        case Fraction:
             delete m_fractionCoordinates;
             m_fractionCoordinates=nullptr;
             break;
        default:
             break;
    }
}
CCartesianCoordinates* CCoordinateSet::cartesianCoordinates() //const
{
    if(this->m_type==Cartesian)
        return this->m_cartesianCordinates;
    return nullptr;
}
CInternalCoordinates*  CCoordinateSet::internalCoordinates()// const
{
    if(this->m_type==Internal)
        return this->m_internalCoordinates;
    return nullptr;
}
CFractionCoordinates*  CCoordinateSet::fractionCoordinates()// const
{
    if(this->m_type==Fraction)
        return this->m_fractionCoordinates;
    return nullptr;
}
Point3 CCoordinateSet::position(size_t index) const
{
    if(m_type == Cartesian)
        return m_cartesianCordinates->position(index);
    else if(m_type == Internal){
        boost::scoped_ptr<CCartesianCoordinates> coord1(m_internalCoordinates->toCartesianCoordinates());
        return coord1->position(index);
    }else if(m_type == Fraction){
        boost::scoped_ptr<CCartesianCoordinates> coord2(m_fractionCoordinates->toCartesianCoordinates());
        return coord2->position(index);
    }else
        return Point3();
}
CCoordinateSet& CCoordinateSet::operator=(CCoordinateSet &other)
{
    if(this!= &other)
    {
        if(other.type()==CCoordinateSet::None)
           clear();
        else if(other.type()==Cartesian)
            setCoordinates(new CCartesianCoordinates( *other.cartesianCoordinates()));
        else if(other.type()==Internal)
            setCoordinates(new CInternalCoordinates( *other.internalCoordinates()));
        else if(other.type()==Fraction)
            setCoordinates(new CFractionCoordinates( *other.fractionCoordinates()));
    }
    return *this;
}
CCoordinateSet::~CCoordinateSet()
{
    this->clear();
}



}
