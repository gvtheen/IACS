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
#include <boost/scoped_ptr.hpp>
#include "CCoordinateSet.h"
#include "CCartesianCoordinates.h"
#include "CInternalCoordinates.h"
#include "CFractionCoordinates.h"
#include "../Util/log.hpp"

using util::Log;
namespace CATAZJUT{

CCoordinateSet::CCoordinateSet()
{
    m_type = CATAZJUT::DEFINED::None;
    m_cartesianCordinates = nullptr;
}
CCoordinateSet::CCoordinateSet(CCartesianCoordinates *coordinates)
{
    this->m_type=CATAZJUT::DEFINED::Cartesian;
    this->m_cartesianCordinates = coordinates;
}
CCoordinateSet::CCoordinateSet(CInternalCoordinates  *coordinates)
{
    this->m_type=CATAZJUT::DEFINED::Internal;
    this->m_internalCoordinates = coordinates;
}
CCoordinateSet::CCoordinateSet(CFractionCoordinates  *coordinates)
{
//    #ifdef DEBUG
//        Log ::Debug<<"CCoordinateSet(CFractionCoordinates  *coordinates) "<<std::endl;
//    #endif
    this->m_type=CATAZJUT::DEFINED::Fraction;
    this->m_fractionCoordinates = coordinates;
}
CATAZJUT::DEFINED::CoordinateType CCoordinateSet::type() const
{
    return this->m_type;
}
size_t CCoordinateSet::size() const
{
    switch (m_type)
    {
        case CATAZJUT::DEFINED::Cartesian:
             return m_cartesianCordinates->size();
             break;
        case CATAZJUT::DEFINED::Internal:
             return m_internalCoordinates->size();
             break;
        case CATAZJUT::DEFINED::Fraction:
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
    this->m_type=CATAZJUT::DEFINED::Cartesian;
}
void CCoordinateSet::setCoordinates(CInternalCoordinates *coordinates)
{
    this->clear();
    this->m_internalCoordinates=coordinates;
    this->m_type=CATAZJUT::DEFINED::Internal;
}
void CCoordinateSet::setCoordinates(CFractionCoordinates *coordinates)
{
    this->clear();
    this->m_fractionCoordinates=coordinates;
    this->m_type=CATAZJUT::DEFINED::Fraction;
}
void CCoordinateSet::clear()
{
    switch (m_type)
    {
        case CATAZJUT::DEFINED::Cartesian:
             delete m_cartesianCordinates;
             m_cartesianCordinates=nullptr;
             break;
        case CATAZJUT::DEFINED::Internal:
             delete m_internalCoordinates;
             m_internalCoordinates=nullptr;
             break;
        case CATAZJUT::DEFINED::Fraction:
             delete m_fractionCoordinates;
             m_fractionCoordinates=nullptr;
             break;
        default:
             break;
    }
}
CCartesianCoordinates* CCoordinateSet::cartesianCoordinates() //const
{
    if(this->m_type==CATAZJUT::DEFINED::Cartesian)
        return this->m_cartesianCordinates;
    return nullptr;
}
CInternalCoordinates*  CCoordinateSet::internalCoordinates()// const
{
    if(this->m_type==CATAZJUT::DEFINED::Internal)
        return this->m_internalCoordinates;
    return nullptr;
}
CFractionCoordinates*  CCoordinateSet::fractionCoordinates()// const
{
    if(this->m_type==CATAZJUT::DEFINED::Fraction)
        return this->m_fractionCoordinates;
    return nullptr;
}
Point3 CCoordinateSet::position(size_t index) const
{
    if(m_type == CATAZJUT::DEFINED::Cartesian)
        return m_cartesianCordinates->position(index);
    else if(m_type == CATAZJUT::DEFINED::Internal){
        boost::scoped_ptr<CCartesianCoordinates> coord1(m_internalCoordinates->toCartesianCoordinates());
        return coord1->position(index);
    }else if(m_type == CATAZJUT::DEFINED::Fraction){
        boost::scoped_ptr<CCartesianCoordinates> coord2(m_fractionCoordinates->toCartesianCoordinates());
        return coord2->position(index);
    }else
        return Point3();
}
CCoordinateSet& CCoordinateSet::operator=(CCoordinateSet &other)
{
    if(this!= &other)
    {
        if(other.type()==CATAZJUT::DEFINED::None)
           clear();
        else if(other.type()==CATAZJUT::DEFINED::Cartesian)
            setCoordinates(new CCartesianCoordinates( *other.cartesianCoordinates()));
        else if(other.type()==CATAZJUT::DEFINED::Internal)
            setCoordinates(new CInternalCoordinates( *other.internalCoordinates()));
        else if(other.type()==CATAZJUT::DEFINED::Fraction)
            setCoordinates(new CFractionCoordinates( *other.fractionCoordinates()));
    }
    return *this;
}
CCoordinateSet::~CCoordinateSet()
{
    this->clear();
}



}
