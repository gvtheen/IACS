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
#include<cmath>
#include "CFractionCoordinates.h"
#include "CCartesianCoordinates.h"
#include "CPeriodicFramework.h"
#include "Geometry.h"
#include "CUnitCell.h"
#include "../Util/log.hpp"

using util::Log;
namespace CATAZJUT{

CFractionCoordinates::CFractionCoordinates()
  :m_pPeriodicFramework(nullptr)
{
}
//
CFractionCoordinates::CFractionCoordinates(CPeriodicFramework* mPeriodicFramework)
  :m_pPeriodicFramework(mPeriodicFramework)
{
}
//
CFractionCoordinates::CFractionCoordinates(CPeriodicFramework* mconf,size_t size_m)
:m_pPeriodicFramework(mconf)
{
   this->m_coordinates.insert(this->m_coordinates.begin(),size_m,Point3().setZero());
}
//
CFractionCoordinates::CFractionCoordinates(const CFractionCoordinates& mht)
{
    this->m_pPeriodicFramework = mht.m_pPeriodicFramework;
    this->m_coordinates.assign(mht.m_coordinates.begin(),mht.m_coordinates.end());
}

CFractionCoordinates::~CFractionCoordinates()
{
    clear();
}

void CFractionCoordinates::resize(size_t size)
{
    m_coordinates.resize(size);
}
size_t CFractionCoordinates::size() const
{
    return m_coordinates.size();
}
bool CFractionCoordinates::isEmpty() const
{
    return m_coordinates.empty();
}
// setting and return the periodic framework
void CFractionCoordinates::setPeriodicFramework(CPeriodicFramework* tempPeridFramekwork)
{
    m_pPeriodicFramework = tempPeridFramekwork;
}
CPeriodicFramework* CFractionCoordinates::PeriodicFramework() const
{
    return this->m_pPeriodicFramework;
}

//
void CFractionCoordinates::setPosition(size_t index, const Point3 &position)
{
    m_coordinates[index]=position;
}
void CFractionCoordinates::setPosition(size_t index, double x, double y, double z)
{
    Point3 p3;
    p3<<x,y,z;
    setPosition(index,p3);
}
Point3 CFractionCoordinates::position(size_t index) const
{
    return m_coordinates[index];
}
void CFractionCoordinates::setValue(size_t row, size_t column, double value)
{
    assert(row<size());
    assert(column<3);
    m_coordinates[row][column]=value;
}
double CFractionCoordinates::value(size_t row, size_t column) const
{
    assert(row<size());
    assert(column<3);
    return m_coordinates[row][column];
}
void CFractionCoordinates::append(const Point3 &position)
{
    m_coordinates.push_back(position);
}
void CFractionCoordinates::append(double x, double y, double z)
{
    Point3 p1;
    p1<<x,y,z;
    append(p1);
}
void CFractionCoordinates::insert(size_t index, const Point3 &position)
{
    setPosition(index,position);
}
void CFractionCoordinates::insert(size_t index, double x, double y, double z)
{
    setPosition(index,x,y,z);
}
void CFractionCoordinates::remove(size_t index)
{
    m_coordinates.erase(m_coordinates.begin() + index);
}
void CFractionCoordinates::clear()
{
    m_coordinates.clear();
}

// if the points is out of lattice, move it to the lattice.
double CFractionCoordinates::movePointIntoCell(double temp)
{
    if( temp >0 && temp < 1)
        return temp;
    else if(std::fabs(temp)<1 && temp<0)
        return -1*temp;
    else
        return std::fabs(temp) - std::floor(std::fabs(temp));

}
void CFractionCoordinates::adjustCoordinateIntoCell()
{
    for(size_t i=0;i<m_coordinates.size();i++)
        for( size_t j=0;j<3;j++)
            m_coordinates[i][j] = movePointIntoCell( m_coordinates[i][j] );
}

// geometrical property
double CFractionCoordinates::distance(size_t i, size_t j)
{
    Point3 p1,p2;
    p1=toCartesian(m_coordinates[i]);
    p2=toCartesian(m_coordinates[j]);
    return CATAZJUT::Geometry::distance(p1,p2);
}
double CFractionCoordinates::angle(size_t i, size_t j, size_t k)
{
    Point3 p1,p2,p3;
    p1=toCartesian(m_coordinates[i]);
    p2=toCartesian(m_coordinates[j]);
    p3=toCartesian(m_coordinates[k]);
    return CATAZJUT::Geometry::angle(p1,p2,p3);
}
double CFractionCoordinates::angleRadians(size_t i, size_t j, size_t k)
{
    Point3 p1,p2,p3;
    p1=toCartesian(m_coordinates[i]);
    p2=toCartesian(m_coordinates[j]);
    p3=toCartesian(m_coordinates[k]);
    return CATAZJUT::Geometry::angleRadians(p1,p2,p3);
}
double CFractionCoordinates::torsionAngle(size_t i, size_t j, size_t k, size_t l)
{
    Point3 p1,p2,p3,p4;
    p1=toCartesian(m_coordinates[i]);
    p2=toCartesian(m_coordinates[j]);
    p3=toCartesian(m_coordinates[k]);
    p4=toCartesian(m_coordinates[l]);
    return CATAZJUT::Geometry::torsionAngle(p1,p2,p3,p4);
}
double CFractionCoordinates::torsionAngleRadians(size_t i, size_t j, size_t k, size_t l)
{
    Point3 p1,p2,p3,p4;
    p1=toCartesian(m_coordinates[i]);
    p2=toCartesian(m_coordinates[j]);
    p3=toCartesian(m_coordinates[k]);
    p4=toCartesian(m_coordinates[l]);
    return CATAZJUT::Geometry::torsionAngleRadians(p1,p2,p3,p4);
}
double CFractionCoordinates::wilsonAngle(size_t i, size_t j, size_t k, size_t l)
{
    Point3 p1,p2,p3,p4;
    p1=toCartesian(m_coordinates[i]);
    p2=toCartesian(m_coordinates[j]);
    p3=toCartesian(m_coordinates[k]);
    p4=toCartesian(m_coordinates[l]);
    return CATAZJUT::Geometry::wilsonAngle(p1,p2,p3,p4);
}
double CFractionCoordinates::wilsonAngleRadians(size_t i, size_t j, size_t k, size_t l)
{
    Point3 p1,p2,p3,p4;
    p1=toCartesian(m_coordinates[i]);
    p2=toCartesian(m_coordinates[j]);
    p3=toCartesian(m_coordinates[k]);
    p4=toCartesian(m_coordinates[l]);
    return CATAZJUT::Geometry::wilsonAngleRadians(p1,p2,p3,p4);
}

Point3 CFractionCoordinates::center()
{
    Point3 centerP = toCartesianCoordinates()->center();
    return toFractionFromCart(centerP);
}
Point3 CFractionCoordinates::weightedCenter(const std::vector<double> &weights)
{
   Point3 centerP = toCartesianCoordinates()->weightedCenter(weights);
   return toFractionFromCart(centerP);
}
void CFractionCoordinates::moveBy(const Vector3 &m_vector)
{
    CCartesianCoordinates* res = toCartesianCoordinates();
    res->moveBy(toCartesian(m_vector));
    for(size_t i=0;i<m_coordinates.size();i++)
        this->insert(i,toFractionFromCart(res->position(i)));
}
void CFractionCoordinates::moveBy(double x, double y, double z)
{
    Vector3 vec1;
     vec1<<x,y,z;
    this->moveBy(vec1);
}

void CFractionCoordinates::rotate(const Vector3 &axis, double angle)
{

    CCartesianCoordinates* res = this->toCartesianCoordinates();
    res->rotate(toCartesian(axis),angle);
    for(size_t i=0;i<m_coordinates.size();i++)
        insert(i,toFractionFromCart(res->position(i)));
}

//conversion
Point3 CFractionCoordinates::toCartesian(const size_t& index) const
{
    return toCartesian(m_coordinates[index]);
}
Point3 CFractionCoordinates::toCartesian(const Point3& position) const
{
    assert(m_pPeriodicFramework);

    //double scalingFactor = m_pPeriodicFramework->m_pUnitCell->scalingFactor();
    return (m_pPeriodicFramework->unitcell()->MatrixOfBravaisLattice())*position;
}
Point3 CFractionCoordinates::toCartesian(const double x, const double y, const double z) const
{
    Point3 p1;
    p1<<x,y,z;
    return toCartesian(p1);
}
Point3 CFractionCoordinates::toFractionFromCart(const Point3& position)
{
   assert(m_pPeriodicFramework);

   Point3 tmpP;
   tmpP<<position[0],
         position[1],
         position[2];
   return (m_pPeriodicFramework->unitcell()->MatrixOfBravaisLattice().inverse())*tmpP;
}
Point3 CFractionCoordinates::toFractionFromCart(const double x, const double y, const double z)
{
    Point3 p1;
    p1<<x,y,z;
    return toFractionFromCart(p1);
}
CCartesianCoordinates* CFractionCoordinates::toCartesianCoordinates() const
{
   CCartesianCoordinates* temp_CartesianCoordinates = new CCartesianCoordinates(this->m_pPeriodicFramework);
   for(size_t i=0;i<this->m_coordinates.size();i++)
   {
       Point3 temP = std::move(toCartesian(i));
       temp_CartesianCoordinates->append(temP);
   }

   return temp_CartesianCoordinates;
}
//
CInternalCoordinates*  CFractionCoordinates::toInternalCoordinates()  const
 {
     return toCartesianCoordinates()->toInternalCoordinates();
 }
CFractionCoordinates& CFractionCoordinates::operator = (const CFractionCoordinates& othr)
{
     this->m_pPeriodicFramework=othr.m_pPeriodicFramework;
     this->m_coordinates.assign(othr.m_coordinates.begin(),othr.m_coordinates.end());
     return *this;
}
CFractionCoordinates CFractionCoordinates::add(const CFractionCoordinates& othr) const
{
     if(this->size()==othr.size() &&
        *(m_pPeriodicFramework->unitcell()) == *(othr.m_pPeriodicFramework->unitcell()))
     {
         CFractionCoordinates resultFractCoord(m_pPeriodicFramework);
         Point3 resultP;
         for(size_t i=0;i<size();i++)
         {
             resultP = m_coordinates[i] + othr.m_coordinates[i];
             resultFractCoord.append(resultP);
         }
         return resultFractCoord;
     }else
         return nullptr;

}
CFractionCoordinates CFractionCoordinates::substract(const CFractionCoordinates& othr)const
{
    if(this->size()==othr.size() &&
        *(m_pPeriodicFramework->m_pUnitCell) == *(othr.m_pPeriodicFramework->m_pUnitCell))
    {
             CFractionCoordinates resultFractCoord(m_pPeriodicFramework);
             Point3 resultP;
             for(size_t i=0;i<size();i++)
             {
                 resultP = m_coordinates[i] - othr.m_coordinates[i];
                 resultFractCoord.append(resultP);
             }
             return resultFractCoord;
    }else
       return nullptr;

}
Eigen::Matrix<double,3,3> CFractionCoordinates::multiply(const CFractionCoordinates& othr)const
{
    Eigen::Matrix<double,3,3> resultMat;
    resultMat = (toCartesianCoordinates()->toMatrix().transpose())*(othr.toCartesianCoordinates()->toMatrix().transpose());
    return resultMat;
}
CFractionCoordinates  CFractionCoordinates::operator + (const CFractionCoordinates& othr) const
{
   return add(othr);
}
CFractionCoordinates  CFractionCoordinates::operator - (const CFractionCoordinates& othr) const
{
   return substract(othr);
}
const Point3& CFractionCoordinates::operator[](size_t index) const
{
      if(index >= this->m_coordinates.size()){
        Log::Error<<"Index exceed the size of vector in CFractionCoordinates::operator!!!"<<std::endl;
        boost::throw_exception(std::runtime_error("Index exceed the size of vector in CFractionCoordinates::operator!!!"));
      }
      return m_coordinates[index];
}



}// namespace of CATAZJUT

