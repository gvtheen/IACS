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
#include <boost/bind.hpp>
#include "CAtom.h"
#include "CConfigurationBase.h"
#include "CConfigurationPrivateData.h"
#include "CCartesianCoordinates.h"
#include "CBond.h"
#include "Geometry.h"
#include "CElement.h"
#include "../Util/foreach.h"
#include "CFragment.h"

namespace CATAZJUT{

CAtom::CAtom(CConfigurationBase* tempConf,size_t tempIndex)
:m_pConfiguration(tempConf),m_index(tempIndex)
{

}
CAtom::~CAtom()
{
    //dtor
}
void CAtom::SetPosition(const Point3 &position)
{
    this->m_pConfiguration->coordinates()->setPosition(m_index,position);
}
void CAtom::SetPosition(double x, double y, double z)
{
    Point3 temp;
    temp<<x,y,z;
    this->SetPosition(temp);
}
Point3 CAtom::position() const
{
    return this->m_pConfiguration->coordinates()->position(m_index);
}
double CAtom::x() const
{
    return this->position().x();
}
double CAtom::y() const
{
    return this->position().y();
}
double CAtom::z() const
{
    return this->position().z();
}
//void CAtom::setMassNumber(MassNumberType massNumber)
//{
//
//}
void CAtom::setPartialCharge(double charge)
{
    m_pConfiguration->m_pData->partialCharges[m_index] = charge;
}
double CAtom::PartialCharge() const
{
    return m_pConfiguration->m_pData->partialCharges[m_index];
}
double CAtom::distance(const CAtom *atom) const
{
    return CATAZJUT::Geometry::distance(position(),atom->position());
}
std::string CAtom::Symbol() const
{
    return m_pConfiguration->m_Element[this->m_index].symbol();
}
std::string CAtom::Name() const
{
    return m_pConfiguration->m_Element[this->m_index].name();
}
double CAtom::Mass() const
{
    return element().mass();
}
double CAtom::CovalentRadius() const
{
    return element().covalentRadius();
}
double CAtom::VanDerWaalsRadius() const
{
    return element().vanDerWaalsRadius();
}
double CAtom::AtomNumber() const
{
    return element().atomicNumber();
}
int CAtom::ExpectedValence() const
{
    return element().expectedValence();
}
double CAtom::Electronegativity() const
{
    return element().electronegativity();
}
size_t CAtom::maxCoordinationNum()const
{
    return element().maxCoordinationNum();
}
std::string CAtom::valentConfigurationStr()const
{
    return element().valentConfigurationStr();
}


CAtom::NeighborRange CAtom::neighbors() const
{
    const std::vector<CBond *> &bonds = m_pConfiguration->m_pData->atomBonds[m_index];

    return boost::make_iterator_range(
                boost::make_transform_iterator(
                    bonds.begin(), boost::bind(&CBond::otherAtom, _1, this)),
                boost::make_transform_iterator(
                      bonds.end(), boost::bind(&CBond::otherAtom, _1, this)));
}
CAtom::BondRange CAtom::bonds() const
{
    const std::vector<CBond *> &bonds = m_pConfiguration->m_pData->atomBonds[m_index];
    return boost::make_iterator_range(bonds.begin(), bonds.end());
}
CBond* CAtom::bond(size_t index) const
{
    return bonds()[index];
}
size_t CAtom::bondCount() const
{
    return bonds().size();
}
CBond* CAtom::bondTo(const CAtom *atom) const
{
    foreach(CBond* mbond, bonds())
       if(mbond->otherAtom(this)==atom)
          return mbond;
    //if not exist
    return nullptr;

}
CAtom* CAtom::neighbor(size_t index) const
{
    return neighbors()[index];
}
size_t CAtom::neighborCount() const
{
    return bondCount();
}
size_t CAtom::neighborCount(const CElement &element) const
{
    size_t num=0;
    foreach(CBond *mbond,bonds())
       if(mbond->otherAtom(this)->IsElement(element))
          num++;
    return num;
}
bool CAtom::isBondedTo(const CAtom *atom) const
{
    return bondTo(atom) != nullptr;
}
bool CAtom::isBondedTo(const CElement &element) const
{
    foreach(CBond *mbond,bonds())
       if(mbond->otherAtom(this)->IsElement(element))
          return true;
    return false;
}
bool CAtom::isTerminal() const
{
    return neighborCount() == 1;
}
bool CAtom::isTerminalHydrogen() const
{
    return isTerminal() && this->IsElement(*(new CElement(Hydrogen)));
}
CFragment* CAtom::fragment() const
{
    return this->m_pConfiguration->fragmentForAtom(this);
}
bool CAtom::isBondSaturated()const
{
    return element().maxCoordinationNum() <= bondCount();
}
Vector3 CAtom::NewBondingVect()const
{
    Vector3 res(0,0,0);
    foreach(CBond *bond,this->bonds())
         res = res + (bond->bondVector(this)).normalized();
    return -1* res.normalized();
}

}
