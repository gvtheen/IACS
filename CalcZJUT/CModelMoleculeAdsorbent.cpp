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
#include "CModelMoleculeAdsorbent.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../CataZJUT/CAtom.h"
#include "../IACS.h"
#include "../Util/foreach.h"
#include "../Util/log.hpp"

using util::Log;
namespace CALCZJUT{


CModelMoleculeAdsorbent::CModelMoleculeAdsorbent(CATAZJUT::CConfigurationBase* mpConf,
                                                 Bitset& molIndexBit)
:m_pConfiguration(mpConf),m_AtomicBits(molIndexBit)
{
    assert(m_pConfiguration);
    #ifdef DEBUG
        Log::Debug<<"CModelMoleculeAdsorbent::CModelMoleculeAdsorbent" << std::endl;
    #endif
    foreach(CATAZJUT::CAtom* atom,m_pConfiguration->atoms())
          if(m_AtomicBits.test(atom->index())==1)
             m_Atoms.push_back(atom);
}
CModelMoleculeAdsorbent::CModelMoleculeAdsorbent(CATAZJUT::CConfigurationBase* mpConf,
                                                 std::vector<CATAZJUT::CAtom*>& molAtoms)
:m_pConfiguration(mpConf),m_Atoms(molAtoms)
{
      m_AtomicBits.resize(m_pConfiguration->atomCount(),false);
      foreach(CATAZJUT::CAtom* atom, m_Atoms){
          m_AtomicBits.set(atom->index(),true);
      }
}
CModelMoleculeAdsorbent::~CModelMoleculeAdsorbent()
{
    //dtor
}

CATAZJUT::CConfigurationBase* CModelMoleculeAdsorbent::configuration() const
{
    return this->m_pConfiguration;
}
void CModelMoleculeAdsorbent::setConfiguration(CATAZJUT::CConfigurationBase* mbf)
{
    this->m_pConfiguration=mbf;
    foreach(CATAZJUT::CAtom* atom,m_pConfiguration->atoms())
          if(m_AtomicBits.test(atom->index())==1)
             m_Atoms.push_back(atom);
}
CModelMoleculeAdsorbent::AtomRange CModelMoleculeAdsorbent::atoms() const
{
    return boost::make_iterator_range(this->m_Atoms);
}
CATAZJUT::CAtom* CModelMoleculeAdsorbent::atom(size_t index) const
{
    return m_Atoms[index];
}
Bitset CModelMoleculeAdsorbent::bitSet()const
{
     return this->m_AtomicBits;
}
size_t CModelMoleculeAdsorbent::atomCount() const
{
    return m_Atoms.size();
}
void CModelMoleculeAdsorbent::moveBy( const Point3 vect )
{
   Point3 p1;
   foreach( CATAZJUT::CAtom* atom, atoms()){
       p1=atom->position()+vect;
       atom->SetPosition(p1);
   }
}
void CModelMoleculeAdsorbent::moveBy(double x,double y,double z)
{
    Point3 vect;
    vect<<x,y,z;
    moveBy(vect);
}
void CModelMoleculeAdsorbent::rotate(const Point3& vect, const double& angle_Radian)
{
   //double angle_Radian = angle * CATAZJUT::constants::DegreesToRadians;

   Point3 displaceVector = -1*gravityCentre();
   this->moveBy(displaceVector);

   Eigen::Matrix<double, 3, 1> axisVector(vect.x(), vect.y(), vect.z());
   Eigen::Transform<double, 3, 3> transform(Eigen::AngleAxis<double>(angle_Radian, axisVector));
   Point3 p1;
   foreach( CATAZJUT::CAtom* atom, atoms()){
       p1=transform*(atom->position());
       atom->SetPosition(p1);
   }

   this->moveBy(-1*displaceVector);

}
Point3 CModelMoleculeAdsorbent::gravityCentre()
{
   Point3 temp;
   temp<<0,0,0;
   foreach(const CATAZJUT::CAtom* atom, atoms()){
      temp=temp+atom->position();
   }
   return temp/atomCount();
}


}
