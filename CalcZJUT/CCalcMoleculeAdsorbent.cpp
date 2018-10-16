#include "CCalcMoleculeAdsorbent.h"
#include "CConfigurationBase.h"
#include "CAtom.h"
#include "foreach.h"

namespace CALCZJUT{


CCalcMoleculeAdsorbent::CCalcMoleculeAdsorbent(CConfigurationBase* mpConf,Bitset& molIndexBit)
:m_pConfiguration(mpConf),m_AtomicBits(molIndexBit)
{
    assert(m_pConfiguration);

    for(size_t i=0;i<m_AtomicBits.size();i++)
        if(m_AtomicBits.test(i)==1){
          m_Atoms.push_back(m_pConfiguration->m_Atom[i]);
        }

}
CCalcMoleculeAdsorbent::CCalcMoleculeAdsorbent(CConfigurationBase* mpConf,std::vector<CAtom*>& molAtoms)
:m_pConfiguration(mpConf),m_Atoms(molAtoms);
{
      m_AtomicBits.resize(m_pPeriodicFramework->atomCount(),false);
      foreach(CAtom* atom, m_Atoms){
          m_AtomicBits.set(atom->index(),true);
      }
}
CCalcMoleculeAdsorbent::~CCalcMoleculeAdsorbent()
{
    //dtor
}

CConfigurationBase* CCalcMoleculeAdsorbent::configuration() const
{
    return this->m_pConfiguration;
}
void CCalcMoleculeAdsorbent::setConfiguration(CATAZJUT::CConfigurationBase* mbf)
{
    this->m_pConfiguration=mbf;
    for(size_t i=0;i<m_AtomicBits.size();i++)
        if(m_AtomicBits.test(i)==1){
          m_Atoms.push_back(m_pConfiguration->m_Atom[i]);
        }
}
CCalcMoleculeAdsorbent::AtomRange CCalcMoleculeAdsorbent::atoms() const
{
    return boost::make_iterator_range(this->m_Atoms);
}
CAtom* CCalcMoleculeAdsorbent::atom(size_t index) const
{
    return m_Atoms[index];
}
Bitset CCalcMoleculeAdsorbent::bitSet()const
{
     return this->m_AtomicBits;
}
size_t CCalcMoleculeAdsorbent::atomCount() const
{
    return m_Atoms.size();
}
void CCalcMoleculeAdsorbent::moveBy( const Point3 vect )
{
   Point3 p1;
   foreach(const CAtom* atom, atoms()){
       p1=atom->position()+vect;
       atom->SetPosition(p1);
   }
}
void CCalcMoleculeAdsorbent::moveBy(double x,double y,double z)
{
    Point3 vect;
    vect<<x,y,z;
    moveBy(vect);
}
void CCalcMoleculeAdsorbent::rotate(const Point3& vect, const double& angle_Radian)
{
   //double angle_Radian = angle * CATAZJUT::constants::DegreesToRadians;

   Point3 displaceVector = -1*gravityCentre();
   this->moveBy(displaceVector);

   Eigen::Matrix<double, 3, 1> axisVector(vect.x(), vect.y(), vect.z());
   Eigen::Transform<double, 3, 3> transform(Eigen::AngleAxis<double>(angle_Radian, axisVector));
   Point3 p1;
   foreach(const CAtom* atom, atoms()){
       p1=transform*(atom->position());
       atom->SetPosition(p1);
   }

   this->moveBy(-1*displaceVector);

}
Point3 CCalcMoleculeAdsorbent::gravityCentre()
{
   Point3 temp;
   temp<<0,0,0;
   foreach(cost CAtom* atom, atoms()){
      temp=temp+atom->position();
   }
   return temp/atomCount();
}


}
