#include "CCalcSupportBase.h"
#include "CConfigurationBase.h"
#include "CPeriodicFramework.h"
#include "CCartesianCoordinates.h"
#include "foreach.h"
#include "CAtom.h"
namespace CALCZJUT{

CCalcSupportBase::CCalcSupportBase(CATAZJUT::CPeriodicFramework* mpconf, Bitset& supportBit)
:m_pPeriodicFramework(mpconf),m_AtomicBits(supportBit)
{
      for(size_t i=0;i<m_AtomicBits.size();i++)
         if(m_AtomicBits.test(i)==1){
            this->m_Atoms.push_back(m_pPeriodicFramework->m_Atom[i]);
         }
}
CCalcSupportBase::CCalcSupportBase(CATAZJUT::CPeriodicFramework* mpconf, std::vector<CATAZJUT::CAtom*> support)
:m_pPeriodicFramework(mpconf),m_Atoms(support)
{
      m_AtomicBits.resize(m_pPeriodicFramework->atomCount(),false);
      foreach(CATAZJUT::CAtom* atom, m_Atoms){
          m_AtomicBits.set(atom->index(),true);
      }
}

size_t CCalcSupportBase::size() const
{
    return this->atoms().size();
}
bool CCalcSupportBase::isEmpty() const
{
    return size()==0;
}
CATAZJUT::CAtom* CCalcSupportBase::atom(size_t index) const
{
    if(m_AtomicBits.test(index)!=1)
        return nullptr;
    else
        return m_pPeriodicFramework->m_Atom[index];
}
CATAZJUT::AtomRange CCalcSupportBase::atoms() const
{
    return boost::make_iterator_range(this->m_Atoms);
}
size_t CCalcSupportBase::atomCount() const
{
    return m_Atoms.size();
}
CATAZJUT::CConfigurationBase* CCalcSupportBase::configuration() const
{
    return this->m_pPeriodicFramework;
}

void CCalcSupportBase::setConfiguration(CATAZJUT::CPeriodicFramework* mbf)
{
    this->m_pPeriodicFramework=mbf;
    for(size_t i=0;i<m_AtomicBits.size();i++)
         if(m_AtomicBits.test(i)==1){
            this->m_Atoms.push_back(m_pPeriodicFramework->m_Atom[i]);
         }
}
CATAZJUT::CCartesianCoordinates* CCalcSupportBase::coordinates()
{
    return m_pPeriodicFramework->coordinates();
}
Bitset CalcSupportBase::bitSet()const
{
    return this->m_AtomicBits;
}
CCalcSupportBase::~CCalcSupportBase()
{
    //dtor
}

}
