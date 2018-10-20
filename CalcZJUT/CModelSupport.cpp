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
#include "CModelSupport.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../CataZJUT/CCartesianCoordinates.h"
#include "../Util/foreach.h"
#include "../CataZJUT/CAtom.h"
namespace CALCZJUT{

CModelSupport::CModelSupport(CATAZJUT::CPeriodicFramework* mpconf, Bitset& supportBit)
:m_pPeriodicFramework(mpconf),m_AtomicBits(supportBit)
{
      for(size_t i=0;i<m_AtomicBits.size();i++)
         if(m_AtomicBits.test(i)==1){
            this->m_Atoms.push_back(m_pPeriodicFramework->m_Atom[i]);
         }
}
CModelSupport::CModelSupport(CATAZJUT::CPeriodicFramework* mpconf, std::vector<CATAZJUT::CAtom*> support)
:m_pPeriodicFramework(mpconf),m_Atoms(support)
{
      m_AtomicBits.resize(m_pPeriodicFramework->atomCount(),false);
      foreach(CATAZJUT::CAtom* atom, m_Atoms){
          m_AtomicBits.set(atom->index(),true);
      }
}

size_t CModelSupport::size() const
{
    return this->atoms().size();
}
bool CModelSupport::isEmpty() const
{
    return size()==0;
}
CATAZJUT::CAtom* CModelSupport::atom(size_t index) const
{
    if(m_AtomicBits.test(index)!=1)
        return nullptr;
    else
        return m_pPeriodicFramework->m_Atom[index];
}
CModelSupport::AtomRange CModelSupport::atoms() const
{
    return boost::make_iterator_range(this->m_Atoms);
}
size_t CModelSupport::atomCount() const
{
    return m_Atoms.size();
}
CATAZJUT::CConfigurationBase* CModelSupport::configuration() const
{
    return this->m_pPeriodicFramework;
}

void CModelSupport::setConfiguration(CATAZJUT::CPeriodicFramework* mbf)
{
    this->m_pPeriodicFramework=mbf;
    for(size_t i=0;i<m_AtomicBits.size();i++)
         if(m_AtomicBits.test(i)==1){
            this->m_Atoms.push_back(m_pPeriodicFramework->m_Atom[i]);
         }
}
CATAZJUT::CCartesianCoordinates* CModelSupport::coordinates()
{
    return m_pPeriodicFramework->coordinates();
}
Bitset CModelSupport::bitSet()const
{
    return this->m_AtomicBits;
}
CModelSupport::~CModelSupport()
{
    //dtor
}

}
