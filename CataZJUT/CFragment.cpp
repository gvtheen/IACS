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
#include "CFragment.h"
#include <cassert>
#include "CAtom.h"
#include "foreach.h"
#include "Cbond.h"
#include "CConfigurationBase.h"

namespace CATAZJUT{

CFragment::CFragment(CConfigurationBase* mConfig,const Bitset& bit_set)
:m_pConfiguration(mConfig),m_bitset(bit_set)
{
    //ctor
}
CFragment::~CFragment()
{
    //dtor
}
std::vector<CBond*> CFragment::bonds() const
{
    std::vector<CBond*> bonds;
    foreach(CAtom* atom, atoms()){
        foreach(CBond* bond, atom->bonds()){
           if(std::find(bonds.begin(),bonds.end(),bond)==bonds.end())
               bonds.push_back(bond);
        }
    }
    return bonds;
}
size_t CFragment::bondCount() const
{
   return bonds().size();
}
bool CFragment::contains(const CBond *bond) const
{
   return bond->atom1()->fragment()==this;
}
Bitset CFragment::bitSet()const
{
    return this->m_bitset;
}
std::vector<Point3> CFragment::coordinates() const
{
    std::vector<Point3> coordinate_res;
    foreach(CAtom* atom, atoms())
        coordinate_res.push_back(atom->position());
    return coordinate_res;
}
void CFragment::move(const Vector3& vect_R)
{
    foreach(CAtom* atom, atoms())
       atom->SetPosition( atom->position() + vect_R );
}
bool CFragment::isBondTo(const CFragment* otherFragment)
{
    assert(otherFragment);

    if(this->m_pConfiguration != otherFragment->Configuration())
        return false;

    foreach(CAtom* atom_s, this->atoms())
      foreach(CAtom* atom_m, otherFragment->atoms())
      {
          if(this->m_pConfiguration->isBondBetween(atom_s,atom_m))
             return true;
      }
    return false;
}


}
