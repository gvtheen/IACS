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
#include <boost/make_shared.hpp>
#include "CConfigurationBase.h"
#include "CConfigurationPrivateData.h"
#include "CFractionCoordinates.h"
#include "CUnitCell.h"
#include "../CalcZJUT/CParameter.h"
#include "CMolecularSymmetry.h"
#include "CSpaceGroup.h"

namespace CATAZJUT{

CConfigurationBase::CConfigurationBase(CALCZJUT::CParameter* mpa)
:CConfigurationBase(mpa)
{
    this->m_pUnitCell = new CUnitCell();
}
CConfigurationBase::CConfigurationBase(CConfigurationBase& other)
:CConfigurationBase(other)
{
    if(other.dimensionalType()==CATAZJUT::DEFINED::Periodic)
    {
        CUnitCell* temp=other.unitcell();
        this->m_pUnitCell = new CUnitCell(*temp);
    }else
        this->m_pUnitCell = nullptr;
}
CConfigurationBase::CConfigurationBase(CConfigurationBase& other)
:CConfigurationBase(other)
{
       this->m_pUnitCell = nullptr;
}
CConfigurationBase::~CConfigurationBase()
{

    delete m_pUnitCell;
}
CConfigurationBase* CConfigurationBase::clone()
{
    return new CConfigurationBase(*this);
}
CUnitCell* CConfigurationBase::unitcell()
{
    if(this->m_pUnitCell==nullptr)
        this->m_pUnitCell=new CUnitCell();

    return this->m_pUnitCell;
}


std::string CConfigurationBase::SymmetrySymbol()
{
    if(m_pParameter->simulationMode==CALCZJUT::CParameter::CLUSTER ||
       m_pParameter->simulationMode == CALCZJUT::CParameter::MOL_CLUSTER ){
        CMolecularSymmetry* molsym = new CMolecularSymmetry(this);
        std::string res("Molecular symmetry: point group ");
        res=res + molsym->GetPointGroup();
        delete molsym;
        return res;
    }else{
        CSpaceGroup* spaceGroup = new CSpaceGroup(this);
        std::string res("Crystallography symmetry: space group ");
        res=res + spaceGroup->GetSpaceGroup();
        delete spaceGroup;
        return res;
    }
}
void CConfigurationBase::perceiveBonds()
{
    CConfigurationBase::perceiveBonds();
}



}
