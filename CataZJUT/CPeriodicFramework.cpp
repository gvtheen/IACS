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
#include "CPeriodicFramework.h"
#include "CConfigurationPrivateData.h"
#include "CFractionCoordinates.h"
#include "CUnitCell.h"
#include "../CalcZJUT/CParameter.h"
#include "CMolecularSymmetry.h"
#include "CSpaceGroup.h"

namespace CATAZJUT{

CPeriodicFramework::CPeriodicFramework(CALCZJUT::CParameter* mpa)
:CConfigurationBase(mpa)
{
    this->m_pUnitCell = new CUnitCell();
}
CPeriodicFramework::CPeriodicFramework(CPeriodicFramework& other)
:CConfigurationBase(other)
{
    if(other.dimensionalType()==CATAZJUT::DEFINED::Periodic)
    {
        CUnitCell* temp=other.unitcell();
        this->m_pUnitCell = new CUnitCell(*temp);
    }else
        this->m_pUnitCell = nullptr;
}
CPeriodicFramework::CPeriodicFramework(CConfigurationBase& other)
:CConfigurationBase(other)
{
       this->m_pUnitCell = nullptr;
}
CPeriodicFramework::~CPeriodicFramework()
{

    delete m_pUnitCell;
}
CPeriodicFramework* CPeriodicFramework::clone()
{
    return new CPeriodicFramework(*this);
}
CUnitCell* CPeriodicFramework::unitcell()
{
    return this->m_pUnitCell;
}

CFractionCoordinates* CPeriodicFramework::Fractioncoordinates()
{
    int returnindex=-1;
    if(m_pData->coordinateSets.empty()|| m_pData->coordinateSets.front()->type()\
                                               == CCoordinateSet::None){
            // if no, construct new coordinate
            CFractionCoordinates* temp_Internal = new CFractionCoordinates(this,atomCount());
            m_pData->coordinateSets.push_back(boost::make_shared<CCoordinateSet>(temp_Internal));
            returnindex = 0;
    }else{
        for(size_t i=0;i<m_pData->coordinateSets.size();i++)
          if(m_pData->coordinateSets[i]->type() == CCoordinateSet::Fraction){
             returnindex = i;
             break;
          }
        if(returnindex == -1)
        {
            CFractionCoordinates* temp_Internal= new CFractionCoordinates(this,atomCount());
            m_pData->coordinateSets.push_back(boost::make_shared<CCoordinateSet>(temp_Internal));
            returnindex = m_pData->coordinateSets.size() - 1;
        }
    }
    return m_pData->coordinateSets[returnindex]->fractionCoordinates();
}
std::string CPeriodicFramework::SymmetrySymbol()
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

}
