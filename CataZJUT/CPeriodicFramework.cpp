#include <boost/scoped_ptr.hpp>
#include <boost/make_shared.hpp>
#include "CPeriodicFramework.h"
#include "CConfigurationPrivateData.h"
#include "CFractionCoordinates.h"
#include "CUnitCell.h"

namespace CATAZJUT{

CPeriodicFramework::CPeriodicFramework(CALCZJUT::CParameter* mPa)
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
            CFractionCoordinates* temp_Internal= new CFractionCoordinates(this,atomCount());
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


}
