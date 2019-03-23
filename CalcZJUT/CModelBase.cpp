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
#include <string>
#include <sstream>
#include "CModelBase.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../CataZJUT/CatalystUniverseDefine.h"
#include "../CataZJUT/CUnitCell.h"
#include "../CataZJUT/CAtom.h"
#include "../Util/log.hpp"
#include "../Util/foreach.h"
#include "../Util/Point-Vector.h"
#include "CIOCar.h"

using util::Log;
using util::Point3;

namespace CALCZJUT{

CModelBase::CModelBase(CParameter* temParameter,size_t index)
:m_pParameter(temParameter),m_index(index)
{
    m_IsNeedRandomInit=true;
}
CModelBase::CModelBase(CModelBase& obj)
{
   this->m_pParameter=obj.m_pParameter;

//    #ifdef DEBUG
//        Log::Debug<<"CModelBase::CModelBase(CModelBase& obj)" << std::endl;
//    #endif // DEBU
   this->m_pPeriodicFramework=new CATAZJUT::CConfigurationBase(*(obj.m_pPeriodicFramework));
   this->m_ppBackupPeriodicFramework=obj.m_ppBackupPeriodicFramework;
}
CModelBase::~CModelBase()
{
    //dtor
}
void CModelBase::init()
{

}
CATAZJUT::CConfigurationBase* CModelBase::periodicFramework()
{
    return m_pPeriodicFramework;
}
void CModelBase::setPeriodicFramekwork(CATAZJUT::CConfigurationBase* mbf)
{
    m_pPeriodicFramework = mbf;
}
void CModelBase::createSupport(const Bitset & mt)
{
    // no doing
}
void CModelBase::createMoleAdsorb(const Bitset & mt)
{
    // no doing
}
//void CModelBase::createStructureAtGene()
//{
//
//}
void CModelBase::setRandomInitState(const bool& mht)
{
    this->m_IsNeedRandomInit=mht;
}
bool CModelBase::RandomInitState()
{
    return this->m_IsNeedRandomInit;
}
//void CModelBase::removeStructureOfGene()
//{
//   for(size_t i=0;i<m_PopuPeriodicFramework.size();i++)
//       delete m_PopuPeriodicFramework[i];
//   m_PopuPeriodicFramework.clear();
//}

std::vector<std::pair<std::string,size_t>>& CModelBase::chemicalFormula()
{
     return this->m_chemicalFormula;
}
void CModelBase::setChemicalFormula(const std::vector<std::pair<std::string,size_t>>& mth)
{

}

size_t CModelBase::index()
{
    return this->m_index;
}
void CModelBase::setIndex(size_t index)
{
     this->m_index=index;
}
void CModelBase::outputStructureToFile()
{
   CIOCar* io = new CIOCar(this->m_pPeriodicFramework);
   std::stringstream str ;
   str<<"Struct_Gene_" << m_pParameter->currentGenerationNum()<<"_Pop_" << this->index();
   std::string filename;
   str>>filename;
   io->output(filename);
   delete io;
}
void CModelBase::standardOutput(size_t type)
{
/** \brief   output the most stable structure of current population.
 *
 * \type = 1  : output the most stable structure of current population
 * \type = 0  : output each  structure of current population
 * \type = 2  : output the most stable structure after the whole calculation
 *
 */
    switch (type){
      case 1:
         Log::Output<<"*****The most stable structure in ";
         Log::Output<<m_pParameter->currentGenerationNum()<<"th Generation"<<"*****"<<std::endl;
         break;
      case 0:
         Log::Output<<"*****The relaxed structure in " << m_pParameter->popNum() << "th Population ";
         Log::Output<<m_pParameter->currentGenerationNum()<<"th Generation"<<"*****"<<std::endl;
         break;
      case 2:
         Log::Output<<"*****The most stable structures after calculation"<<"*****"<<std::endl;
         break;
      default:
         break;
    }
    if(m_pPeriodicFramework->dimensionalType()==CATAZJUT::DEFINED::Periodic){
       Eigen::Matrix<double, 3, 3> mat = m_pPeriodicFramework->unitcell()->MatrixOfBravaisLattice();
       Log::Output<<"********************Unit cell parameter****************************"<<std::endl;
       for(size_t i=0;i<3;i++)
          Log::Output<<mat(i,0)<<"    "<<mat(i,1)<<"   "<<mat(i,2)<<std::endl;
       Log::Output<<"*******************************************************************"<<std::endl;
    }
    Log::Output<<"-------------------------------------------------------------------"<<std::endl;
    Log::Output<<"Element               X                Y              Z            "<<std::endl;
    Log::Output<<"-------------------------------------------------------------------"<<std::endl;
    Point3 poi;
    foreach(CATAZJUT::CAtom* atom,m_pPeriodicFramework->atoms()){
       poi = atom->position();
       Log::Output<<atom->Symbol()<<"    ";
       Log::Output<<poi[0]<<"    "<<poi[1]<<"    "<<poi[2]<<std::endl;
    }
    Log::Output<<"-------------------------------------------------------------------"<<std::endl;
    Log::Output<<m_pPeriodicFramework->SymmetrySymbol()<<std::endl;
    Log::Output<<std::endl;
}

}//namespace
