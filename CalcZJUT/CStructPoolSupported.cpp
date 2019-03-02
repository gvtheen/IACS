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
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <strings.h>
#include <cstring>
#include <math.h>
#include "../CataZJUT/CConfigurationBase.h"
#include "CModelCluster.h"
#include "CParameter.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CCartesianCoordinates.h"
#include "../Util/log.hpp"
#include "../Util/foreach.h"
#include "../Util/Bitset.h"
#include "../Util/Point-Vector.h"
#include "../CataZJUT/Constant.h"
#include "../CataZJUT/CFragment.h"
#include "../Util/utilFunction.h"
#include "../Util/CRandomgenerator.h"
#include "../GaZJUT/CGaparameter.h"
#include "../CataZJUT/CElement.h"
#include "CIOBase.h"
#include "CIOMol.h"
#include "CIOCar.h"
#include "CIOCif.h"
#include "CIOGjf.h"
#include "CIOXyz.h"
#include "CIOCellFile.h"
#include "CIOPoscar.h"
#include "CStructPoolSupported.h"
#include "CModel2DSupport.h"
#include "CModelClusterSupport.h"
#include "CModelClusterLoaded2DSupport.h"

using util::Log;

namespace CALCZJUT{

CStructPoolSupported::CStructPoolSupported(CParameter* oth)
:CStructPoolBase(oth)
{
    for(size_t i=0;i<this->m_pParameter->GaParameter()->PopNum();i++){
        switch ( (int)(m_pParameter->simulationMode) ){
            case CParameter::MOL_2DMATERIAL:
                m_CalcStructPool.push_back(new CModel2DSupport(this->m_pParameter,&copy_pPeriodicFramework,i));
                break;
            case CParameter::MOL_CLUSTER2DMATERIAL:
                m_CalcStructPool.push_back(new CModelClusterLoaded2DSupport(this->m_pParameter,&copy_pPeriodicFramework,i));
                break;
            case CParameter::MOL_CLUSTER:
                m_CalcStructPool.push_back(new CModelClusterSupport(this->m_pParameter,&copy_pPeriodicFramework,i));
            default:
                break;
        }
        m_CalcStructPool[m_CalcStructPool.size()-1]->periodicFramework()->setExcludeBond(m_pParameter->excludeBond);
        m_CalcStructPool[m_CalcStructPool.size()-1]->periodicFramework()->setTolerancefactor(m_pParameter->bondToleranceFactor);
    }

}

CStructPoolSupported::~CStructPoolSupported()
{
    //dtor
}
void CStructPoolSupported::init()
{
       if(this->m_pParameter->adso_supp_Input_File.size()!= 0){
           //read coordinate of mixed adso_supp_Struct, add the pointer of m_pCalcModeStruct;
           CATAZJUT::CFragment *m1,*m2;
           for(size_t i=0;i<m_pParameter->adso_supp_Input_File.size();i++){
               if(i >= m_CalcStructPool.size())
                  break;

                this->getIO(*(m_pParameter->adso_supp_Input_File[i]),m_CalcStructPool[i]->periodicFramework());
                this->m_IO->input(*(m_pParameter->adso_supp_Input_File[i]));

                m_CalcStructPool[i]->periodicFramework()->perceiveBonds();
                m_CalcStructPool[i]->periodicFramework()->perceiveFragments();
                this->m_IO->output("moltest-iacs");
                #ifdef DEBUG
                     Log::Debug<<"*fragment num:" <<m_CalcStructPool[i]->periodicFramework()->fragmentNum()<< std::endl;
                #endif // DEBU
                m1 = m_CalcStructPool[i]->periodicFramework()->fragment(0);
                m2 = m_CalcStructPool[i]->periodicFramework()->fragment(1);
               if(m1== nullptr || m2 ==nullptr){
                  Log::Error<<"Support and adsorbent setting of input file are error! CStructPoolSupported::init().\n";
                  boost::throw_exception(std::runtime_error("Support and adsorbent setting are error! Check the file: CStructPoolSupported::init()."));
               }
               if(m1->atomCount() > m2->atomCount())
               {
                  m_CalcStructPool[i]->createSupport(m1->bitSet());
                  #ifdef DEBUG
                     Log::Debug<<"1-CStructPoolSupported::init()" << std::endl;
                  #endif // DEBU
                  m_CalcStructPool[i]->createMoleAdsorb(m2->bitSet());
               }else{
                  #ifdef DEBUG
                     Log::Debug<<"2-CStructPoolSupported::init()" << std::endl;
                  #endif // DEBU
                  m_CalcStructPool[i]->createSupport(m2->bitSet());
                  m_CalcStructPool[i]->createMoleAdsorb(m1->bitSet());
               }
               #ifdef DEBUG
                     Log::Debug<<"CStructPoolSupported::init()" << std::endl;
               #endif // DEBU
               m_CalcStructPool[i]->setRandomInitState(false);
               m_CalcStructPool[i]->setIndex(i);
          }
          this->copy_pPeriodicFramework = m_CalcStructPool[0]->periodicFramework()->clone();

          for(size_t i=m_pParameter->adso_supp_Input_File.size();i<m_CalcStructPool.size();i++){
               delete m_CalcStructPool[i];
               m_CalcStructPool[i]=m_CalcStructPool[0]->clone();
               // get support and adsorption molecular, then re-set new structure.
               #ifdef DEBUG
                     Log::Debug<<"4-CStructPoolSupported::init()" << std::endl;
               #endif // DEBU
               m_CalcStructPool[i]->setRandomInitState(true);
               m_CalcStructPool[i]->setIndex(i);
          }
      }else{ // treat cluster-support
            //read coordinate of mixed adso_supp_Struct, add the pointer of m_pCalcModeStruct;
           //
           this->getIO(m_pParameter->supportStructFile,m_CalcStructPool[0]->periodicFramework());
           this->m_IO->input(m_pParameter->supportStructFile);


           this->getIO(m_pParameter->adsorbentStructFile,m_CalcStructPool[0]->periodicFramework());
           Bitset Bit_adsorbent = this->m_IO->input(m_pParameter->adsorbentStructFile,CParameter::MOL_CLUSTER);
           // get opposite bit for support
           m_CalcStructPool[0]->createSupport(~Bit_adsorbent);
           //set molecular adsorbent bits
           m_CalcStructPool[0]->createMoleAdsorb(Bit_adsorbent);
           m_CalcStructPool[0]->setRandomInitState(false);

           this->copy_pPeriodicFramework = m_CalcStructPool[0]->periodicFramework()->clone();

           for(size_t i=1;i<m_CalcStructPool.size();i++){
               delete m_CalcStructPool[i];
               m_CalcStructPool[i]=m_CalcStructPool[0]->clone();
               // get support and adsorption molecular, then re-set new structure.
               m_CalcStructPool[i]->setRandomInitState(true);
               m_CalcStructPool[i]->setIndex(i);
           }
     }
      #ifdef DEBUG
                     Log::Debug<<"5-CStructPoolSupported::init()" << std::endl;
      #endif // DEBU
}
void CStructPoolSupported::getIO(std::string &file_name,CATAZJUT::CConfigurationBase* currentPeriodicFramework)
{
    std::vector<std::string> vectstr;
    boost::algorithm::split(vectstr,file_name,boost::algorithm::is_any_of("."),boost::algorithm::token_compress_on);
    boost::algorithm::trim(vectstr[1]);
    if(vectstr[1]=="mol")
        this->m_IO = new CIOMol(currentPeriodicFramework);
    else if(vectstr[1]=="car")
        this->m_IO = new CIOCar(currentPeriodicFramework);
    else if(vectstr[1]=="poscar")
        this->m_IO = new CIOPoscar(currentPeriodicFramework);
    else if(vectstr[1]=="gjf")
        this->m_IO = new CIOGjf(currentPeriodicFramework);
    else if(vectstr[1]=="cif")
        this->m_IO = new CIOCif(currentPeriodicFramework);
    else if(vectstr[1]=="cell")
        this->m_IO = new CIOCellFile(currentPeriodicFramework);
    else if(vectstr[1]=="xyz")
        this->m_IO = new CIOXyz(currentPeriodicFramework);

    vectstr.clear();

}
void CStructPoolSupported::GeneVARRange(std::vector<GeneVAR>& mht)
{
     this->m_CalcStructPool[0]->GeneVARRange(mht);
}



}//namespace
