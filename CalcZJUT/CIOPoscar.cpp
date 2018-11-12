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
#include<iostream>
#include "unistd.h"
#include<fstream>
#include <boost/algorithm/string.hpp>
#include<string>
#include "../Util/Point-Vector.h"
#include "../Util/log.hpp"
#include "../GaZJUT/GaUtilityFunction.h"
#include "CIOPoscar.h"
#include "../CataZJUT/CUnitCell.h"
#include "CParameter.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../CataZJUT/CCartesianCoordinates.h"
#include "../CataZJUT/CFractionCoordinates.h"
#include "../CataZJUT/CConfigurationPrivateData.h"
using util::Log;
using util::Vector3;
using util::Point3;
using CATAZJUT::CFractionCoordinates;
using CATAZJUT::CCartesianCoordinates;
namespace CALCZJUT{

CIOPoscar::CIOPoscar(CATAZJUT::CPeriodicFramework* mth)
:CIOBase(mth)
{

}
/*
Format of POSCAR file:
1st line:   name of system
2nd line:   scaling factor
3rd line:   vector a
4th line:   vector b
5th line:   vector c
6th line:   element name
7th line:   atomic number
8th line:   coordinate mode   D: direct fraction coordination; C or K:  Carteancoordinate

*/
Bitset CIOPoscar::input(std::string file,CParameter::SIMULATION_MODE temp)
{
     assert(m_pPeriodicFramework);

     size_t n1 = m_pPeriodicFramework->atomCount();
     this->input(file);
     size_t n2 = m_pPeriodicFramework->atomCount();
     Bitset res(n2);

     for(size_t i=n1;i<res.size();i++)
        res.set(i);
     return res;
}
void CIOPoscar::input(std::string file)
{
     assert(m_pPeriodicFramework);

     if (file=="")
         file="POSCAR";

     if(access(file.c_str(),F_OK) != 0 )
      {
           Log::Error<<"POSCAR file is no exist! input_CIOPoscar!\n";
           boost::throw_exception(std::runtime_error("POSCAR file is no exist! Check the file: Error_information.txt."));
      }
      std::ifstream *in;
      std::string str;
      try{
          in= new std::ifstream(file,std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );
          //name of compound
          std::getline(*in,str,'\n');
          boost::algorithm::trim(str);
          m_pPeriodicFramework->m_pData->name = str;    //name of compound
          //// scaling factor
          std::getline(*in,str,'\n');
          boost::algorithm::trim(str);
          m_pPeriodicFramework->unitcell()->setscalingFactor(std::stod(str));
          // lattice vector
          Vector3 tmpvect;
          std::vector<std::string> vecStr;
          for(size_t i=0;i<3;i++)
          {
              std::getline(*in,str,'\n');
              boost::algorithm::trim(str);
              boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
              tmpvect<<std::stod(vecStr[0]),std::stod(vecStr[1]),std::stod(vecStr[2]);
              m_pPeriodicFramework->unitcell()->setVec(i,tmpvect);
              vecStr.clear();
          }
          //setting dimensional type.
          m_pPeriodicFramework->setDimensionalType(CATAZJUT::DEFINED::Periodic);

          //atomic name
          std::vector<std::string> AtomicName;
          std::getline(*in,str,'\n');
          boost::algorithm::trim(str);
          boost::algorithm::split(AtomicName,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
          //atomic number
          std::vector<size_t>  AtomicNum;
          std::getline(*in,str,'\n');
          boost::algorithm::trim(str);

          boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
          for(size_t i=0;i<vecStr.size();i++)
             AtomicNum.push_back(std::stoi(vecStr[i]));

          if(AtomicName.size()!=AtomicNum.size())
          {
               Log::Error<<"Atomic setting of POSCAR is error! input_CIOPoscar!\n";
               boost::throw_exception(std::runtime_error("Atomic setting is error! Check the file: Error_information.txt."));
          }
          // coordinate mode:  cartesian fraction
          std::getline(*in,str,'\n');
          boost::algorithm::trim(str);
          if(str[0]=='c' || str[0]=='C' || str[0]=='k' || str[0]=='K')
             m_pPeriodicFramework->setCoordinateType(CATAZJUT::DEFINED::Cartesian);
          else if(str[0]=='D' || str[0]=='d' )
             m_pPeriodicFramework->setCoordinateType(CATAZJUT::DEFINED::Fraction);

          size_t atom_Number=0;
          size_t index;
          while(!in->eof())
          {
              std::getline(*in,str,'\n');
              boost::algorithm::trim(str);
              if(str=="")
                 continue;
              boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
              if(vecStr.size()==6){
                 atom_Number++;
                 tmpvect<<std::stod(vecStr[0]),std::stod(vecStr[1]),std::stod(vecStr[2]);
                 if(m_pPeriodicFramework->coordinateType()==CATAZJUT::DEFINED::Cartesian){
                      // Oblique coordinate cartesian coordinate
                      tmpvect =m_pPeriodicFramework->unitcell()->NormilizedBravaisMatrix()*tmpvect;

                      index = this->atomicIndex(AtomicNum,atom_Number);
                      m_pPeriodicFramework->addAtom(AtomicName[index],tmpvect);
                 }else{
                      index = this->atomicIndex(AtomicNum,atom_Number);
                      m_pPeriodicFramework->addAtom(AtomicName[index],tmpvect);
                 }

              }else{
                 Log::Error<<"Coordinate in POSCAR file is error! input_CIOPoscar!\n";
                 boost::throw_exception(std::runtime_error("Coordinate in POSCAR file is error! input_CIOPoscar!\n"));
              }

          }
      }catch(const std::ifstream::failure& e){
          Log::Error<<e.what()<<"Input_CIOPoscar!\n";
          exit(-1);
      }
      in->close();
}
size_t CIOPoscar::atomicIndex(std::vector<size_t>& mht,size_t& atomNum)
{
    std::vector<size_t> tmp;
    size_t num=0,res;

    for(size_t i=0;i<mht.size();i++)
    {
        num=num + mht[i];
        tmp.push_back(num);
    }
    for(size_t i=0;i<tmp.size()-1;i++){
         if(i==0 && atomNum<=tmp[0]){
            res = 0;
            break;
         }else if( i > 0 && ( atomNum > tmp[i] && atomNum <=tmp[i+1])){
            res = i;
            break;
         }
    }
    return res;
}
void CIOPoscar::output(const std::string& file)
{
    assert(m_pPeriodicFramework);

    if( m_pPeriodicFramework->m_DimensionalType == CATAZJUT::DEFINED::Molecule )
    {
        Log::Error<<"Dimensional Type is error! CIOPoscar::output!\n";
        boost::throw_exception(std::runtime_error("Dimensional Type is error!! Check the file: CIOPoscar::output."));
    }
    m_pPeriodicFramework->sortAtomsViaElements();

    std::ofstream out(file,std::ios::app);
    out.setf(std::ios::fixed, std::ios::floatfield);
    out.precision(10);
    if(out.is_open())
    {
       out<<m_pPeriodicFramework->formula()<<std::endl;      //name
       out<<m_pPeriodicFramework->unitcell()->scalingFactor()<<std::endl;   // scaling factor
       out<<m_pPeriodicFramework->unitcell()->avec().transpose()<<std::endl;
       out<<m_pPeriodicFramework->unitcell()->bvec().transpose()<<std::endl;
       out<<m_pPeriodicFramework->unitcell()->cvec().transpose()<<std::endl;
       // atomic list of compound
       std::vector<std::pair<std::string,size_t>> compos = m_pPeriodicFramework->composition();

       for( size_t i=0;i<compos.size();i++)
           out<<compos[i].first<<"   ";
       out<<std::endl;
       for( size_t i=0;i<compos.size();i++)
           out<<compos[i].second<<"   ";
       out<<std::endl;
       compos.clear();
       //coordinate type: cart direct
       if(m_pPeriodicFramework->coordinateType()!= CATAZJUT::DEFINED::Fraction)
       {
           out<<"Direct"<<std::endl;
           CFractionCoordinates* coord = m_pPeriodicFramework->Fractioncoordinates();
           for(size_t i=0;coord->size();i++)
               out<<"  "<<coord->position(i).transpose()<<"  T   T   T"<<std::endl;
           out<<std::endl;
           out<<std::endl;
       }else{
           out<<"Cartesian"<<std::endl;
           CCartesianCoordinates* coord = m_pPeriodicFramework->coordinates();
           Point3 tmp;
           for(size_t i=0;coord->size();i++)
           {
               tmp=m_pPeriodicFramework->unitcell()->NormilizedBravaisMatrix()*(coord->position(i));
               out<<"  "<<tmp.transpose()<<"  T   T   T"<<std::endl;
           }
           out<<std::endl;
           out<<std::endl;
       }
       out.close();
    }
}

CIOPoscar::~CIOPoscar()
{
    //dtor
}



}//namespace
