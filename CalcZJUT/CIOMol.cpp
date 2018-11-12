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
#include<fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <string>
#include "unistd.h"
#include "stdlib.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "CIOMol.h"
#include "foreach.h"
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CBond.h"
#include "../CataZJUT/CConfigurationPrivateData.h"
#include "../GaZJUT/GaUtilityFunction.h"
#include "../Util/log.hpp"

using util::Log;

namespace CATAZJUT{
   class CAtom;
   class CBond;
}

namespace CALCZJUT{


CIOMol::CIOMol(CATAZJUT::CPeriodicFramework* mpa)
:CIOBase(mpa)
{
    //ctor
}
void CIOMol::output(const std::string& fileName)
{
    assert(m_pPeriodicFramework);

    if( this->m_pPeriodicFramework->m_DimensionalType != CATAZJUT::DEFINED::Molecule )
    {
        Log::Error<<"Dimensional Type is error! CIOMol::output!\n";
        boost::throw_exception(std::runtime_error("Dimensional Type is error! CIOMol::output!\n"));
    }

    std::string file_Name;
    if(boost::algorithm::contains(fileName,".mol")>0)
        file_Name = fileName;
    else
        file_Name = fileName + ".mol";

    std::ofstream out(file_Name.c_str(),std::ios::out);
    out.setf(std::ios::fixed, std::ios::floatfield);
    out.precision(10);
    if(out.is_open())
    {
        out<<m_pPeriodicFramework->formula()<<std::endl;      //name
        out<<"GACatalyst program ( the author: Gvtheen )"<<std::endl;
        out<<std::endl;
        out<<m_pPeriodicFramework->atomCount()<<"  "<<m_pPeriodicFramework->atomCount() \
           <<"  0  0  0  0  0  0  0  0999 V2000"<<std::endl;
        foreach(CATAZJUT::CAtom* atom,m_pPeriodicFramework->atoms())
        {
           out<<atom->position().transpose()<<" "<<atom->Symbol()<<" " \
              <<"  0  0  0  0  0  0  0  0  0  0  0  0"<<std::endl;
        }
        foreach(CATAZJUT::CBond* bond,m_pPeriodicFramework->bonds())
        {
           out<<bond->atom1()<<"  "<<bond->atom2()<<"  2  0 "<<std::endl;
        }
        out<<"M  END"<<std::endl;
    }
     out.close();

}
Bitset CIOMol::input(std::string file,CParameter::SIMULATION_MODE mode)
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
void  CIOMol::input(std::string filename)
{
     assert(m_pPeriodicFramework);

     if(access(filename.c_str(),F_OK) != 0 )
      {
           Log::Error<<filename <<" file is no exist! input_CIOMol!\n";
           boost::throw_exception(std::runtime_error(filename +" file is no exist! input_CIOMol!\n"));
      }
      std::ifstream *in;
      std::string str;
      try{
          in= new std::ifstream(filename.c_str(),std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );

          std::getline(*in,str,'\n');
          boost::algorithm::trim(str);
          m_pPeriodicFramework->m_pData->name = str;    //name of compound

          std::getline(*in,str,'\n');   //comment line
          std::getline(*in,str,'\n');   // blank line

          std::getline(*in,str,'\n');   //count line
          boost::algorithm::trim(str);
          std::vector<std::string> vecStr;
          boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);

          size_t atomicNum=std::stoi(vecStr[0]);
          for(size_t i=0;i<atomicNum;i++)
          {
              std::getline(*in,str,'\n');   //count line
              boost::algorithm::trim(str);
              boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
              boost::algorithm::trim(vecStr[3]);
              m_pPeriodicFramework->addAtom(vecStr[3],std::stod(vecStr[0]),std::stod(vecStr[1]),std::stod(vecStr[2]));
          }
      }catch(const std::ifstream::failure& e){
          Log::Error<<e.what() <<"input_CIOMol!\n";
      }
      in->close();
}
CIOMol::~CIOMol()
{
    //dtor
}


}
