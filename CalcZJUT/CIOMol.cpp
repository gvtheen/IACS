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
#include "CConfigurationBase.h"
#include "CIOMol.h"
#include "foreach.h"
#include "CAtom.h"
#include "CBond.h"
#include "CConfigurationPrivateData.h"
#include "GaUtilityFunction.h"
namespace CATAZJUT{
   class CAtom;
   class CBond;
}
using GAZJUT::ERROR_OUTPUT;

namespace CALCZJUT{


CIOMol::CIOMol(CATAZJUT::CConfigurationBase* mpa)
:CIOBase(mpa)
{
    //ctor
}
void CIOMol::output(std::string& file_Name)
{
    if( m_pConfigurationBase->m_DimensionalType != CATAZJUT::DEFINED::Molecule )
    {
        ERROR_OUTPUT("Dimensional Type is error!","CIOMol::output");
        boost::throw_exception(std::runtime_error("Dimensional Type is error!! Check the file: Error_information.txt."));
    }
    file_Name = file_Name + ".mol";
    std::ofstream out(file_Name,std::ios::app);
    out.setf(std::ios::fixed, std::ios::floatfield);
    out.precision(10);
    if(out.is_open())
    {
        out<<m_pConfigurationBase->formula()<<std::endl;      //name
        out<<"GACatalyst program ( the author: Gvtheen )"<<std::endl;
        out<<std::endl;
        out<<m_pConfigurationBase->atomCount()<<"  "<<m_pConfigurationBase->atomCount() \
           <<"  0  0  0  0  0  0  0  0999 V2000"<<std::endl;
        foreach(CATAZJUT::CAtom* atom,m_pConfigurationBase->atoms())
        {
           out<<atom->position().transpose()<<" "<<atom->Symbol()<<" " \
              <<"  0  0  0  0  0  0  0  0  0  0  0  0"<<std::endl;
        }
        foreach(CATAZJUT::CBond* bond,m_pConfigurationBase->bonds())
        {
           out<<bond->atom1()<<"  "<<bond->atom2()<<"  2  0 "<<std::endl;
        }
        out<<"M  END"<<std::endl;
    }
     out.close();

}
Bitset CIOMol::input(std::string file,CParameter::SIMULATION_MODE mode)
{
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
     if(access(filename.c_str(),F_OK) != 0 )
      {
           ERROR_OUTPUT(filename + " file is no exist!","input","CIOMol");
           boost::throw_exception(std::runtime_error(filename +" file is no exist! Check the file: Error_information.txt."));
      }
      std::ifstream *in;
      std::string str;
      try{
          in= new std::ifstream(filename.c_str(),std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );

          std::getline(*in,str,'\n');
          boost::algorithm::trim(str);
          m_pConfigurationBase->m_pData->name = str;    //name of compound

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
              m_pConfigurationBase->addAtom(vecStr[3],std::stod(vecStr[0]),std::stod(vecStr[1]),std::stod(vecStr[2]));
          }
      }catch(const std::ifstream::failure& e){
          ERROR_OUTPUT(e.what(),"input","CIOMol");
          exit(-1);
      }
      in->close();
}
CIOMol::~CIOMol()
{
    //dtor
}


}
