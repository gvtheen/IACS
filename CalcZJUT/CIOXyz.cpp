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
#include <string>
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../Util/Point-Vector.h"
#include "../Util/log.hpp"
#include "../IACS.h"
#include "../Util/foreach.h"
#include "CIOXyz.h"

using util::Vector3;
using util::Point3;
using util::Log;
using util::Bitset;

namespace CALCZJUT{

CIOXyz::CIOXyz(CATAZJUT::CConfigurationBase* mth)
:CIOBase(mth)
{
    //ctor
}

CIOXyz::~CIOXyz()
{
    //dtor
}
Bitset CIOXyz::input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE mode)
{
     size_t n1 = m_pPeriodicFramework->atomCount();
     this->input(file);
     size_t n2 = m_pPeriodicFramework->atomCount();
     Bitset res(n2);

     for(size_t i=n1;i<res.size();i++)
        res.set(i);
     return res;
}
void CIOXyz::output(const std::string& file)
{
    assert(m_pPeriodicFramework);
    // check whether this file is exit
    if( m_pPeriodicFramework->dimensionalType() == CATAZJUT::DEFINED::Periodic )
    {
        Log::Error<<"Periodic structure is not saved as xyz file! CIOXyz::output!\n";
        boost::throw_exception(std::runtime_error("Periodic structure is not saved as xyz file! CIOXyz::output!"));
    }
    //check file format
    if(boost::algorithm::contains(file,".xyz")==0){
           Log::Error<<file <<" file is no xyz format! CIOXyz::output"<<std::endl;
           boost::throw_exception(std::runtime_error(file + "file is no no xyz format! Check the file!"));
    }

    std::ofstream out(file,std::ios::app);
    out.setf(std::ios::fixed, std::ios::floatfield);
    out.precision(10);
    if(out.is_open()){
       out<<m_pPeriodicFramework->atomCount()<<std::endl;
       out<<m_pPeriodicFramework->formula()<<std::endl;      //name

       foreach(CATAZJUT::CAtom* atom, m_pPeriodicFramework->atoms())
            out<<atom->Symbol()<<"  "<<atom->position().transpose()<<std::endl;
    }
    out.close();
}
void  CIOXyz::input(std::string file)
{
     assert(m_pPeriodicFramework);

     if (file==""){
        file=m_pPeriodicFramework->formula();
        file=file + ".xyz";
     }

     if(access(file.c_str(),F_OK) != 0 ){
           Log::Error<<file <<" file is no exist! input_CIOXyz!\n";
           boost::throw_exception(std::runtime_error("XYZ file is no exist! Check the file: input_CIOXyz."));
      }
      if(boost::algorithm::contains(file,".xyz")==0){
           Log::Error<<file <<" file is no xyz format! CIOXyz::input"<<std::endl;
           boost::throw_exception(std::runtime_error(file + "file is no no xyz format! Check the code CIOXyz::input!"));
      }
      std::ifstream *in;
      std::string str;
      try{
          in= new std::ifstream(file,std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );
          //1st line
          std::getline(*in,str,'\n');
          size_t atom_Num=std::stod(str);

          //2nd line
          std::getline(*in,str,'\n');

          //From 3rd line  on
          std::vector<std::string> vecStr;
          Vector3 tmpvect;
          for(size_t i=0;i<atom_Num;i++){
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(str==""){
                Log::Error<<i+2 <<"th rows is no exist! input_CIOXyz!\n";
                boost::throw_exception(std::runtime_error("XYZ file is no exist! Check the file: input_CIOXyz."));
             }
             boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
             if(vecStr.size()!=4){
                Log::Error<<i+2 <<"th rows is ERROR! input_CIOXyz!\n";
                boost::throw_exception(std::runtime_error("XYZ file have one ERROR! Check the file: input_CIOXyz."));
             }
             tmpvect<<std::stod(vecStr[1]),std::stod(vecStr[2]),std::stod(vecStr[3]);
             m_pPeriodicFramework->addAtom(vecStr[0],tmpvect);
          }
      }catch(const std::ifstream::failure& e){
          Log::Error<<e.what()<<"Input_CIOPoscar!\n";
          boost::throw_exception(std::runtime_error(e.what()));
      }
      in->close();

      m_pPeriodicFramework->setDimensionalType(CATAZJUT::DEFINED::Molecule);
}

}
