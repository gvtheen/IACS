#include<iostream>
#include "unistd.h"
#include<fstream>
#include <boost/algorithm/string.hpp>
#include <string>
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../Util/Point-Vector.h"
#include "../Util/log.hpp"
#include "../GACatalyst.h"
#include "../Util/foreach.h"
#include "CIOXyz.h"

using util::Vector3;
using util::Point3;
using util::Log;
using util::Bitset;

namespace CALCZJUT{

CIOXyz::CIOXyz(CATAZJUT::CPeriodicFramework* mth)
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

    if( m_pPeriodicFramework->dimensionalType() == CATAZJUT::DEFINED::Periodic )
    {
        Log::Error<<"Periodic structure is not saved as xyz file! CIOXyz::output!\n";
        boost::throw_exception(std::runtime_error("Periodic structure is not saved as xyz file! CIOXyz::output!"));
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

     if (file=="")
         file="POSCAR";

     if(access(file.c_str(),F_OK) != 0 ){
           Log::Error<<file <<" file is no exist! input_CIOXyz!\n";
           boost::throw_exception(std::runtime_error("XYZ file is no exist! Check the file: input_CIOXyz."));
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
}

}
