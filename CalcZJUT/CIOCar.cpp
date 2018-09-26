#include <ctime>
#include <fstream>
#include <iostream>
#include "unistd.h"
#include <boost/algorithm/string.hpp>
#include "GaUtilityFunction.h"
#include "CUnitCell.h"
#include "Point-Vector.h"
#include "foreach.h"
#include "CIOCar.h"
#include "CAtom.h"
#include "CUnitCell.h"
#include "CElement.h"
#include "CPeriodicFramework.h"
#include "CElement.h"

using GAZJUT::ERROR_OUTPUT;
using CATAZJUT::Vector3;

namespace CALCZJUT{

CIOCar::CIOCar(CATAZJUT::CPeriodicFramework* mpa)
:Cios(mpa)
{
    //ctor
}

CIOCar::~CIOCar()
{
    //dtor
}
void CIOCar::output(std::string& file)
{
    std::ofstream out(file,std::ios::app);
    out.setf(std::ios::fixed, std::ios::floatfield);
    out.precision(10);
    if(out.is_open())
    {
        out<<"!BIOSYM archive 3"<<std::endl;
        if(m_pPeriodicFramework->dimensionalType()==CATAZJUT::DEFINED::Periodic)
            out<<"PBC=ON"<<std::endl;
        else
            out<<"PBC=OFF"<<std::endl;
        out<<"GACatalyst Program generate structure of " << m_pPeriodicFramework->formula()<<std::endl; //name
        std::time_t result = std::time(nullptr);
        out<<"!DATE "<<std::asctime(std::localtime(&result))<<std::endl;
        if(m_pPeriodicFramework->dimensionalType()==CATAZJUT::DEFINED::Periodic)
            out<<"PBC  "<<m_pPeriodicFramework->unitcell()->a() <<"  "<<   \
                          m_pPeriodicFramework->unitcell()->b() <<"  "<<   \
                          m_pPeriodicFramework->unitcell()->c() <<"  "<<   \
                       m_pPeriodicFramework->unitcell()->alfa() <<"  "<<   \
                       m_pPeriodicFramework->unitcell()->beta() <<"  "<<   \
                       m_pPeriodicFramework->unitcell()->gama() <<"(P1)"<<std::endl;
         foreach(CATAZJUT::CAtom* atom, m_pPeriodicFramework->atoms()){
            out<<atom->element().symbol()<<atom->index()<<"  "<<atom->position().transpose()
               <<"  XXXX 1    XX"<<atom->element().symbol()<<"    0.0"<<std::endl;
         }
         out<<"end"<<std::endl;
         out<<"end"<<std::endl;
         out.close();
    }
}
void CIOCar::input(std::string  file)
{
     bool isPBC = false;
     if(access(file.c_str(),F_OK) != 0 )
      {
           ERROR_OUTPUT(file + " file is no exist!","input","CIOCar");
           boost::throw_exception(std::runtime_error(file + "  file is no exist! Check the file: Error_information.txt."));
      }
      std::ifstream *in;
      std::string str;
      Vector3 tmpvect;
      std::vector<std::string> vecStr;
      try{
          in= new std::ifstream(file,std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );
          std::getline(*in,str,'\n');      // 1st line

          std::getline(*in,str,'\n');      // 2nd line
          boost::algorithm::trim(str);
          boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of("="),boost::algorithm::token_compress_on);
          if(vecStr.size()!=2 || vecStr[0]!="PBC"){
            ERROR_OUTPUT(file + " file has NO PBC setting!","input","CIOCar");
            boost::throw_exception(std::runtime_error(file + " file has NO PBC setting! Check the file: Error_information.txt."));
          }else{
             boost::algorithm::trim(vecStr[1]);
             if(vecStr[1]=="ON")
                isPBC = true;
             else if(vecStr[1]=="OFF")
                isPBC = false;
             else{
                ERROR_OUTPUT(file + " file has NO PBC setting!","input","CIOCar");
                boost::throw_exception(std::runtime_error(file + " file has NO PBC setting! Check the file: Error_information.txt."));
             }
           }

           std::getline(*in,str,'\n');   //3rd line
           std::getline(*in,str,'\n');   //4th line
           if(isPBC == true){
              m_pPeriodicFramework->setDimensionalType(CATAZJUT::DEFINED::Periodic);
              std::getline(*in,str,'\n');   //5th line
              boost::algorithm::trim(str);
              boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
              if(vecStr.size()==8){
                 m_pPeriodicFramework->unitcell()->fromCellPara(std::stod(vecStr[1]),std::stod(vecStr[2]),std::stod(vecStr[3]),\
                                                               std::stod(vecStr[4]),std::stod(vecStr[5]),std::stod(vecStr[6]),vecStr[7][1]);
              }else{
                 ERROR_OUTPUT(file + "PBC setting is error!","input","CIOCar");
                 boost::throw_exception(std::runtime_error(file + " PBC setting is error! Check the file: Error_information.txt."));
              }
           }else
              m_pPeriodicFramework->setDimensionalType(CATAZJUT::DEFINED::Molecule);
           while(1)
           {
               std::getline(*in,str,'\n');   //ith line
               boost::algorithm::trim(str);
               if(str=="end") break;
               boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
               tmpvect<<std::stod(vecStr[1]),std::stod(vecStr[2]),std::stod(vecStr[3]);
               if(isPBC == true){
                  tmpvect=m_pPeriodicFramework->unitcell()->NormilizedBravaisMatrix()*tmpvect;
               }
               m_pPeriodicFramework->addAtom(vecStr[7],tmpvect);
           }
      }catch(const std::ifstream::failure& e){
          ERROR_OUTPUT(e.what(),"Input","CIOCar");
          exit(-1);
      }
      in->close();
}
Bitset CIOCar::input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE mode)
{
     size_t n1 = m_pPeriodicFramework->atomCount();
     this->input(file);
     size_t n2 = m_pPeriodicFramework->atomCount();
     Bitset res(n2);

     for(size_t i=n1;i<res.size();i++)
        res.set(i);
     return res;
}



}
