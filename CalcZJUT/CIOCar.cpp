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
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "unistd.h"
#include <boost/algorithm/string.hpp>
#include "../GaZJUT/GaUtilityFunction.h"
#include "../CataZJUT/CUnitCell.h"
#include "../Util/Point-Vector.h"
#include "../Util/foreach.h"
#include "../Util/log.hpp"
#include "CIOCar.h"
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CElement.h"
#include "../CataZJUT/CConfigurationBase.h"

using util::Vector3;
using util::Log;
namespace CALCZJUT{

CIOCar::CIOCar(CATAZJUT::CConfigurationBase* mpa)
:CIOBase(mpa)
{
    //ctor
}

CIOCar::~CIOCar()
{
    //dtor
}
void CIOCar::output(const std::string& file)
{
    Point3 tp;

    std::ofstream out(file,std::ios::out);

    //out.setf(std::ios_base::right, std::ios_base::adjustfield);
    out.setf(std::ios_base::fixed, std::ios_base::floatfield);
    out.precision(9);
    if(out.is_open())
    {
        out<<"!BIOSYM archive 3"<<std::endl;
        if(m_pPeriodicFramework->dimensionalType()==CATAZJUT::DEFINED::Periodic)
            out<<"PBC=ON"<<std::endl;
        else
            out<<"PBC=OFF"<<std::endl;
        out<<"IACS Program generate structure of " << m_pPeriodicFramework->formula()<<std::endl; //name
        std::time_t result = std::time(nullptr);
        out<<"!DATE "<<std::asctime(std::localtime(&result));
        if(m_pPeriodicFramework->dimensionalType()==CATAZJUT::DEFINED::Periodic){
            out<<"PBC";
            out<<std::setw(10)<<m_pPeriodicFramework->unitcell()->a(); /*<<"  "<<   */
            out<<std::setw(10)<<m_pPeriodicFramework->unitcell()->b();  /*<<"  "<<   */
            out<<std::setw(10)<<m_pPeriodicFramework->unitcell()->c();  /*<<"  "<<   */\
            out<<std::setw(10)<<m_pPeriodicFramework->unitcell()->alfa(); /*<<"  "<<   */
            out<<std::setw(10)<<m_pPeriodicFramework->unitcell()->beta(); /*<<"  "<<   */
            out<<std::setw(10)<<m_pPeriodicFramework->unitcell()->gama()<<"(P1)"<<std::endl;
        }
         std::string tmp;
         foreach(CATAZJUT::CAtom* atom, m_pPeriodicFramework->atoms()){
            tmp=std::to_string(atom->index()+1);
            out.setf(std::ios_base::left, std::ios_base::adjustfield);
            tmp=atom->element().symbol()+tmp;
            out<<std::setw(5)<<tmp<<"  ";
            out.setf(std::ios_base::right, std::ios_base::adjustfield);
            tp=atom->position();
            out<<std::setw(13)<<tp(0)<<"  ";
            out<<std::setw(13)<<tp(1)<<"  ";
            out<<std::setw(13)<<tp(2)<<" ";
            out<<"XXXX 1      xx      "<<atom->element().symbol()<<"  0.000"<<std::endl;
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
           Log::Error<<file <<" file is no exist! input_CIOCar!\n";
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
            Log::Error<< file<<" file has NO PBC setting! input_CIOCar"<<std::endl;
            boost::throw_exception(std::runtime_error(file + " file has NO PBC setting! Check the file: Error_information.txt."));
          }else{
             boost::algorithm::trim(vecStr[1]);
             if(vecStr[1]=="ON")
                isPBC = true;
             else if(vecStr[1]=="OFF")
                isPBC = false;
             else{
                Log::Error<<file <<" file has NO PBC setting!_input_CIOCar!"<<std::endl;
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
                 Log::Error<<file<< " PBC setting is error! input_CIOCar!"<<std::endl;
                 boost::throw_exception(std::runtime_error(file + " PBC setting is error! Check the file: Error_information.txt."));
              }
           }else
              m_pPeriodicFramework->setDimensionalType(CATAZJUT::DEFINED::Molecule);
           bool endbol=false;
           while(!in->eof())
           {
               std::getline(*in,str,'\n');   //ith line
               boost::algorithm::trim(str);
               if(str=="end") {
                  if(endbol==true)
                      break;
                  else{
                      endbol=true;
                      continue;
                  }
               }else
                  endbol=false;

               boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
               tmpvect<<std::stod(vecStr[1]),std::stod(vecStr[2]),std::stod(vecStr[3]);
               if(isPBC == true){
                  tmpvect=m_pPeriodicFramework->unitcell()->NormilizedBravaisMatrix()*tmpvect;
               }
               m_pPeriodicFramework->addAtom(vecStr[7],tmpvect);
           }
      }catch(const std::ifstream::failure& e){
          Log::Error<<e.what()<<" Input_CIOCar!"<<std::endl;
          boost::throw_exception(std::runtime_error(e.what()));
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
