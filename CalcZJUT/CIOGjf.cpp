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
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <iostream>
#include <iomanip>
#include "unistd.h"
#include "foreach.h"
#include "CIOGjf.h"
#include "../GaZJUT/GaUtilityFunction.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../CataZJUT/CCartesianCoordinates.h"
#include "../CataZJUT/CFractionCoordinates.h"
#include "../CataZJUT/CConfigurationPrivateData.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../CataZJUT/CAtom.h"
#include "../Util/Point-Vector.h"
#include "../Util/log.hpp"
#include "../IACS.h"

using util::Log;
using util::Point3;
using util::Point3i;
namespace CATAZJUT{
   class CAtom;
}
namespace CALCZJUT{

CIOGjf::CIOGjf(CATAZJUT::CConfigurationBase* mpa)
:CIOBase(mpa)
{
    //ctor
}
void CIOGjf::output(const std::string& fileName)
{
    assert(m_pPeriodicFramework);

    Point3 tp;

    if( this->m_pPeriodicFramework->m_DimensionalType != CATAZJUT::DEFINED::Molecule )
    {
        Log::Error<<"Dimensional Type is error! CIOGjf::output!\n";
        boost::throw_exception(std::runtime_error("Dimensional Type is error! CIOGjf::output!\n"));
    }
    if(this->m_Commandline.size()==0){
        Log::Error<<"Computational setting is not included! CIOGjf::output!\n";
        boost::throw_exception(std::runtime_error("Computational setting is not included! CIOGjf::output!\n"));
    }
    std::string file_Name=fileName;

    std::ofstream out(file_Name.c_str(),std::ios::app);
    out.setf(std::ios_base::right, std::ios_base::adjustfield);
    out.setf(std::ios_base::fixed, std::ios_base::floatfield);
    out.precision(SIGN_DIGIT_NUMBER);

    if(out.is_open())
    {
       for(size_t i=0;i<this->m_Commandline.size();i++)
           out<<this->m_Commandline[i]<<std::endl;

       foreach(CATAZJUT::CAtom* atom_s,m_pPeriodicFramework->atoms()){
          tp=atom_s->position();
          out<<atom_s->Symbol()<<"        ";
          out<<std::setw(SIGN_DIGIT_LENGTH)<<tp(0)<<"    ";
          out<<std::setw(SIGN_DIGIT_LENGTH)<<tp(1)<<"    ";
          out<<std::setw(SIGN_DIGIT_LENGTH)<<tp(2)<<"    ";
          out<<std::endl;
       }
       out<<std::endl;
       out<<std::endl;
    }
    out.close();
}
void CIOGjf::input(std::string file)
{
     if(access(file.c_str(),F_OK) != 0 )
      {
         Log::Error<<"gjf file is no exist! input_CIOGjf!\n";
         boost::throw_exception(std::runtime_error("gjf file is no exist! input_CIOGjf!"));
      }
      std::ifstream *in;
      std::string str;
      try{
         in= new std::ifstream(file.c_str(),std::ifstream::in);
         in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );

         // gaussian command lines
         while(true)
         {
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(str=="")
                break;
             else
                m_Commandline.push_back(str);
         }
         // blank line
         m_Commandline.push_back(" ");

         // comment
         size_t blanklineNum=0;
         while(true)
         {
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             blanklineNum++;
             if(str=="")
                break;
             else
                m_Commandline.push_back(str);
         }
         if(blanklineNum==1){
            Log::Error<<file<<" file format is error! input_CIOGjf!\n";
            boost::throw_exception(std::runtime_error(file + " file format is error!input_CIOGjf!"));
         }

         // blank line
         m_Commandline.push_back(" ");
         // charge multiplicity
         std::getline(*in,str,'\n');
         boost::algorithm::trim(str);
         m_Commandline.push_back(str);

         // via the coordinate of 1st atom, judge the coordinate mode: Cartesian, Internal
         std::getline(*in,str,'\n');
         boost::algorithm::trim(str);
         std::vector<std::string> vecStr;
         boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);

         m_pPeriodicFramework->setDimensionalType(CATAZJUT::DEFINED::Molecule);

         if(vecStr.size()==1){
            //1st atom
            m_pPeriodicFramework->setCoordinateType(CATAZJUT::DEFINED::Internal);
            m_pPeriodicFramework->clear();
            m_pPeriodicFramework->addAtom(vecStr[0]);

            //2nd atom
            std::getline(*in,str,'\n');
            boost::algorithm::trim(str);
            boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
            m_pPeriodicFramework->addAtom(vecStr[0],std::stoi(vecStr[1]),std::stod(vecStr[2]));

            //3rd atom
            std::getline(*in,str,'\n');
            boost::algorithm::trim(str);
            boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
            m_pPeriodicFramework->addAtom(vecStr[0],std::stoi(vecStr[1]),std::stod(vecStr[2]),\
                                          std::stoi(vecStr[3]),std::stod(vecStr[4]));
            Point3i connection;
            Point3  coordinate;
            while(!in->eof())
            {
                std::getline(*in,str,'\n');
                boost::algorithm::trim(str);
                boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
                if(vecStr.size()==7 ||vecStr.size()==8 ){
                  connection<<std::stoi(vecStr[1]),std::stoi(vecStr[3]),std::stoi(vecStr[5]);
                  coordinate<<std::stod(vecStr[2]),std::stod(vecStr[4]),std::stod(vecStr[6]);
                  m_pPeriodicFramework->addAtom(vecStr[0],connection,coordinate);
                }else{
                  Log::Error<<"gjf file has wrong format! input_CIOGjf!\n";
                  boost::throw_exception(std::runtime_error("gjf file has wrong format! input_CIOGjf!\n"));
                }
            }
         }else if(vecStr.size()==4){
            m_pPeriodicFramework->setCoordinateType(CATAZJUT::DEFINED::Cartesian);
            m_pPeriodicFramework->clear();
            Point3 coord;
            coord<<std::stod(vecStr[1]),std::stod(vecStr[2]),std::stod(vecStr[3]);
            m_pPeriodicFramework->addAtom(vecStr[0],coord);

            while(!in->eof())
            {
                std::getline(*in,str,'\n');
                boost::algorithm::trim(str);
                boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);

                if(vecStr.size()==4){
                  coord<<std::stod(vecStr[1]),std::stod(vecStr[2]),std::stod(vecStr[3]);
                  m_pPeriodicFramework->addAtom(vecStr[0],coord);
                }else{
                  Log::Error<<file<<" file has wrong format! input_CIOGjf!\n";
                  boost::throw_exception(std::runtime_error(file+" file has wrong format! input_CIOGjf!\n"));
                }
            }
         }else if(vecStr.size()==5){
            m_pPeriodicFramework->setCoordinateType(CATAZJUT::DEFINED::Cartesian);
            m_pPeriodicFramework->clear();
            //
            Point3 coord;
            coord<<std::stod(vecStr[2]),std::stod(vecStr[3]),std::stod(vecStr[4]);
            //
            m_pPeriodicFramework->addAtom(vecStr[0],coord);
            m_Constranit.push_back(std::stoi(vecStr[1]));

            while(!in->eof())
            {
                std::getline(*in,str,'\n');
                boost::algorithm::trim(str);
                boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);

                if(vecStr.size()==5){
                  coord<<std::stod(vecStr[2]),std::stod(vecStr[3]),std::stod(vecStr[4]);
                  m_Constranit.push_back(std::stoi(vecStr[1]));
                  m_pPeriodicFramework->addAtom(vecStr[0],coord);
                }else{
                  Log::Error<<file<<" file has wrong format!input_CIOGjf!\n";
                  boost::throw_exception(std::runtime_error(file+" file has wrong format!input_CIOGjf!\n"));
                }
            }
         }else{
            Log::Error<<file<<" file has wrong format!input_CIOGjf!\n";
            boost::throw_exception(std::runtime_error(file+" file has wrong format!input_CIOGjf!\n"));
         }

      }catch(const std::ifstream::failure& e){
          Log::Error<<file<<" file has wrong format!input_CIOGjf!\n";
          boost::throw_exception(std::runtime_error(file+"file has wrong format!input_CIOGjf!\n"));
      }
      in->close();
}
Bitset CIOGjf::input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE tmp)
{
     size_t n1 = m_pPeriodicFramework->atomCount();
     this->input(file);
     size_t n2 = m_pPeriodicFramework->atomCount();
     Bitset res(n2);

     for(size_t i=n1;i<res.size();i++)
        res.set(i);
     return res;
}
void CIOGjf::setCommandlines(const std::vector<std::string>& tmpCmdLine)
{
    this->m_Commandline=tmpCmdLine;
}
std::vector<std::string>  CIOGjf::Commandlines() const
{
    return this->m_Commandline;
}

CIOGjf::~CIOGjf()
{
    //dtor
}



}
