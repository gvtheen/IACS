#include<fstream>
#include <boost/algorithm/string.hpp>
#include<string>
#include<iostream>
#include "unistd.h"
#include "Point-Vector.h"
#include "CIOGjf.h"
#include "GaUtilityFunction.h"
#include "CPeriodicFramework.h"
#include "CCartesianCoordinates.h"
#include "CFractionCoordinates.h"
#include "CConfigurationPrivateData.h"
#include "CConfigurationBase.h"
#include "Point-Vector.h"

using GAZJUT::ERROR_OUTPUT;
using CATAZJUT::Point3;
using CATAZJUT::Point3i;

namespace CALCZJUT{

CIOGjf::CIOGjf(CCATAZJUT::CPeriodicFramework* mpa)
:Cios(mpa)
{
    //ctor
}
void CIOGjf::output(std::string& file)
{

}
void CIOGjf::input(std::string file)
{
     if(access(file.c_str(),F_OK) != 0 )
      {
         ERROR_OUTPUT("gjf file is no exist!","input","CIOGjf");
         boost::throw_exception(std::runtime_error("gjf file is no exist! Check the file: Error_information.txt."));
      }
      std::ifstream *in;
      std::string str;
      try{
         in= new std::ifstream("POSCAR",std::ifstream::in);
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
                  ERROR_OUTPUT("gjf file has wrong format!","input","CIOGjf");
                  boost::throw_exception(std::runtime_error("gjf file has wrong format! Check the file: Error_information.txt."));
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
                  ERROR_OUTPUT("gjf file has wrong format!","input","CIOGjf");
                  boost::throw_exception(std::runtime_error("gjf file has wrong format! Check the file: Error_information.txt."));
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
                  ERROR_OUTPUT("gjf file has wrong format!","input","CIOGjf");
                  boost::throw_exception(std::runtime_error("gjf file has wrong format! Check the file: Error_information.txt."));
                }
            }
         }else{
            ERROR_OUTPUT("gjf file has wrong format!","input","CIOGjf");
            boost::throw_exception(std::runtime_error("gjf file has wrong format! Check the file: Error_information.txt."));
         }

      }catch(const std::ifstream::failure& e){
          ERROR_OUTPUT("gjf file has wrong format!","input","CIOGjf");
          boost::throw_exception(std::runtime_error("gjf file has wrong format! Check the file: Error_information.txt."));
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
CIOGjf::~CIOGjf()
{
    //dtor
}



}
