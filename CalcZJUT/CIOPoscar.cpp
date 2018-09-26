#include<iostream>
#include "unistd.h"
#include<fstream>
#include <boost/algorithm/string.hpp>
#include<string>
#include "Point-Vector.h"
#include "GaUtilityFunction.h"
#include "CIOPoscar.h"
#include "CUnitCell.h"
#include "CParameter.h"
#include "CPeriodicFramework.h"
#include "CCartesianCoordinates.h"
#include "CFractionCoordinates.h"
#include "CConfigurationPrivateData.h"
using GAZJUT::ERROR_OUTPUT;
using CATAZJUT::Vector3;
using CATAZJUT::Point3;
using CATAZJUT::CFractionCoordinates;
using CATAZJUT::CCartesianCoordinates;
namespace CALCZJUT{

CIOPoscar::CIOPoscar(CATAZJUT::CPeriodicFramework* mth)
:Cios(mth)
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
           ERROR_OUTPUT("POSCAR file is no exist!","input","CIOPoscar");
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
          for(size_t i=0;i<3;i++)
          {
              std::vector<std::string> vecStr;
              std::getline(*in,str,'\n');
              boost::algorithm::trim(str);
              boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
              tmpvect<<std::stod(vecStr[0]),std::stod(vecStr[1]),std::stod(vecStr[2]);
              m_pPeriodicFramework->unitcell()->setVec(i,tmpvect);
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
          std::vector<std::string> vecStr;
          boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
          for(size_t i=0;i<vecStr.size();i++)
             AtomicNum.push_back(std::stoi(vecStr[i]));

          if(AtomicName.size()!=AtomicNum.size())
          {
               ERROR_OUTPUT("Atomic setting of POSCAR is error!","input", "CIOPoscar");
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
                      // °ÑOblique coordinate ×ª»»Îª cartesian coordinate
                      tmpvect =m_pPeriodicFramework->unitcell()->NormilizedBravaisMatrix()*tmpvect;

                      index = this->atomicIndex(AtomicNum,atom_Number);
                      m_pPeriodicFramework->addAtom(AtomicName[index],tmpvect);
                 }else{
                      index = this->atomicIndex(AtomicNum,atom_Number);
                      m_pPeriodicFramework->addAtom(AtomicName[index],tmpvect);
                 }

              }else
                 ERROR_OUTPUT("Coordinate POSCAR is error!","input", "CIOPoscar");
          }
      }catch(const std::ifstream::failure& e){
          ERROR_OUTPUT(e.what(),"Input","CIOPoscar");
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
void CIOPoscar::output(std::string& file)
{
    assert(m_pPeriodicFramework);

    if( m_pPeriodicFramework->m_DimensionalType == CATAZJUT::DEFINED::Molecule )
    {
        ERROR_OUTPUT("Dimensional Type is error!","CIOPoscar::output");
        boost::throw_exception(std::runtime_error("Dimensional Type is error!! Check the file: Error_information.txt."));
    }
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
       std::vector<std::pair<std::string,size_t>>* compos = m_pPeriodicFramework->composition();

       for( size_t i=0;i<compos->size();i++)
           out<<(*compos)[i].first<<"   ";
       out<<std::endl;
       for( size_t i=0;i<compos->size();i++)
           out<<(*compos)[i].second<<"   ";
       out<<std::endl;
       delete compos;
       //coordinate type: cart direct
       if(m_pPeriodicFramework->coordinateType()!= CATAZJUT::DEFINED::Fraction)
       {
           out<<"Direct"<<std::endl;
           CFractionCoordinates* coord = m_pPeriodicFramework->Fractioncoordinates();
           for(size_t i=0;coord->size();i++)
              out<<"  "<<coord->position(i).transpose()<<"  T   T   T"<<std::endl;
           out<<std::endl<<std::endl;
       }else{
           out<<"Cartesian"<<std::endl;
           CCartesianCoordinates* coord = m_pPeriodicFramework->coordinates();
           Point3 tmp;
           for(size_t i=0;coord->size();i++)
           {
               tmp=m_pPeriodicFramework->unitcell()->NormilizedBravaisMatrix()*(coord->position(i));
               out<<"  "<<tmp.transpose()<<"  T   T   T"<<std::endl;
           }
           out<<std::endl<<std::endl;
       }
       out.close();
    }
}

CIOPoscar::~CIOPoscar()
{
    //dtor
}



}//namespace
