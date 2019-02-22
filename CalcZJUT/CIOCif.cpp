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
#include <iostream>
#include "unistd.h"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <string>
#include "../Util/Point-Vector.h"
#include "CIOCif.h"
#include "../CataZJUT/CPeriodicFramework.h"
#include "../CataZJUT/CCartesianCoordinates.h"
#include "../CataZJUT/CFractionCoordinates.h"
#include "../CataZJUT/CUnitCell.h"
#include "../Util/foreach.h"
#include "../Util/log.hpp"

using util::Log;
using util::Point3;
namespace CALCZJUT{

CIOCif::CIOCif(CATAZJUT::CPeriodicFramework* mpa)
:CIOBase(mpa)
{
    //ctor
}

CIOCif::~CIOCif()
{
    //dtor
}
void CIOCif::output(const std::string& file_name)
{
    assert(m_pPeriodicFramework);

    if( this->m_pPeriodicFramework->dimensionalType() == CATAZJUT::DEFINED::Molecule )
    {
        Log::Error<<"Dimensional Type is error! CIOCif::output!\n";
        boost::throw_exception(std::runtime_error("Dimensional Type is error! CIOCif::output!\n"));
    }

    std::string file_Name;
    if(boost::algorithm::contains(file_name,".cif")>0)
        file_Name = file_name;
    else
        file_Name = file_name + ".cif";

    std::ofstream out(file_Name.c_str(),std::ios::app);
    out.setf(std::ios::fixed, std::ios::floatfield);
    out.precision(10);
    if(out.is_open()){
       out<<"data_"<<m_pPeriodicFramework->formula()<<std::endl;
       out<<"_audit_creation_method     "<<"IACS_out"<<std::endl;
       out<<"_cell_length_a             "<<m_pPeriodicFramework->unitcell()->a()<<std::endl;
       out<<"_cell_length_b             "<<m_pPeriodicFramework->unitcell()->b()<<std::endl;
       out<<"_cell_length_c             "<<m_pPeriodicFramework->unitcell()->c()<<std::endl;
       out<<"_cell_angle_alpha          "<<m_pPeriodicFramework->unitcell()->alfa()<<std::endl;
       out<<"_cell_angle_beta           "<<m_pPeriodicFramework->unitcell()->beta()<<std::endl;
       out<<"_cell_angle_gamma          "<<m_pPeriodicFramework->unitcell()->gama()<<std::endl;
       out<<"_symmetry_space_group_H-M     'P1'"<<std::endl;
       out<<"_symmetry_Int_Tables_number   '1'"<<std::endl;
       out<<"_symmetry_cell_setting        'triclinic'"<<std::endl;
       out<<std::endl;

       out<<"loop_"<<std::endl;
       out<<"_symmetry_equiv_pos_as_xyz"<<std::endl;
       out<<"x,y,z"<<std::endl;
       out<<std::endl;

       out<<"loop_"<<std::endl;
       out<<"_atom_site_label"<<std::endl;
       out<<"_atom_site_type_symbol"<<std::endl;
       out<<"_atom_site_occupancy"<<std::endl;
       out<<"_atom_site_fract_x"<<std::endl;
       out<<"_atom_site_fract_y"<<std::endl;
       out<<"_atom_site_fract_z"<<std::endl;
       out<<"_atom_site_U_iso_or_equiv"<<std::endl;

       // Here donot delete temp pointer.
       CATAZJUT::CCartesianCoordinates* temp = m_pPeriodicFramework->coordinates();
       CATAZJUT::CFractionCoordinates *res = temp->toFractionCoordinates(m_pPeriodicFramework);
       std::vector<std::pair<std::string,size_t>> tempComposition = m_pPeriodicFramework->composition();
       std::map<std::string,size_t> label_Map;
       for(size_t i=0;i<tempComposition.size();i++)
           label_Map[tempComposition[i].first] = 1;

       Point3 pos;
       foreach(CATAZJUT::CAtom* atom, m_pPeriodicFramework->atoms()){
          pos = (*res)[atom->index()];
          out<<atom->Symbol()<<label_Map[atom->Symbol()]<<"    ";
          out<<atom->Symbol()<<"     1     "<<pos[0]<<"     "<<pos[1]<<"     "<<pos[2];
          out<<"        0"<<std::endl;
          label_Map[atom->Symbol()]=label_Map[atom->Symbol()]+1;
       }
       delete res;
       out<<std::endl;
       out<<std::endl;
       out.close();
    }

}
void CIOCif::input(std::string filename)
{
      assert(m_pPeriodicFramework);
      //check whether the file is exist.
      if(access(filename.c_str(),F_OK) != 0 ){
           Log::Error<<filename <<" file is no exist! CIOCellFile::input"<<std::endl;
           boost::throw_exception(std::runtime_error(filename + "file is no exist! Check the file CIOCif::input!"));
      }
      // check the file format
      if(boost::algorithm::contains(filename,".cif")==0){
           Log::Error<<filename <<" file is no cif format! CIOCif::input"<<std::endl;
           boost::throw_exception(std::runtime_error(filename + "file is no no cell format! Check the file CIOCif::input!"));
      }

      std::ifstream *in;
      std::string str;
      std::vector<std::string> vecStr;
      std::vector<double> celldata;
      //clear all
      this->m_pPeriodicFramework->clear();

      try{
          in= new std::ifstream(filename.c_str(),std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );
          while(1){
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str,"_cell_length_a")>0)
                 break;
          }
          for(size_t i=0;i<6;i++){
            boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
            celldata.push_back(std::stod(vecStr[1]));
            std::getline(*in,str,'\n');
            boost::algorithm::trim(str);
          }
          m_pPeriodicFramework->unitcell()->fromCellPara(celldata[0],celldata[1],celldata[2],
                                                         celldata[3],celldata[4],celldata[5],'p');
          while(1){
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str,"loop_")>0){
                std::getline(*in,str,'\n');
                if(boost::algorithm::contains(str,"_atom_site_label")>0)
                    break;
             }
          }
          for(size_t i=0;i<6;i++)
             std::getline(*in,str,'\n');
          //set coordination type
          m_pPeriodicFramework->setCoordinateType(CATAZJUT::DEFINED::Fraction);
          while(1){
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(str=="")
                break;
             boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
             if(vecStr.size()!=7){
                Log::Error<<filename <<" file format is error! CIOCif::input"<<std::endl;
                boost::throw_exception(std::runtime_error(filename + "file format is error! Check the file CIOCif::input!"));
             }
             m_pPeriodicFramework->addAtom(vecStr[1],std::stod(vecStr[3]),
                                           std::stod(vecStr[4]),std::stod(vecStr[5]));
          }
          m_pPeriodicFramework->perceiveBonds();

      }catch(const std::ifstream::failure& e){
          Log::Error<< e.what() <<"CIOCif::input"<<std::endl;
          boost::throw_exception(std::runtime_error(filename + "file is no exist! Check the file: CIOCif::input."));
      }
      in->close();
}
Bitset CIOCif::input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE mode)
{
     assert(m_pPeriodicFramework);

     //check whether the file is exist.
     if(access(file.c_str(),F_OK) != 0 ){
           Log::Error<<file <<" file is no exist! CIOCellFile::input"<<std::endl;
           boost::throw_exception(std::runtime_error(file + "file is no exist! Check the file CIOCellFile::input!"));
     }
      // check the file format
     if(boost::algorithm::contains(file,".cell")==0){
           Log::Error<<file<<" file is no cell format! CIOCellFile::input"<<std::endl;
           boost::throw_exception(std::runtime_error(file + "file is no no cell format! Check the file CIOCellFile::input!"));
     }
     size_t n1 = m_pPeriodicFramework->atomCount();
     this->input(file);
     size_t n2 = m_pPeriodicFramework->atomCount();
     Bitset res(n2);

     for(size_t i=n1;i<res.size();i++)
        res.set(i);
     return res;
}

}
