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
#include<string>
#include "unistd.h"
#include "stdlib.h"
#include "../Util/foreach.h"
#include "../Util/log.hpp"
#include "../GaZJUT/GaUtilityFunction.h"
#include "../CataZJUT/CAtom.h"
#include "../CataZJUT/CElement.h"
#include "../GaZJUT/GaUtilityFunction.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "CIOCellFile.h"
#include "../CataZJUT/CUnitCell.h"
#include "../CataZJUT/CFractionCoordinates.h"
#include "../Util/Point-Vector.h"

using util::Log;

namespace CALCZJUT{

CIOCellFile::CIOCellFile(CATAZJUT::CConfigurationBase* mpa)
:CIOBase(mpa)
{
    //ctor
}
void CIOCellFile::output(const std::string& file_name)
{
   assert(m_pPeriodicFramework);

   if( m_pPeriodicFramework->m_DimensionalType == CATAZJUT::DEFINED::Molecule )
    {
        Log::Error<<"Dimensional Type is error! CIOCellFile_output"<< std::endl;
        boost::throw_exception(std::runtime_error("Dimensional Type is error!! Check the file: output file."));
    }

    std::string filename;
    if(boost::algorithm::contains(file_name,".cell")>0)
        filename = file_name;
    else
        filename = file_name + ".cell";

    std::ofstream out(filename,std::ios::app);
    out.setf(std::ios::fixed, std::ios::floatfield);
    out.precision(10);
    if(out.is_open()){
        out<<"%BLOCK LATTICE_CART"<<std::endl;
        double scalingfactor = m_pPeriodicFramework->unitcell()->scalingFactor();
        out<<scalingfactor*(m_pPeriodicFramework->unitcell()->avec().transpose())<<std::endl;
        out<<scalingfactor*(m_pPeriodicFramework->unitcell()->bvec().transpose())<<std::endl;
        out<<scalingfactor*(m_pPeriodicFramework->unitcell()->cvec().transpose())<<std::endl;
        out<<"%ENDBLOCK LATTICE_CART"<<std::endl;
        out<<std::endl;
        //Fractional coordinate
        out<<"%BLOCK POSITIONS_FRAC"<<std::endl;
        foreach(CATAZJUT::CAtom* atom, m_pPeriodicFramework->atoms())
        {
          out<<atom->Symbol()<<" "<<m_pPeriodicFramework->Fractioncoordinates()->position(atom->index()).transpose()<<std::endl;
        }
        out<<"%ENDBLOCK POSITIONS_FRAC"<<std::endl;
        //K-point
        out<<std::endl;
        out<<"%BLOCK KPOINTS_LIST"<<std::endl;
        out<<"0.0000000000000000   0.0000000000000000   0.0000000000000000       1.000000000000000"<<std::endl;
        out<<"%ENDBLOCK KPOINTS_LIST"<<std::endl;
        //
        out<<std::endl;
        out<<"FIX_COM : false"<<std::endl;
        out<<"%BLOCK IONIC_CONSTRAINTS"<<std::endl;
        out<<"%ENDBLOCK IONIC_CONSTRAINTS"<<std::endl;

        out<<std::endl;
        out<<"%BLOCK EXTERNAL_EFIELD"<<std::endl;
        out<<"0.0000000000     0.0000000000     0.0000000000"<<std::endl;
        out<<"%ENDBLOCK EXTERNAL_EFIELD"<<std::endl;

        out<<std::endl;
        std::vector<std::pair<std::string,size_t>> tmp = m_pPeriodicFramework->composition();
        std::vector<std::pair<std::string,size_t>>::iterator iter;
        out<<"%BLOCK SPECIES_MASS"<<std::endl;
        for(iter = tmp.begin();iter!= tmp.end(); iter++)
           out<<iter->first<<"   "<< (new CATAZJUT::CElement(iter->first))->exactMass() <<std::endl;

        out<<"%ENDBLOCK SPECIES_MASS"<<std::endl;
        out<<std::endl;

        out<<"%BLOCK SPECIES_POT"<<std::endl;
        if(m_pseudoPotentialFiles.size()==0){
            for(iter = tmp.begin();iter!= tmp.end(); iter++)
                out<<iter->first<<"  "<<iter->first<<"_00PBE.usp" <<std::endl;
        }else{
           std::map<std::string,std::string>::iterator it;
           for(it=m_pseudoPotentialFiles.begin();it!=m_pseudoPotentialFiles.end();it++)
                out<<"      "<<it->first<<"   "<<it->second<<std::endl;
        }
        out<<"%ENDBLOCK SPECIES_POT"<<std::endl;
        out<<std::endl;

        out<<"%BLOCK SPECIES_LCAO_STATES"<<std::endl;
        for(iter = tmp.begin();iter!= tmp.end(); iter++)
              out<<iter->first<<"   "<<iter->second<<std::endl;
        out<<"%ENDBLOCK SPECIES_LCAO_STATES"<<std::endl;
        out<<std::endl;
        out<<std::endl;
        out.close();
    }
}
Bitset CIOCellFile::input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE tmp)
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
void CIOCellFile::input(std::string filename)
{
     //check whether the file is exist.
      if(access(filename.c_str(),F_OK) != 0 ){
           Log::Error<<filename <<" file is no exist! CIOCellFile::input"<<std::endl;
           boost::throw_exception(std::runtime_error(filename + "file is no exist! Check the file!"));
      }
      // check the file format
      if(boost::algorithm::contains(filename,".cell")==0){
           Log::Error<<filename <<" file is no cell format! CIOCellFile::input"<<std::endl;
           boost::throw_exception(std::runtime_error(filename + "file is no no cell format! Check the file!"));
      }
      std::ifstream *in;
      std::string str;
      try{
          in= new std::ifstream(filename.c_str(),std::ifstream::in);
          in->exceptions ( std::ifstream::failbit | std::ifstream::badbit );

          std::getline(*in,str,'\n');
          std::vector<std::string> vecStr;
          Vector3 tmpVet;
          for(size_t i=0;i<3;i++){
            std::getline(*in,str,'\n');
            boost::algorithm::trim(str);
            boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
            tmpVet<<std::stod(vecStr[0]),std::stod(vecStr[1]),std::stod(vecStr[2]);
            m_pPeriodicFramework->unitcell()->setVec(i,tmpVet);
          }
          //set the scalingFactor to 1.0.
          m_pPeriodicFramework->unitcell()->setscalingFactor(1.0);
          //set the coordinate Type
          m_pPeriodicFramework->setCoordinateType(CATAZJUT::DEFINED::Fraction);
          m_pPeriodicFramework->setDimensionalType(CATAZJUT::DEFINED::Periodic);
          std::getline(*in,str,'\n');
          std::getline(*in,str,'\n');
          std::getline(*in,str,'\n');
          while(1){
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str,"%ENDBLOCK POSITIONS_FRAC")>0)
                 break;
             boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
             m_pPeriodicFramework->addAtom(vecStr[0],std::stod(vecStr[1]),std::stod(vecStr[2]),std::stod(vecStr[3]));
          }
          while(1){
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str,"%BLOCK SPECIES_POT")>0)
                 break;
          }

          while(1){
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(str=="")
                continue;
             if(boost::algorithm::contains(str,"%ENDBLOCK SPECIES_POT")>0)
                 break;
             boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
             if(vecStr.size()==2)
                m_pseudoPotentialFiles[vecStr[0]]=vecStr[1];
             else{
                Log::Error<<filename <<"Pseudopotential file is no error! CIOCellFile::input"<<std::endl;
                boost::throw_exception(std::runtime_error(filename + "Pseudopotential file is no error! CIOCellFile::input!"));
             }
          }
      }catch(const std::ifstream::failure& e){
          Log::Error<< e.what() <<"input_CIOCellFile"<<std::endl;
          boost::throw_exception(std::runtime_error(filename + "file is no exist! Check the file: Error_information.txt."));
      }
      in->close();
}
CIOCellFile::~CIOCellFile()
{
    m_pseudoPotentialFiles.clear();
}

}
