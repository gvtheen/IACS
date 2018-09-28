#include<fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include<string>
#include "unistd.h"
#include "stdlib.h"
#include "../Util/foreach.h"
#include "../Util/log.hpp"
#include "GaUtilityFunction.h"
#include "CAtom.h"
#include "GaUtilityFunction.h"
#include "CPeriodicFramework.h"
#include "CIOCellFile.h"
#include "CUnitCell.h"
#include "CFractionCoordinates.h"
#include "Point-Vector.h"

using util::Log;

namespace CALCZJUT{

CIOCellFile::CIOCellFile(CATAZJUT::CPeriodicFramework* mpa)
:Cios(mpa)
{
    //ctor
}
void CIOCellFile::output(std::string& filename)
{
   if( m_pPeriodicFramework->m_DimensionalType == CATAZJUT::DEFINED::Molecule )
    {
        Log::OutputToFile<<"ERROR: " <<"Dimensional Type is error! CIOCellFile_output"<< std::endl;
        boost::throw_exception(std::runtime_error("Dimensional Type is error!! Check the file: output file."));
    }
    filename = filename + ".cell";
    std::ofstream out(filename,std::ios::app);
    out.setf(std::ios::fixed, std::ios::floatfield);
    out.precision(10);

    out<<"%BLOCK LATTICE_CART"<<std::endl;
    double scalingfactor = m_pPeriodicFramework->unitcell()->scalingFactor();
    out<<scalingfactor*(m_pPeriodicFramework->unitcell()->avec().transpose())<<std::endl;
    out<<scalingfactor*(m_pPeriodicFramework->unitcell()->bvec().transpose())<<std::endl;
    out<<scalingfactor*(m_pPeriodicFramework->unitcell()->cvec().transpose())<<std::endl;
    out<<"%ENDBLOCK LATTICE_CART"<<std::endl;
    out<<std::endl;
    //Fractional coordinate
    out<<"%BLOCK POSITIONS_FRAC"<<std::endl;
    foreach(CAtom* atom, m_pPeriodicFramework->atoms())
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
    std::vector<std::pair<std::string,size_t>>* tmp = m_pPeriodicFramework->composition();
    std::vector<std::pair<std::string,size_t>>::iterator iter;
    out<<"%BLOCK SPECIES_MASS"<<std::endl;
    for(iter = tmp->begin();iter!= tmp->end(); iter++)
    {
       out<<iter->first<<"   "<< (new CElement(iter->first))->exactMass() <<std::endl;
    }
    out<<"%ENDBLOCK SPECIES_MASS"<<std::endl;
    out<<std::endl;

    out<<"%BLOCK SPECIES_POT"<<std::endl;
    for(iter = tmp->begin();iter!= tmp->end(); iter++)
    {
        out<<iter->first<<"  "<<iter->first<<"_00PBE.usp" <<std::endl;
    }
    out<<"%ENDBLOCK SPECIES_POT"<<std::endl;
    out<<std::endl;

    out<<"%BLOCK SPECIES_LCAO_STATES"<<std::endl;
    for(iter = tmp->begin();iter!= tmp->end(); iter++)
          out<<iter->first<<"   "<<iter->second<<std::endl;
    out<<"%ENDBLOCK SPECIES_LCAO_STATES"<<std::endl;
    out<<std::endl;

    out.close();
}
Bitset CIOCellFile::input(std::string file,CALCZJUT::CParameter::SIMULATION_MODE tmp)
{
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
     if(access(filename.c_str(),F_OK) != 0 )
      {
           Log.OutputToFile<<filename <<" file is no exist! input_CIOPoscar"<<std::endl;
           boost::throw_exception(std::runtime_error(filename + "file is no exist! Check the file!"));
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
          m_pPeriodicFramework->unitcell()->setscalingFactor(1.0);
          std::getline(*in,str,'\n');
          std::getline(*in,str,'\n');
          std::getline(*in,str,'\n');
          while(1){
             std::getline(*in,str,'\n');
             boost::algorithm::trim(str);
             if(boost::algorithm::contains(str,"%ENDBLOCK")>0)
                 break;
             boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);
             m_pPeriodicFramework->addAtom(vecStr[0],std::stod(vecStr[1]),std::stod(vecStr[2]),std::stod(vecStr[3]));
          }
      }catch(const std::ifstream::failure& e){
          Log::OutputToFile<<"ERROR: "<< e.what() <<"input_CIOCellFile"<<std::endl;
          boost::throw_exception(std::runtime_error(filename + "file is no exist! Check the file: Error_information.txt."));
      }
      in->close();
}
CIOCellFile::~CIOCellFile()
{
    //dtor
}

}
