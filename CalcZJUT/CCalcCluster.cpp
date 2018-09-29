#include <boost/algorithm/string.hpp>
#include <string>
#include <math.h>
#include "../CataZJUT/CPeriodicFramework.h"
#include "CCalcCluster.h"
#include "../CataZJUT/CConfigurationBase.h"
#include "../CataZJUT/CElement.h"
#include "../Util/log.hpp"
#include "../Util/Bitset.h"
#include "../Util/Point-Vector.h"
#include "../CataZJUT/Constant.h"

using util::Log;

namespace CALCZJUT{

CCalcCluster::CCalcCluster(CParameter* mPara)
:CCalcModeStruct(mPara)
{
    CATAZJUT::CConfigurationBase* temp = new CATAZJUT::CConfigurationBase(mPara);

    this->m_pPeriodicFramework = new CATAZJUT::CPeriodicFramework(*temp);
     // In this situation, this->m_pPeriodicFramework->m_pUnitCell is empty pointer NULL;
    delete temp;

}

CCalcCluster::~CCalcCluster()
{
    //dtor
}
void CCalcCluster::setGeneValueToStruct(const std::vector<double>& realValueOfgene)
{

}
std::vector<double>*  CCalcCluster::getGeneValuefromStruct()const
{

}
std::vector<GENEVAR>* CCalcCluster::GeneVarRange()
{

}
void CCalcCluster::Initialization(const std::string& mth)
{
       // initialize from chemical
       // chemical formula = C10 H10 O12
/*
    e.g.:   mth = C10 N S10        // delimiter = blank " "
*/
    std::vector<std::pair<std::string,size_t>> chemicalFormula;
    //dealwith chemical formula by blank
    std::string str=mth;
    std::vector<std::string> vecStr;
    boost::algorithm::trim(str);
    boost::algorithm::split(vecStr,str,boost::algorithm::is_any_of(" "),boost::algorithm::token_compress_on);

    std::pair<std::string,size_t> ChemicalStr;
    for(size_t i=0;i<vecStr.size();i++){
        ChemicalStr.first  = boost::algorithm::trim_copy_if(vecStr[i], boost::algorithm::is_digit());
        //check whether the label is OK.
        if( CATAZJUT::CElement::isValidSymbol(ChemicalStr.first) == false ){
           Log::Error<<vecStr[i]<< " is error. Initialization_CCalcCluster!\n";
           boost::throw_exception(std::runtime_error(vecStr[i] + " is error. Initialization_CCalcCluster!\n"));//ERROR TREATMENT;
        }
        str = boost::algorithm::trim_copy_if(vecStr[i], boost::algorithm::is_alpha());
        if ( str=="" )
            ChemicalStr.second = 1;
        else
            ChemicalStr.second = std::stoi(str);
        chemicalFormula.push_back(ChemicalStr);
    }


}
void CCalcCluster::Initialization(char* mth)
{      // initialize from chemical formula
     std::string tmp(mth);
     this->Initialization(tmp);
}
void CCalcCluster::Initialization(const CATAZJUT::CConfigurationBase& mth)  // initialize from exit
{

}
void CCalcCluster::RandomBuildFromChemicalFormula(std::vector<std::pair<std::string,size_t>> mht)
{
     std::vector<CATAZJUT::CElement*> res;
     size_t atom_Sum=0;
     for(size_t i=0;i<mht.size();i++){
        res.push_back(new CATAZJUT::CElement(mht[i].first));
        atom_Sum = atom_Sum + mht[i].second;
     }
     util::Vector3 polar_coord, basic_coord;
     util::Point3 coordinate;
     double cosTheta,sinTheta,cosPhi,sinPhi;
     size_t cluster_type = ClusterType(res);
     if(cluster_type==1){        //pure metal
         basic_coord<< std::pow(atom_Sum,1.0/3),180.0,360.0;
         CATAZJUT::CPeriodicFramework* new_Struct = new CATAZJUT::CPeriodicFramework(this->m_pParameter);
         for(size_t i=0;i<mht.size();i++)
            for(size_t j=0;j<mht[i].second;j++){
                    // using polar coordinate to predict it
                polar_coord = util::Vector3::Random().normalized() * basic_coord.transpose();
                   // transfer polar coordinate to cartesian coordinate.
                cosTheta= std::cos(polar_coord(1)*CATAZJUT::constants::DegreesToRadians);
                sinTheta= std::sin(polar_coord(1)*CATAZJUT::constants::DegreesToRadians);
                cosPhi = std::cos(polar_coord(2)*CATAZJUT::constants::DegreesToRadians);
                sinPhi = std::sin(polar_coord(2)*CATAZJUT::constants::DegreesToRadians);
                coordinate<<polar_coord(0)*cosTheta*cosPhi,polar_coord(0)*cosTheta*sinPhi,polar_coord(0)*sinTheta;
                new_Struct->addAtom(mht[i].first,coordinate);
            }
     }else if(cluster_type==2){  //pure nonmetal

     }else if(cluster_type==3){   //pure mixed nonmetal

     }
}
size_t CCalcCluster::ClusterType(std::vector<CATAZJUT::CElement*>& mht)
{
   // size_t M_Pointer=0,NoM_Pointer=0, H_Pointer=0;
    util::Bitset Element_Bit(2);
     //set all bit to be 1;
    Element_Bit.set();
    for(size_t i=0;i<mht.size();i++)
    {
       if( mht[i]->isMetal() )    Element_Bit.set(0,false);
       if( mht[i]->isNonmetal() ) Element_Bit.set(1,false);
    }
    if(Element_Bit.to_ulong()== BOOST_BINARY(10))
        return 1;
    if(Element_Bit.to_ulong()== BOOST_BINARY(01))
        return 2;
    if(Element_Bit.to_ulong()== BOOST_BINARY(00))
        return 3;
    return 0;
}


}
