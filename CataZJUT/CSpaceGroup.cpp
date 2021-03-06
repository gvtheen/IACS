#include "CSpaceGroup.h"
#include "../Util/foreach.h"
#include "CAtom.h"
#include "CConfigurationBase.h"
#include "CUnitCell.h"
#include "../Spacegroup/spglib.h"
#include "../Util/log.hpp"
#include "../Util/Point-Vector.h"

using util::Point3;
using util::Log;

namespace CATAZJUT{

CSpaceGroup::CSpaceGroup(CConfigurationBase* mbf)
:m_pCConfigurationBase(mbf)
{
    //ctor
    m_pSpaceGroupName = new char[21];
}

CSpaceGroup::~CSpaceGroup()
{
    delete []m_pSpaceGroupName;
}
char* CSpaceGroup::GetSpaceGroup()
{
 // double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
  double (*lattice)[3];
  double (*position)[3];
//  {
//      {0, 0, 0},
//      {0.5, 0.5, 0.5},
//      {0.3, 0.3, 0},
//      {0.7, 0.7, 0},
//      {0.2, 0.8, 0.5},
//      {0.8, 0.2, 0.5},
//  };
  //int types[] = {1, 1, 2, 2, 2, 2};
  int num_spg, num_atom;
  num_atom = this->m_pCConfigurationBase->atomCount();

  lattice = (double(*)[3])(new double[3*3]);
  /*
Transtion matrix A
ax   bx   cx
ay   by   cy
az   bz   cz
*/
  Eigen::Matrix<double, 3, 3> Latt = m_pCConfigurationBase->unitcell()->MatrixOfBravaisLattice();
  Latt.transpose();
  for(size_t i=0;i<3;i++)
      for(size_t j=0;j<3;j++)
        lattice[i][j]=Latt(i,j);

  int *types = new int[num_atom];
  size_t index=0;
  position =(double(*)[3])(new double[3*num_atom]);
  std::vector<std::pair<std::string,size_t>> composition = m_pCConfigurationBase->composition();

  Point3 pos;
  foreach(CAtom* atom, this->m_pCConfigurationBase->atoms()){
    for(size_t i=0;i<composition.size();i++)
       if(atom->Symbol() == composition[i].first ){
           types[index] = i+1;
           break;
       }
    pos = atom->position();
    position[index][0]=pos[0];
    position[index][1]=pos[1];
    position[index][2]=pos[2];
  }

  num_spg = spg_get_international(m_pSpaceGroupName, lattice, position, types, num_atom, 1e-5);

  delete position;
  delete types;
  delete lattice;

  if ( num_spg > 0 )
    return m_pSpaceGroupName;
   else{
    Log::Error<<"Space group could not be found."<<std::endl;
    return nullptr;
  }

}
bool CSpaceGroup::Find_primitive()
{
   double (*lattice)[3];
   double (*position)[3];

   int  num_atom;
   num_atom = this->m_pCConfigurationBase->atomCount();

   lattice = (double(*)[3])(new double[3*3]);

   Eigen::Matrix<double, 3, 3> Latt = m_pCConfigurationBase->unitcell()->MatrixOfBravaisLattice();
   Latt.transpose();
   for(size_t i=0;i<3;i++)
      for(size_t j=0;j<3;j++)
        lattice[i][j]=Latt(i,j);

   int *types = new int[num_atom];
   size_t index=0;
   position =(double(*)[3])(new double[3*num_atom]);
   std::vector<std::pair<std::string,size_t>> composition = m_pCConfigurationBase->composition();

   Point3 pos;
   foreach(CAtom* atom, this->m_pCConfigurationBase->atoms()){
    for(size_t i=0;i<composition.size();i++)
       if(atom->Symbol() == composition[i].first ){
           types[index] = i+1;
           break;
       }
    pos = atom->position();
    position[index][0]=pos[0];
    position[index][1]=pos[1];
    position[index][2]=pos[2];
  }
  int num_primitive_atom;
  double symprec = 1e-5;
  num_primitive_atom = spg_find_primitive(lattice, position, types, num_atom, symprec);

  if (num_primitive_atom != 0){
      m_pCConfigurationBase->clear();
      for(size_t i=0;i<3;i++){
         pos<<lattice[i][0],lattice[i][1],lattice[i][2];
         m_pCConfigurationBase->unitcell()->setVec(i,pos);
      }
      m_pCConfigurationBase->unitcell()->setscalingFactor(1.0);

      for(int i=0;i<num_primitive_atom;i++){
         pos<<position[i][0],position[i][1],position[i][2];
         m_pCConfigurationBase->addAtom(composition[types[i]-1].first,pos);
      }
      m_pCConfigurationBase->perceiveBonds();
      return true;
  }else
      return false;

}


}
