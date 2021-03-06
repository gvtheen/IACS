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
#include<cmath>
#include "../IACS.h"
#include "CAtom.h"
#include "CCrystalPlanes.h"
#include "CConfigurationBase.h"
#include "../Util/log.hpp"
#include "../Util/Bitset.h"
#include "../Util/Point-Vector.h"

using util::Point3;
using util::Log;

namespace CATAZJUT{

CCrystalPlanes::CCrystalPlanes()
{
    //ctor
}
CCrystalPlanes::CCrystalPlanes(std::vector<size_t>& index_Vect,CConfigurationBase* tmpConfigurationBase)
{
    assert(index_Vect.size()>=3);
    assert(tmpConfigurationBase);

    this->m_pCConfigurationBase=tmpConfigurationBase;

    this->m_indexInConfig.assign(index_Vect.begin(),index_Vect.end());

    m_PointsMat= new Eigen::MatrixXd(index_Vect.size(),3);
    Point3 tmpPoint;
    for(size_t i=0;i<m_indexInConfig.size();i++){
         tmpPoint=this->m_pCConfigurationBase->atom(m_indexInConfig[i])->position();

         (*m_PointsMat)(i,0)=tmpPoint(0);
         (*m_PointsMat)(i,1)=tmpPoint(1);
         (*m_PointsMat)(i,2)=tmpPoint(2);
    }


}
CCrystalPlanes::CCrystalPlanes(CCrystalPlanes& tmpCryPlane)
{
   #ifdef DEBUG
        Log::Debug<<"CCrystalPlanes::CCrystalPlanes(CCrystalPlanes& tmpCryPlane)" << std::endl;
   #endif // DEBU
   this->SetDistanceCutoff(tmpCryPlane.DistanceCutoff());
   this->SetPointsMat(tmpCryPlane.PointsOfMat());
   this->SetLatticPlane(tmpCryPlane.LatticePlane());
   this->setConfigurationBase(tmpCryPlane.ConfigurationBase());
   #ifdef DEBUG
        Log::Debug<<"2-CCrystalPlanes::CCrystalPlanes(CCrystalPlanes& tmpCryPlane)" << std::endl;
   #endif
}
CCrystalPlanes* CCrystalPlanes::Clone()
{
    return new CCrystalPlanes(*this);
}
void CCrystalPlanes::SetLatticPlane(std::vector<CPlane*> &mpt)
{
    if(this->m_Plane.size()!=0){
       for(size_t i=0;i<this->m_Plane.size();i++)
           delete m_Plane[i];
       this->m_Plane.clear();
    }

    for(size_t i=0;i<mpt.size();i++)
        this->m_Plane.push_back(new CPlane(*(mpt[i])));
}
Eigen::MatrixXd* CCrystalPlanes::PointsOfMat()
{
    return this->m_PointsMat;
}
void CCrystalPlanes::SetPointsMat(Eigen::MatrixXd* tmpMatPoint)
{
    size_t num=tmpMatPoint->rows();

    m_PointsMat = new Eigen::MatrixXd(num,3);
    for(size_t i=0;i<num;i++)
    {
       (*m_PointsMat)(i,0)= (*tmpMatPoint)(i,0);
       (*m_PointsMat)(i,1)= (*tmpMatPoint)(i,1);
       (*m_PointsMat)(i,2)= (*tmpMatPoint)(i,2);
    }
}
void CCrystalPlanes::SetDistanceCutoff(double mValue)
{
    this->mDistance_Cutoff = mValue;
}
double CCrystalPlanes::DistanceCutoff()
{
    return this->mDistance_Cutoff;
}
CConfigurationBase* CCrystalPlanes::ConfigurationBase()
{
    return this->m_pCConfigurationBase;
}
void CCrystalPlanes::setConfigurationBase( CConfigurationBase* mbf)
{
   this->m_pCConfigurationBase=mbf;
}

std::vector<CPlane*>& CCrystalPlanes::LatticePlane()
{
    return this->m_Plane;
}
CPlane* CCrystalPlanes::operator[](size_t index)
{
    if( index >= crystalPlaneNum()){
       Log::Error<<"Crystal plane index is out of range! CCrystalPlanes::operator[]!\n";
       boost::throw_exception(std::runtime_error("Crystal plane index is out of range! CCrystalPlanes::operator[]!"));
    }
    return this->m_Plane[index];
}
size_t CCrystalPlanes::crystalPlaneNum()
{
    return this->m_Plane.size();
}
void CCrystalPlanes::CreateCrystalPlane()
{
    size_t rowN = m_PointsMat->rows();
    util::Bitset matIndex(rowN);
    matIndex.set(false);

    Eigen::MatrixXd *pointsMat = new (Eigen::MatrixXd)(3,3);
    #ifdef DEBUG
        Log::Debug<<"CCrystalPlanes::CreateCrystalPlane() ROWN:"<<rowN << std::endl;
    #endif // DEBU
    std::vector<size_t> selectedPoints;
    //this->m_Plane = new (std::vector<CPlane*>);
    //this->m_pPointsInIndividualPlanes = new (std::vector<Eigen::MatrixXd*>);
    for(size_t i=0;i<rowN;i++)
      for(size_t j=i+1;j<rowN;j++)
       for(size_t k=j+1;k<rowN;k++){

              if(IsInclude3PointsInExitPlane(i,j,k))
                 continue;

              if(IsAdjecentPointer(i,j,k)==false)
                 continue;

              matIndex.set(false);

              selectedPoints.clear();
              selectedPoints.push_back(i);
              selectedPoints.push_back(j);
              selectedPoints.push_back(k);
              matIndex.set(i,true);
              matIndex.set(j,true);
              matIndex.set(k,true);
              (*pointsMat).row(0)=m_PointsMat->row(i);
              (*pointsMat).row(1)=m_PointsMat->row(j);
              (*pointsMat).row(2)=m_PointsMat->row(k);

           if(IsOnLine(*pointsMat)==true)
              continue;

           CPlane *newPlane = new CPlane();
           newPlane->CreatePlane(pointsMat);

           if(CheckIsCrystalPlane(*newPlane,selectedPoints,matIndex)==true){
               //collecting the points for newPlane.
               int rowsNum=selectedPoints.size();
               Eigen::MatrixXd* tempMat = new (Eigen::MatrixXd)(rowsNum,3);
               //for(int j=0;j<pointsMat->rows();j++)
                //  tempMat->row(j) = pointsMat->row(j);

              // int start = pointsMat->rows();
               for(size_t j=0;j<selectedPoints.size();j++){
                  //matIndex.at(selectedPoints[j])= true;
                  tempMat->row(j) = m_PointsMat->row(selectedPoints[j]);
               }

               bool isNewPlane=true;
               for(size_t i=0;i<m_Plane.size();i++)
                   if( (*m_Plane[i]) == (*newPlane) ){
                        isNewPlane=false;
                        break;
                   }
               if(isNewPlane){
                  m_Plane.push_back(newPlane);
                  m_PointsInIndividualPlanes.push_back(tempMat);
               }else{
                  delete newPlane;
                  delete tempMat;
               }
               //clear the content of selectedPoints.
               selectedPoints.clear();
           }else
               delete newPlane;
     }
}
bool CCrystalPlanes::IsInclude3PointsInExitPlane(size_t m,size_t n,size_t h)
{
    Eigen::MatrixXd* data=new (Eigen::MatrixXd)(3,3);
    (*data).row(0)=(*m_PointsMat).row(m);
    (*data).row(1)=(*m_PointsMat).row(n);
    (*data).row(2)=(*m_PointsMat).row(h);
//    #ifdef DEBUG
//        Log::Debug<<"CCrystalPlanes::IsInclude3PointsInExitPlane" << std::endl;
//    #endif // DEBU
    util::Bitset bol(3);
    bol.set();
    size_t pos=0;

    for(size_t i=0;i<this->m_PointsInIndividualPlanes.size();i++){
        bol.set();
        pos=0;
        for(size_t j=0;j<this->m_PointsInIndividualPlanes[i]->rows();j++){
            if((*data)(pos,0)==(*m_PointsInIndividualPlanes[i])(j,0)&&
               (*data)(pos,1)==(*m_PointsInIndividualPlanes[i])(j,1)&&
               (*data)(pos,2)==(*m_PointsInIndividualPlanes[i])(j,2)){
                bol.set(pos,false);
                pos = bol.find_next(pos);
                if(pos==util::Bitset::npos){
                    delete data;
                    return true;
                }
            }
        }
    }
//    #ifdef DEBUG
//        Log::Debug<<"2-CCrystalPlanes::IsInclude3PointsInExitPlane" << std::endl;
//    #endif // DEBU
    delete data;
    return false;
}
bool CCrystalPlanes::IsOnLine(Eigen::MatrixXd& mMat)
{
   Point3 a,b,c;
   double eps = 0.5;
   a = mMat.col(0);
   b = mMat.col(1);
   c = mMat.col(2);
   if(std::fabs((b - a).cross(c - b).norm())< eps)
      return true;
   else
      return false;
}
bool CCrystalPlanes::IsAdjecentPointer(size_t i,size_t j,size_t k)
{
   CAtom  *atom_i,*atom_j,*atom_k;
   atom_i=m_pCConfigurationBase->atom(m_indexInConfig[i]);
   atom_j=m_pCConfigurationBase->atom(m_indexInConfig[j]);
   atom_k=m_pCConfigurationBase->atom(m_indexInConfig[k]);

   if(atom_i->isBondedTo(atom_j)){
      if(atom_j->isBondedTo(atom_k)||atom_i->isBondedTo(atom_k))
         return true;
      else
         return false;
   }else if(atom_j->isBondedTo(atom_k)&& atom_i->isBondedTo(atom_k))
      return true;
   else
      return false;
}
bool CCrystalPlanes::CheckIsCrystalPlane(CPlane& tempPlane,std::vector<size_t> &pIindex,util::Bitset& matIndex)
{
     int rowN = m_PointsMat->rows();
     Point3 temP;
     bool IsCrystalPlane=true;
     double temp_value;
     std::vector<size_t> AddPInthePlane;
     double convergence=0.1,Last_Distance=0;

     for(size_t i=0;i<rowN;i++)
     {
         if(matIndex[i]==1)
             continue;

               temP = m_PointsMat->row(i);
         temp_value = tempPlane.Distance(temP);

         if(std::fabs(temp_value) >=convergence )
         {
            if(Last_Distance==0)
               Last_Distance=temp_value;
            else if(Last_Distance*temp_value<0){
               IsCrystalPlane=false;
               break;
            }
         }else
            AddPInthePlane.push_back(i);

     }

     if( IsCrystalPlane == true )
        pIindex.insert(pIindex.end(),AddPInthePlane.begin(),AddPInthePlane.end());
     else
        AddPInthePlane.clear();

     return IsCrystalPlane;
}
Point3 CCrystalPlanes::CartesianCoordinateAtGene(size_t crystal_Plane_index,double height,
                                                 double R_radio,double thea_Radian)
{
     if(crystal_Plane_index>=this->m_Plane.size()){
         Log::Error<<"Crystal plane index is out of range! CCrystalPlanes::cartesianCoordinateAtGene!\n";
         boost::throw_exception(std::runtime_error("Crystal plane index is out of range! CCrystalPlanes::cartesianCoordinateAtGene!"));
     }
     return m_Plane[crystal_Plane_index]->PointInCircleFromGene(height,R_radio,thea_Radian);
}
void CCrystalPlanes::RemoveRow(Eigen::MatrixXd& matrix, size_t rowToRemove)
{
    size_t numRows = matrix.rows()-1;
    size_t numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,\
                                                                  numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void CCrystalPlanes::RemoveColumn(Eigen::MatrixXd& matrix, size_t colToRemove)
{
    size_t numRows = matrix.rows();
    size_t numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,\
                                                                  numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}
void CCrystalPlanes::outputCrystalPlane()
{
    assert(m_Plane.size()>0);

    Log::Info<<"*******Crystal Plane Equations output*******"<<std::endl;
    Log::Info<<" Number of crystal plane: "<<m_Plane.size()<<std::endl;
    for(size_t i=0;i<this->m_Plane.size();i++)
        this->m_Plane[i]->outputPlane();

    Log::Info<<"*******End Crystal Plane Equations*******"<<std::endl;
}
CCrystalPlanes::~CCrystalPlanes()
{
//        if(m_PointsMat.rows()!=0)
//            m_PointsMat.;
        if(m_Plane.size()!=0)
        {
            for(size_t i=0;i<m_Plane.size();i++)
                delete m_Plane[i];
            m_Plane.clear();
        }
        if(m_PointsInIndividualPlanes.size()!=0)
        {
            for(size_t i=0;i<m_PointsInIndividualPlanes.size();i++)
                delete m_PointsInIndividualPlanes[i];
            m_PointsInIndividualPlanes.clear();
        }
        if(m_PointsMat!=nullptr)
            delete m_PointsMat;
}


}
