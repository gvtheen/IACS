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
#include "CCrystalPlanes.h"
#include "../Util/log.hpp"

using util::Log;

namespace CATAZJUT{

CCrystalPlanes::CCrystalPlanes()
{
    //ctor
}
CCrystalPlanes::CCrystalPlanes(Eigen::MatrixXd *tempPoints,double m_Value)
{
    m_PointsMat = *tempPoints;
    mDistance_Cutoff = m_Value;
}
CCrystalPlanes::CCrystalPlanes(CCrystalPlanes& tmpCryPlane)
{
   this->SetDistanceCutoff(tmpCryPlane.DistanceCutoff());
   this->SetPointsMat(tmpCryPlane.PointsOfMat());
   this->SetLatticPlane(tmpCryPlane.LatticePlane());
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
Eigen::MatrixXd& CCrystalPlanes::PointsOfMat()
{
    return this->m_PointsMat;
}
void CCrystalPlanes::SetPointsMat(Eigen::MatrixXd& tmpMatPoint)
{

    m_PointsMat = tmpMatPoint;
}
void CCrystalPlanes::SetDistanceCutoff(double mValue)
{
    this->mDistance_Cutoff = mValue;
}
double CCrystalPlanes::DistanceCutoff()
{
    return this->mDistance_Cutoff;
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
    int rowN;
    rowN = m_PointsMat.rows();
    std::vector<bool> matIndex;
    for(int i=0;i<rowN;i++)
        matIndex.push_back(false);

    Eigen::MatrixXd *pointsMat = new (Eigen::MatrixXd)(3,3);

    std::vector<size_t> selectedPoints;
    int Num=0;

    //this->m_Plane = new (std::vector<CPlane*>);
    //this->m_pPointsInIndividualPlanes = new (std::vector<Eigen::MatrixXd*>);
    for(int i=0;i<rowN;i++)
    {
       if(matIndex.at(i)!=true)
       {
          selectedPoints.push_back(i);
          (*pointsMat).row(Num)=m_PointsMat.row(i);
          Num++;
       }
       if( Num == 3 ){
           if(IsOnLine(*pointsMat)==true){
              Num = Num -1;
              continue;
           }
           CPlane *newPlane = new CPlane();
           newPlane->CreatePlane(pointsMat);

           if(CheckIsPlane(*newPlane,selectedPoints,matIndex)==true){
               //collecting the points for newPlane.
               int rowsNum=pointsMat->rows()+selectedPoints.size();
               Eigen::MatrixXd* tempMat = new (Eigen::MatrixXd)(rowsNum,3);
               for(int j=0;j<pointsMat->rows();j++)
                  tempMat->row(j) = pointsMat->row(j);

               int start = pointsMat->rows();
               for(size_t j=0;j<selectedPoints.size();j++){
                  matIndex.at(selectedPoints[j])= true;
                  tempMat->row(start+j) = m_PointsMat.row(selectedPoints[j]);
               }
               m_PointsInIndividualPlanes.push_back(tempMat);
               m_Plane.push_back(newPlane);
               //clear the content of selectedPoints.
               selectedPoints.clear();
           }else
               delete newPlane;
           Num=0;
       }
    }
}
bool CCrystalPlanes::IsOnLine(Eigen::MatrixXd& mMat)
{
   Point3 a,b,c;
   double eps = 0.01;
   a = mMat.col(0);
   b = mMat.col(1);
   c = mMat.col(2);
   if(std::fabs((b - a).cross(c - b).norm())< eps)
      return true;
   else
      return false;
}
bool CCrystalPlanes::CheckIsPlane(CPlane tempPlane,std::vector<size_t> &pIindex,std::vector<bool>& matIndex)
{
     int rowN = m_PointsMat.rows();
     Point3 temP;
     bool res=true;
     double temp_value;
     std::vector<size_t> AddPInthePlane;
     for(int i=0;i<rowN;i++)
     {
         if(matIndex.at(i)==true)
             continue;

               temP = m_PointsMat.row(i);
         temp_value = tempPlane.Distance(temP);
         if(temp_value > 0.5 || std::fabs(temp_value) > this->mDistance_Cutoff )
         {
              res = false;
              break;
         }
         if(std::fabs(temp_value) <= 0.5)
              AddPInthePlane.push_back(i);
     }
     if( res == true )
        pIindex.insert(pIindex.end(),AddPInthePlane.begin(),AddPInthePlane.end());

     return res;
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
}


}
