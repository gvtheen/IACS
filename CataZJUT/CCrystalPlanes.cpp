#include "CCrystalPlanes.h"
#include<cmath>
namespace CATAZJUT{

CCrystalPlanes::CCrystalPlanes()
{
    //ctor
}
CCrystalPlanes::CCrystalPlanes(Eigen::MatrixXd *tempPoints,double m_Value)
{
    int rowN,colN;
    rowN = tempPoints->rows();
    colN = tempPoints->cols();
     m_pPointsMat = new (Eigen::MatrixXd)(rowN,colN);
    *m_pPointsMat = *tempPoints;
    mDistance_Cutoff = m_Value;

    m_pPlane = nullptr;
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
void CCrystalPlanes::SetLatticPlane(std::vector<CPlane*> *mpt)
{
    if(this->m_pPlane!=nullptr)
        delete this->m_pPlane;
    this->m_pPlane = new (std::vector<CPlane*>);

    for(size_t i=0;i<mpt->size();i++)
        this->m_pPlane->push_back(new CPlane(*(mpt->at(i))));
}
Eigen::MatrixXd* CCrystalPlanes::PointsOfMat()
{
    return this->m_pPointsMat;
}
void CCrystalPlanes::SetPointsMat(Eigen::MatrixXd* tmpMatPoint)
{
    if(m_pPointsMat == nullptr)
    {
        int rowN,colN;
        rowN=tmpMatPoint->rows();
        colN=tmpMatPoint->cols();
        m_pPointsMat = new (Eigen::MatrixXd)(rowN,colN);
    }
    *m_pPointsMat = *tmpMatPoint;
}
void CCrystalPlanes::SetDistanceCutoff(double mValue)
{
    this->mDistance_Cutoff = mValue;
}
double CCrystalPlanes::DistanceCutoff()
{
    return this->mDistance_Cutoff;
}
std::vector<CPlane*>* CCrystalPlanes::LatticePlane()
{
    return this->m_pPlane;
}
void CCrystalPlanes::CreateCrystalPlane()
{
    int rowN;
    rowN = m_pPointsMat->rows();
    std::vector<bool> matIndex;
    for(int i=0;i<rowN;i++)
        matIndex.push_back(false);

    Eigen::MatrixXd *pointsMat = new (Eigen::MatrixXd)(3,3);

    std::vector<size_t> selectedPoints;
    int Num=0;

    this->m_pPlane = new (std::vector<CPlane*>);
    this->m_pPointsInIndividualPlanes = new (std::vector<Eigen::MatrixXd*>);
    for(int i=0;i<rowN;i++)
    {
       if(matIndex.at(i)!=true)
       {
          selectedPoints.push_back(i);
          (*pointsMat).row(Num)=(*m_pPointsMat).row(i);
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
                  tempMat->row(start+j) = (*m_pPointsMat).row(selectedPoints[j]);
               }
               m_pPointsInIndividualPlanes->push_back(tempMat);
               m_pPlane->push_back(newPlane);
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
     int rowN = m_pPointsMat->rows();
     Point3 temP;
     bool res=true;
     double temp_value;
     std::vector<size_t> AddPInthePlane;
     for(int i=0;i<rowN;i++)
     {
          if(matIndex.at(i)==true)
             continue;

               temP = m_pPointsMat->row(i);
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
        if(m_pPointsMat!=nullptr)
            delete m_pPointsMat;
        if(m_pPlane!=nullptr)
        {
            for(size_t i=0;i<m_pPlane->size();i++)
                delete m_pPlane->at(i);
            m_pPlane->clear();
            delete m_pPlane;
        }
        if(m_pPointsInIndividualPlanes!=nullptr)
        {
            for(size_t i=0;i<m_pPointsInIndividualPlanes->size();i++)
                delete m_pPointsInIndividualPlanes->at(i);
            m_pPointsInIndividualPlanes->clear();
            delete m_pPointsInIndividualPlanes;
        }
}


}
