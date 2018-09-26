#include <Eigen/Dense>
#include "Point-Vector.h"
#include "Geometry.h"
#include "foreach.h"
#include "Constant.h"
#include "CCartesianCoordinates.h"
#include "CConfigurationBase.h"
#include "CFractionCoordinates.h"
#include "CInternalCoordinates.h"
#include "CPeriodicFramework.h"
#include "CAtom.h"
#include "CUnitCell.h"


namespace CATAZJUT{

CCartesianCoordinates::CCartesianCoordinates()
{
}

/// Creates a new, empty coordinate matrix with space for \p size
/// points.
CCartesianCoordinates::CCartesianCoordinates(CConfigurationBase* m_Conf,size_t size_m)
    : m_pConfiguration(m_Conf)
{
    m_coordinates.insert(m_coordinates.begin(), size_m, Point3().setZero());
}

CCartesianCoordinates::CCartesianCoordinates(CConfigurationBase* m_Conf)
    : m_pConfiguration(m_Conf)
{

}
/// Creates a new coordinate matrix that is a copy of \p coordinates.
CCartesianCoordinates::CCartesianCoordinates(const CCartesianCoordinates &coordinates)
    : m_coordinates(coordinates.m_coordinates)
{
}

/// Destroys the coordinate matrix.
CCartesianCoordinates::~CCartesianCoordinates()
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the size of the matrix to \p size.
void CCartesianCoordinates::resize(size_t size)
{
    m_coordinates.resize(size);
}

/// Returns the number of coordinates in the matrix.
size_t CCartesianCoordinates::size() const
{
    return m_coordinates.size();
}

/// Returns \c true if the matrix is empty.
bool CCartesianCoordinates::isEmpty() const
{
    return m_coordinates.empty();
}

/// Returns a matrix containing the data in the coordinate matrix.
Matrix CCartesianCoordinates::toMatrix() const
{
    Matrix matrix(size(), 3);

    for(size_t i = 0; i < m_coordinates.size(); i++){
        const Point3 &point = m_coordinates[i];

        matrix(i, 0) = point.x();
        matrix(i, 1) = point.y();
        matrix(i, 2) = point.z();
    }

    return matrix;
}

// --- Coordinates --------------------------------------------------------- //
/// Sets the position at \p index to \p position.
void CCartesianCoordinates::setPosition(size_t index, const Point3 &position)
{
    assert(index < size());

    m_coordinates[index] = position;
}

/// Sets the position at \p index to (\p x, \p y, \p z).
void CCartesianCoordinates::setPosition(size_t index, double x, double y, double z)
{
    setPosition(index, Point3(x, y, z));
}

/// Returns the coordinates at \p index.
Point3 CCartesianCoordinates::position(size_t index) const
{
    assert(index < size());

    return m_coordinates[index];
}

/// Sets the value at \p row and \p column to \p value.
void CCartesianCoordinates::setValue(size_t row, size_t column, double value)
{
    assert(row < size());
    assert(column < 3);

    m_coordinates[row][column] = value;
}

/// Returns the value at \p row and \p column;
double CCartesianCoordinates::value(size_t row, size_t column) const
{
    assert(row < size());
    assert(column < 3);

    return m_coordinates[row][column];
}

/// Appends \p position to the coordinates.
void CCartesianCoordinates::append(const Point3 &position)
{
    m_coordinates.push_back(position);
}

/// Appends the point (\p x, \p y, \p z) to the coordinates.
void CCartesianCoordinates::append(double x, double y, double z)
{
    append(Point3(x, y, z));
}

/// Inserts \p position at \p index.
void CCartesianCoordinates::insert(size_t index, const Point3 &position)
{
    m_coordinates.insert(m_coordinates.begin() + index, position);
}

/// Inserts the point (\p x, \p y, \p z) at \p index.
void CCartesianCoordinates::insert(size_t index, double x, double y, double z)
{
    insert(index, Point3(x, y, z));
}

/// Removes the position at \p index.
void CCartesianCoordinates::remove(size_t index)
{
    m_coordinates.erase(m_coordinates.begin() + index);
}
CConfigurationBase* CCartesianCoordinates::Configuration()
{
    return this->m_pConfiguration;
}
void CCartesianCoordinates::setConfiguration(CConfigurationBase*tempConfiguration)
{
    this->m_pConfiguration = tempConfiguration;
}
// --- Geometry ------------------------------------------------------------ //
/// Returns the distance between the points at \p i and \p j. The
/// returned distance is in Angstroms.
double CCartesianCoordinates::distance(size_t i, size_t j) const
{
    return CATAZJUT::Geometry::distance(position(i), position(j));
}
/// Returns the bond angle between the points at \p i, \p j, and
/// \p k. The returned angle is in degrees.
double CCartesianCoordinates::angle(size_t i, size_t j, size_t k) const
{
    return CATAZJUT::Geometry::angle(position(i), position(j), position(k));
}
/// Returns the bond angle between the points at \p i, \p j, and
/// \p k. The returned angle is in radians.
double CCartesianCoordinates::angleRadians(size_t i, size_t j, size_t k) const
{
    return CATAZJUT::Geometry::angleRadians(position(i), position(j), position(k));
}
/// Returns the torsion angle between the points at \p i, \p j, \p k,
/// and \p l. The returned angle is in degrees.
double CCartesianCoordinates::torsionAngle(size_t i, size_t j, size_t k, size_t l) const
{
    return CATAZJUT::Geometry::torsionAngle(position(i), position(j), position(k), position(l));
}

/// Returns the torsion angle between the points at \p i, \p j, \p k,
/// and \p l. The returned angle is in radians.
double CCartesianCoordinates::torsionAngleRadians(size_t i, size_t j, size_t k, size_t l) const
{
    return CATAZJUT::Geometry::torsionAngleRadians(position(i), position(j), position(k), position(l));
}

/// Returns the wilson angle between the points at \p i, \p j, \p k
/// and \p l. The returned angle is in degrees.
double CCartesianCoordinates::wilsonAngle(size_t i, size_t j, size_t k, size_t l) const
{
    return CATAZJUT::Geometry::wilsonAngle(position(i), position(j), position(k), position(l));
}

/// Returns the wilson angle between the points at \p i, \p j, \p k
/// and \p l. The returned angle is in radians.
double CCartesianCoordinates::wilsonAngleRadians(size_t i, size_t j, size_t k, size_t l) const
{
    return CATAZJUT::Geometry::wilsonAngleRadians(position(i), position(j), position(k), position(l));
}

/// Returns the center of the positions in the coordinates. This is
/// also known as the centroid.
Point3 CCartesianCoordinates::center() const
{
    if(isEmpty()){
        return Point3(0, 0, 0);
    }

    Point3 sum(0, 0, 0);

    for(size_t i = 0; i < size(); i++){
        sum += m_coordinates[i];
    }

    return (1.0 / size()) * sum;
}

/// Returns the center of the coordinates after weighting each
/// position with \p weights.
Point3 CCartesianCoordinates::weightedCenter(const std::vector<double> &weights) const
{
    assert(weights.size() == size());

    if(isEmpty()){
        return Point3();
    }

    // sums for each component
    Point3 sum(0, 0, 0);

    // sum of weights
    double sw = 0;

    for(size_t i = 0; i < size(); i++){
        sum += weights[i] * m_coordinates[i];
        sw += weights[i];
    }

    return (1.0 / sw) * sum;
}

/// Moves all of the coordinates by \p vector.
void CCartesianCoordinates::moveBy(const Vector3 &vector)
{
    foreach(Point3 &point, m_coordinates){
        point += vector;
    }
}

/// Moves all of the coordinates by (\p x, \p y, \p z).
void CCartesianCoordinates::moveBy(double x, double y, double z)
{
    moveBy(Vector3(x, y, z));
}

/// Rotates the coordinates by \p angle degrees around \p axis.
void CCartesianCoordinates::rotate(const Vector3 &axis, double angle)
{
    // convert angle to radians
    angle *= CATAZJUT::constants::DegreesToRadians;

    // build rotation transform
    Eigen::Matrix<double, 3, 1> axisVector(axis.x(), axis.y(), axis.z());
    Eigen::Transform<double, 3, 3> transform(Eigen::AngleAxis<double>(angle, axisVector));

    // rotate each point
    for(size_t i = 0; i < m_coordinates.size(); i++){
        setPosition(i, transform * position(i));
    }
}

/// Returns a matrix containing the distances between each pair of
/// points in the coordinates.
Matrix CCartesianCoordinates::distanceMatrix() const
{
    assert(this->size() == m_pConfiguration->m_Atom.size());

    Matrix matrix(size(), size());

    for(size_t i = 0; i < size(); i++){
        // set diagonal entries to zero
        matrix(i, i) = 0;

        for(size_t j = i + 1; j < size(); j++){
            double d = distance(i, j);

            matrix(i, j) = d;
            matrix(j, i) = d;

        }
    }
    return matrix;
}

// --- Derivatives --------------------------------------------------------- //
/// Returns the gradient of the distance between the points at \p i
/// and \p j.
boost::array<Vector3, 2> CCartesianCoordinates::distanceGradient(size_t i, size_t j) const
{
    return CATAZJUT::Geometry::distanceGradient(position(i), position(j));
}

/// Returns the gradient of the angle between the points at \p i,
/// \p j and \p k.
boost::array<Vector3, 3> CCartesianCoordinates::angleGradient(size_t i, size_t j, size_t k) const
{
    return CATAZJUT::Geometry::angleGradient(position(i),
                                            position(j),
                                            position(k));
}

/// Returns the gradient of the angle between the points at \p i,
/// \p j and \p k.
boost::array<Vector3, 3> CCartesianCoordinates::angleGradientRadians(size_t i, size_t j, size_t k) const
{
    return CATAZJUT::Geometry::angleGradientRadians(position(i),
                                                   position(j),
                                                   position(k));
}

/// Returns the gradient of the torsion angle between the points
/// at \p i, \p j, \p k, and \p l.
boost::array<Vector3, 4> CCartesianCoordinates::torsionAngleGradient(size_t i, size_t j, size_t k, size_t l) const
{
    return CATAZJUT::Geometry::torsionAngleGradient(position(i),
                                                   position(j),
                                                   position(k),
                                                   position(l));
}

/// Returns the gradient of the torsion angle between the points at
/// \p i, \p j, \p k, and \p l.
boost::array<Vector3, 4> CCartesianCoordinates::torsionAngleGradientRadians(size_t i, size_t j, size_t k, size_t l) const
{
    return CATAZJUT::Geometry::torsionAngleGradientRadians(position(i),
                                                          position(j),
                                                          position(k),
                                                          position(l));
}

/// Returns the gradient of the wilson angle between the points at
/// \p i, \p j, \p k, and \p l.
boost::array<Vector3, 4> CCartesianCoordinates::wilsonAngleGradient(size_t i, size_t j, size_t k, size_t l) const
{
    return CATAZJUT::Geometry::wilsonAngleGradient(position(i),
                                                  position(j),
                                                  position(k),
                                                  position(l));
}

/// Returns the gradient of the wilson angle between the points at
/// \p i, \p j, \p k, and \p l.
boost::array<Vector3, 4> CCartesianCoordinates::wilsonAngleGradientRadians(size_t i, size_t j, size_t k, size_t l) const
{
    return CATAZJUT::Geometry::wilsonAngleGradientRadians(position(i),
                                                         position(j),
                                                         position(k),
                                                         position(l));
}

// --- Math ---------------------------------------------------------------- //
/// Returns a new coordinate matrix containing the result of adding
/// the coordinates with \p coordinates.
CCartesianCoordinates CCartesianCoordinates::add(const CCartesianCoordinates &coordinates) const
{
    size_t size = std::min(this->size(), coordinates.size());

    CCartesianCoordinates result(m_pConfiguration,size);

    for(size_t i = 0; i < size; i++){
        const Point3 &a = position(i);
        const Point3 &b = coordinates.position(i);

        result.setPosition(i, a + b);
    }

    return result;
}

/// Returns a new coordinate matrix containing the result of
/// subtracting the coordinates with \p coordinates.
CCartesianCoordinates CCartesianCoordinates::subtract(const CCartesianCoordinates &coordinates) const
{
    size_t size = std::min(this->size(), coordinates.size());

    CCartesianCoordinates result(m_pConfiguration,size);

    for(size_t i = 0; i < size; i++){
        const Point3 &a = position(i);
        const Point3 &b = coordinates.position(i);

        result.setPosition(i, a - b);
    }

    return result;
}

/// Returns the 3x3 matrix product of the transpose of the matrix
/// and \p coordinates.
Eigen::Matrix<double, 3, 3> CCartesianCoordinates::multiply(const CCartesianCoordinates *coordinates) const
{
    assert(coordinates->size() == this->size());

    return toMatrix().transpose() * coordinates->toMatrix();
}

// --- Operators ----------------------------------------------------------- //
CCartesianCoordinates CCartesianCoordinates::operator+(const CCartesianCoordinates &coordinates) const
{
    return add(coordinates);
}

CCartesianCoordinates CCartesianCoordinates::operator-(const CCartesianCoordinates &coordinates) const
{
    return subtract(coordinates);
}

CCartesianCoordinates& CCartesianCoordinates::operator=(const CCartesianCoordinates &coordinates)
{
    if(this != &coordinates){
        m_coordinates = coordinates.m_coordinates;
    }

    return *this;
}

/// Returns the position at \p index.
Point3& CCartesianCoordinates::operator[](size_t index)
{
    return m_coordinates[index];
}

/// \overload
const Point3& CCartesianCoordinates::operator[](size_t index) const
{
    return m_coordinates[index];
}

CInternalCoordinates* CCartesianCoordinates::toInternalCoordinates()
{
    size_t size_Cartesian = this->size();
    CInternalCoordinates* res = new CInternalCoordinates(this->m_pConfiguration,size_Cartesian);
    double dis,angle,torsion;
    if (size_Cartesian>0)
    {
       res->setCoordinates(0,0,0,0);
       res->setConnections(0,0,0,0);
       if(size_Cartesian>1)
       {
           res->setConnections(1,1,0,0);
           dis = CATAZJUT::Geometry::distance(m_coordinates[1],m_coordinates[0]);
           res->setCoordinates(1,dis,0,0);

           if(size_Cartesian>2)
           {
               res->setConnections(2,2,1,0);
                 dis = CATAZJUT::Geometry::distance(m_coordinates[2],m_coordinates[1]);
               angle = CATAZJUT::Geometry::angle(m_coordinates[2],m_coordinates[1],m_coordinates[0]);
               res->setCoordinates(2,dis,angle,0);
           }
       }
    }
    for(size_t i=3;i<size_Cartesian;i++)
    {
        std::vector<size_t> vect (std::move(CheckConnection(i)));
        if(vect.size()<3)
            for(size_t j=0;j<3;j++)
                vect.push_back(i-j);

           dis = CATAZJUT::Geometry::distance(m_coordinates[i],m_coordinates[vect[0]]);
         angle = CATAZJUT::Geometry::angle(m_coordinates[i],m_coordinates[vect[0]],m_coordinates[vect[1]]);
       torsion = CATAZJUT::Geometry::torsionAngle(m_coordinates[i],m_coordinates[vect[0]],\
                                        m_coordinates[vect[1]],m_coordinates[vect[2]]);
       res->setCoordinates(i,dis,angle,torsion);
       res->setConnections(i,vect[0],vect[1],vect[2]);
    }
    return res;
}
std::vector<size_t>  CCartesianCoordinates::CheckConnection(size_t currentIndex)
{
   std::vector<size_t> res;
   bool isBack=false;

   std::vector<CAtom*> &tmpAtom = m_pConfiguration->m_Atom;

   foreach(const CAtom* neighbors_1, tmpAtom[currentIndex]->neighbors())
       if(neighbors_1->index() < currentIndex){
          foreach(const CAtom* neighbors_2,neighbors_1->neighbors())
          {
             if(neighbors_2->index() < currentIndex)
             {
                foreach(const CAtom* neighbors_3,neighbors_2->neighbors())
                {
                    if(neighbors_3->index() < currentIndex)
                    {
                           res.push_back(neighbors_1->index());
                           res.push_back(neighbors_2->index());
                           res.push_back(neighbors_3->index());
                           isBack=true;
                    }
                    if(isBack) break;
                }
            }
            if(isBack) break;
         }
         if(isBack) break;
       }
   return res;
}
CFractionCoordinates* CCartesianCoordinates::toFractionCoordinates(CPeriodicFramework* mPeriodicFramework)
{
    if(mPeriodicFramework->dimensionalType()==CATAZJUT::DEFINED::Molecule)
        return nullptr;

    Eigen::Matrix<double, 3, 3> transitionMat(std::move(mPeriodicFramework->m_pUnitCell->MatrixOfBravaisLattice()));
   // double scalingFactor = mPeriodicFramework->m_pUnitCell->scalingFactor();
         Matrix *cartMat = new Matrix(std::move(this->toMatrix()));

                *cartMat = (transitionMat.inverse())*(cartMat->transpose());
     //(3,N)---> (N,3)
     cartMat->transpose();
     CFractionCoordinates  *mFraction = new CFractionCoordinates(mPeriodicFramework);
     for(size_t i=0;i<cartMat->rows();i++)
        mFraction->append(cartMat->row(i));
     // adjust the coordinate into the cell;
     mFraction->adjustCoordinateIntoCell();

     delete cartMat;
     return mFraction;
}
void CCartesianCoordinates::clear()
{
    m_coordinates->clear();
}


}   //namespace of CATAZJUT
