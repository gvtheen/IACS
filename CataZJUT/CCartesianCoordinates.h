#ifndef CCARTESIANCOORDINATES_H
#define CCARTESIANCOORDINATES_H
#include <Eigen/Dense>
#include <boost/array.hpp>
#include <vector>
#include "Point-Vector.h"

#include "CAtom.h"
namespace CATAZJUT{

class CConfigurationBase;
class CInternalCoordinates;
class CUnitCell;
class CPeriodicFramework;

class CCartesianCoordinates
{
public:
// construction and destruction
    CCartesianCoordinates();
    CCartesianCoordinates(CConfigurationBase* mconf,size_t size_m);
    CCartesianCoordinates(CConfigurationBase* mconf);
    CCartesianCoordinates(const CCartesianCoordinates &coordinates);
    ~CCartesianCoordinates();

    // properties
    void resize(size_t size);
    size_t size() const;
    void clear();
    bool isEmpty() const;
    Matrix toMatrix() const;

    // coordinates
    void setPosition(size_t index, const Point3 &position);
    void setPosition(size_t index, double x, double y, double z);
    Point3 position(size_t index) const;
    void setValue(size_t row, size_t column, double value);
    double value(size_t row, size_t column) const;
    void append(const Point3 &position);
    void append(double x, double y, double z);
    void insert(size_t index, const Point3 &position);
    void insert(size_t index, double x, double y, double z);
    void remove(size_t index);
    CConfigurationBase* Configuration();
    void setConfiguration(CConfigurationBase*temp_pConfiguration);
    // geometry
    double distance(size_t i, size_t j) const;
    double angle(size_t i, size_t j, size_t k) const;
    double angleRadians(size_t i, size_t j, size_t k) const;
    double torsionAngle(size_t i, size_t j, size_t k, size_t l) const;
    double torsionAngleRadians(size_t i, size_t j, size_t k, size_t l) const;
    double wilsonAngle(size_t i, size_t j, size_t k, size_t l) const;
    double wilsonAngleRadians(size_t i, size_t j, size_t k, size_t l) const;
    Point3 center() const;
    Point3 weightedCenter(const std::vector<double> &weights) const;
    void moveBy(const Vector3 &vector);
    void moveBy(double x, double y, double z);
    void rotate(const Vector3 &axis, double angle);
    Matrix distanceMatrix() const;

    // derivatives
    boost::array<Vector3, 2> distanceGradient(size_t i, size_t j) const;
    boost::array<Vector3, 3> angleGradient(size_t i, size_t j, size_t k) const;
    boost::array<Vector3, 3> angleGradientRadians(size_t i, size_t j, size_t k) const;
    boost::array<Vector3, 4> torsionAngleGradient(size_t i, size_t j, size_t k, size_t l) const;
    boost::array<Vector3, 4> torsionAngleGradientRadians(size_t i, size_t j, size_t k, size_t l) const;
    boost::array<Vector3, 4> wilsonAngleGradient(size_t i, size_t j, size_t k, size_t l) const;
    boost::array<Vector3, 4> wilsonAngleGradientRadians(size_t i, size_t j, size_t k, size_t l) const;

    // math
    CCartesianCoordinates add(const CCartesianCoordinates &coordinates) const;
    CCartesianCoordinates subtract(const CCartesianCoordinates &coordinates) const;
    Eigen::Matrix<double, 3, 3> multiply(const CCartesianCoordinates *coordinates) const;

    // operators
    CCartesianCoordinates  operator+(const CCartesianCoordinates &coordinates) const;
    CCartesianCoordinates  operator-(const CCartesianCoordinates &coordinates) const;
    CCartesianCoordinates& operator=(const CCartesianCoordinates &coordinates);
                   Point3& operator[](size_t index);
    const Point3& operator[](size_t index) const;

    CInternalCoordinates* toInternalCoordinates();
    std::vector<size_t>   CheckConnection(size_t);
    CFractionCoordinates* toFractionCoordinates(CPeriodicFramework*);


private:
    std::vector<Point3>        m_coordinates;
    CConfigurationBase        *m_pConfiguration;
};


}
#endif // CCARTESIANCOORDINATES_H
