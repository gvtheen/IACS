#ifndef CFRACTIONCOORDINATES_H
#define CFRACTIONCOORDINATES_H
#include<vector>
#include "Point-Vector.h"

namespace CATAZJUT{

class CInternalCoordinates;
class CCartesianCoordinates;
class CUnitCell;
class CPeriodicFramework;

class CFractionCoordinates
{
    public:
        CFractionCoordinates();
        CFractionCoordinates(CPeriodicFramework*);
        CFractionCoordinates(CPeriodicFramework* mconf,size_t size_m);
        CFractionCoordinates(const CFractionCoordinates& mth);

        virtual ~CFractionCoordinates();

        void resize(size_t size);
        size_t size() const;
        bool isEmpty() const;

        void setPeriodicFramework(CPeriodicFramework*);
        CPeriodicFramework* PeriodicFramework() const;


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
        void clear();
        //
        double movePointIntoCell(double temp);
          void adjustCoordinateIntoCell();

        Point3 toCartesian(const size_t& index)const;
        Point3 toCartesian(const Point3& position)const;
        Point3 toCartesian(const double x, const double y, const double z)const;
        //from cartesion to Fraction
        Point3 toFractionFromCart(const Point3& position);
        Point3 toFractionFromCart(double x, double y, double z);

        double distance(size_t i, size_t j);
        double angle(size_t i, size_t j, size_t k);
        double angleRadians(size_t i, size_t j, size_t k);
        double torsionAngle(size_t i, size_t j, size_t k, size_t l);
        double torsionAngleRadians(size_t i, size_t j, size_t k, size_t l);
        double wilsonAngle(size_t i, size_t j, size_t k, size_t l);
        double wilsonAngleRadians(size_t i, size_t j, size_t k, size_t l);
        Point3 center();
        Point3 weightedCenter(const std::vector<double> &weights);
        void moveBy(const Vector3 &vector);
        void moveBy(double x, double y, double z);
        void rotate(const Vector3 &axis, double angle);

        CFractionCoordinates add(const CFractionCoordinates&) const;
        CFractionCoordinates substract(const CFractionCoordinates&)const;
        Eigen::Matrix<double,3,3> multiply(const CFractionCoordinates&)const;

        CInternalCoordinates*  toInternalCoordinates()  const;
        CCartesianCoordinates* toCartesianCoordinates() const;

        CFractionCoordinates& operator = (const CFractionCoordinates&);
        CFractionCoordinates  operator + (const CFractionCoordinates&) const;
        CFractionCoordinates  operator - (const CFractionCoordinates&) const;

    private:
        CPeriodicFramework  *m_pPeriodicFramework;
        std::vector<Point3>  m_coordinates;
};



}//namespace of CATAZJUT
#endif // CFRACTIONCOORDINATES_H
