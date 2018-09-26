#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/array.hpp>
#include "Constant.h"
#include "../Util/Point-Vector.h"

using util::Point3;
using util::Vector3;
using util::Vector3d;

namespace CATAZJUT{

namespace Geometry{

double distance(const Point3 &a, const Point3 &b);
double distanceSquared(const Point3 &a, const Point3 &b);
double angle(const Vector3 &a, const Vector3 &b);
double angleRadians(const Vector3 &a, const Vector3 &b);
double angle(const Point3 &a, const Point3 &b, const Point3 &c);
double angleRadians(const Point3 &a, const Point3 &b, const Point3 &c);
double torsionAngle(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
double torsionAngleRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
double wilsonAngle(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
double wilsonAngleRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
Point3 midpoint(const Point3 &a, const Point3 &b);
Point3 circumcenter(const Point3 &a, const Point3 &b);
Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c);
Point3 circumcenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
double circumradius(const Point3 &a, const Point3 &b);
double circumradius(const Point3 &a, const Point3 &b, const Point3 &c);
double circumradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
Point3 orthocenter(const Point3 &a, const Point3 &b, double wa, double wb);
Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, double wa, double wb, double wc);
Point3 orthocenter(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, double wa, double wb, double wc, double wd);
double orthoradius(const Point3 &a, const Point3 &b, double wa, double wb);
double orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, double wa, double wb, double wc);
double orthoradius(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d, double wa, double wb, double wc, double wd);
double triangleArea(const Point3 &a, const Point3 &b, const Point3 &c);
double tetrahedronVolume(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
double planeOrientation(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &p);
Vector3 planeNormal(const Point3 &a, const Point3 &b, const Point3 &c);
// --- Derivatives ------------------------------------------------------- //
inline boost::array<Vector3, 2> distanceGradient(const Point3 &a, const Point3 &b);
inline boost::array<Vector3, 3> angleGradient(const Point3 &a, const Point3 &b, const Point3 &c);
inline boost::array<Vector3, 3> angleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c);
inline boost::array<Vector3, 4> torsionAngleGradient(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
inline boost::array<Vector3, 4> torsionAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
inline boost::array<Vector3, 4> wilsonAngleGradient(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
inline boost::array<Vector3, 4> wilsonAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d);
/// Returns the gradient of the distance between points \p a and \p b.
template<typename T> Eigen::Matrix<T, 3, 1> rotate(const Eigen::Matrix<T, 3, 1> &vector, const Eigen::Matrix<T, 3, 1> &axis, T angle);
template<typename T> Eigen::Matrix<T, 3, 1> rotateRadians(const Eigen::Matrix<T, 3, 1> &vector, const Eigen::Matrix<T, 3, 1> &axis, T angle);


inline boost::array<Vector3, 2> distanceGradient(const Point3 &a, const Point3 &b)
{
    boost::array<Vector3, 2> gradient;

    double distance = CATAZJUT::Geometry::distance(a, b);

    gradient[0] = (a - b) / distance;
    gradient[1] = -gradient[0];

    return gradient;
}

/// Returns the gradient of the angle between points \p a, \p b
/// and \p c.
inline boost::array<Vector3, 3> angleGradient(const Point3 &a, const Point3 &b, const Point3 &c)
{
    boost::array<Vector3, 3> gradient = CATAZJUT::Geometry::angleGradientRadians(a, b, c);

    for(size_t i = 0; i < gradient.size(); i++){
        gradient[i] *= CATAZJUT::constants::RadiansToDegrees;
    }

    return gradient;
}

/// Returns the gradient of the angle between points \p a, \p b
/// and \p c.
inline boost::array<Vector3, 3> angleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c)
{
    boost::array<Vector3, 3> gradient;

    double theta = CATAZJUT::Geometry::angleRadians(a, b, c);

    double rab = CATAZJUT::Geometry::distance(a, b);
    double rbc = CATAZJUT::Geometry::distance(b, c);

    gradient[0] = ((((c - b) * rab) - (a - b) * ((b - a).dot(b - c) / rab)) / (pow(rab, 2) * rbc)) / -sin(theta);
    gradient[1] = ((((b - c) + (b - a)) * (rab * rbc) - (((b - a) * (rbc/rab) + (b - c) * (rab/rbc)) * (b - a).dot(b - c))) / pow(rab * rbc, 2)) / -sin(theta);
    gradient[2] = -gradient[0] - gradient[1];

    return gradient;
}

/// Returns the gradient of the torsion angle between the points
/// \p a, \p b, \p c, and \p d.
inline boost::array<Vector3, 4> torsionAngleGradient(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    boost::array<Vector3, 4> gradient = CATAZJUT::Geometry::torsionAngleGradientRadians(a, b, c, d);

    for(size_t i = 0; i < gradient.size(); i++){
        gradient[i] *= CATAZJUT::constants::RadiansToDegrees;
    }

    return gradient;
}

/// Returns the gradient of the torsion angle between the points
/// \p a, \p b, \p c, and \p d.
inline boost::array<Vector3, 4> torsionAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    boost::array<Vector3, 4> gradient;

    double phi = CATAZJUT::Geometry::torsionAngleRadians(a, b, c, d);

    Vector3 ab = b - a;
    Vector3 ac = c - a;
    Vector3 bd = d - b;
    Vector3 cb = b - c;
    Vector3 cd = d - c;

    Vector3 m = ab.cross(cb);
    Vector3 n = cb.cross(cd);

    Vector3 p = ((n / (m.norm() * n.norm())) - ((m / m.squaredNorm()) * cos(phi)));
    Vector3 q = ((m / (m.norm() * n.norm())) - ((n / n.squaredNorm()) * cos(phi)));

    gradient[0] = cb.cross(p) * (1.0 / sin(phi));
    gradient[1] = (ac.cross(p) - cd.cross(q)) * (1.0 / sin(phi));
    gradient[2] = (bd.cross(q) - ab.cross(p)) * (1.0 / sin(phi));
    gradient[3] = cb.cross(q) * (1.0 / sin(phi));

    return gradient;
}

/// Returns the gradient of the wilson angle between the points
/// \p a, \p b, \p c, and \p d.
inline boost::array<Vector3, 4> wilsonAngleGradient(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    boost::array<Vector3, 4> gradient = CATAZJUT::Geometry::wilsonAngleGradientRadians(a, b, c, d);

    for(size_t i = 0; i < gradient.size(); i++){
        gradient[i] *= CATAZJUT::constants::RadiansToDegrees;
    }

    return gradient;
}

/// Returns the gradient of the wilson angle between the points
/// \p a, \p b, \p c, and \p d.
inline boost::array<Vector3, 4> wilsonAngleGradientRadians(const Point3 &a, const Point3 &b, const Point3 &c, const Point3 &d)
{
    Vector3 ba = a - b;
    Vector3 bc = c - b;
    Vector3 bd = d - b;

    double rba = ba.norm();
    double rbc = bc.norm();
    double rbd = bd.norm();

    ba.normalize();
    bc.normalize();
    bd.normalize();

    double theta = acos(ba.dot(bc));

    double w = CATAZJUT::Geometry::wilsonAngleRadians(a, b, c, d);

    boost::array<Vector3, 4> gradient;

    gradient[0] = ((bd.cross(bc) / (cos(w) * sin(theta)) - (ba - bc * cos(theta)) * (tan(w) / pow(sin(theta), 2)))) / rba;
    gradient[2] = ((ba.cross(bd) / (cos(w) * sin(theta)) - (bc - ba * cos(theta)) * (tan(w) / pow(sin(theta), 2)))) / rbc;
    gradient[3] = (bc.cross(ba) / (cos(w) * sin(theta)) - bd * tan(w)) / rbd;
    gradient[1] = -(gradient[0] + gradient[2] + gradient[3]);

    return gradient;
}

// --- Transforms ---------------------------------------------------------- //
template<typename T>
inline Eigen::Matrix<T, 3, 1> rotate(const Eigen::Matrix<T, 3, 1> &vector, const Eigen::Matrix<T, 3, 1> &axis, T angle)
{
    return rotateRadians<T>(vector, axis, angle * CATAZJUT::constants::DegreesToRadians);
}

template<typename T>
inline Eigen::Matrix<T, 3, 1> rotateRadians(const Eigen::Matrix<T, 3, 1> &vector, const Eigen::Matrix<T, 3, 1> &axis, T angle)
{
    return Eigen::Transform<T, 3, 3>(Eigen::AngleAxis<T>(angle, axis)) * vector;
}



}
}

#endif // CGEOMETRY_H
