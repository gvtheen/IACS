#ifndef GEOMETRY_INLINE_H
#define GEOMETRY_INLINE_H

#include <Eigen/Dense>
#include "Point-Vector.h"
#include "Constant.h"
namespace CATAZJUT{

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


#endif
