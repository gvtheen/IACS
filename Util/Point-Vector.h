#ifndef POINT_VECTOR_DEFINE_H
#define POINT_VECTOR_DEFINE_H
#include <Eigen/Core>
#include <Eigen/Dense>
namespace util{

typedef Eigen::Matrix<double, 2, 1> Point2;
typedef Eigen::Matrix<float, 2, 1>  Point2f;
typedef Eigen::Matrix<double, 2, 1> Point2d;

typedef Eigen::Matrix<double, 3, 1> Point3;
typedef Eigen::Matrix<float, 3, 1>  Point3f;
typedef Eigen::Matrix<double, 3, 1> Point3d;
typedef Eigen::Matrix<int, 3, 1>    Point3i;

typedef Eigen::Matrix<double, 2, 1> Vector2;
typedef Eigen::Matrix<float, 2, 1>  Vector2f;
typedef Eigen::Matrix<double, 2, 1> Vector2d;

typedef Eigen::Matrix<double, 3, 1> Vector3;
typedef Eigen::Matrix<float, 3, 1>  Vector3f;
typedef Eigen::Matrix<double, 3, 1> Vector3d;

typedef Eigen::Matrix<double, 4, 1> Vector4;

typedef Eigen::Matrix<double,3,3>   Circle3D;

//typedef Eigen::Matrix<double, 3, 1> Vector3

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;

}
#endif
