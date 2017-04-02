#ifndef SKIN_IMPLICIT_HPP__
#define SKIN_IMPLICIT_HPP__

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

void skin_implicit(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & W,
    const Eigen::MatrixXd & T,
    const std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> & vQ,
    const std::vector<Eigen::Vector3d> & vT,
    Eigen::MatrixXd & U

);


#endif // SKIN_IMPLICIT_HPP__