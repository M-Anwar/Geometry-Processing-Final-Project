#include "skin_implicit.h"

#include <igl/dqs.h>
#include <iostream>

void skin_implicit(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & W,
    const Eigen::MatrixXd & T,
    const std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> & vQ,
    const std::vector<Eigen::Vector3d> & vT,
    Eigen::MatrixXd & U
){
	using namespace std;
	cout << "Running Implicit Skinning Iteration" << endl;

	Eigen::MatrixXd vertices;
    igl::dqs(V,W,vQ,vT,vertices);
	//vertices = M*T;
	U = vertices;
}