#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/readTGF.h>
#include <igl/viewer/Viewer.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

#include "hrbf_core.h"
#include "hrbf_phi_funcs.h"

#include <igl/copyleft/marching_cubes.h>

//typedef
//std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> >
//RotationList;
//
//const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
//Eigen::MatrixXd V, W, C, U, M;
//Eigen::MatrixXi F, BE;
//Eigen::VectorXi P;
//std::vector<RotationList > poses;
//double anim_t = 0.0;
//double anim_t_dir = 0.015;
//bool use_dqs = false;
//bool recompute = true;
//
//bool pre_draw(igl::viewer::Viewer & viewer)
//{
//	using namespace Eigen;
//	using namespace std;
//	if (recompute)
//	{
//		// Find pose interval
//		const int begin = (int)floor(anim_t) % poses.size();
//		const int end = (int)(floor(anim_t) + 1) % poses.size();
//		const double t = anim_t - floor(anim_t);
//
//		// Interpolate pose and identity
//		RotationList anim_pose(poses[begin].size());
//		for (int e = 0;e<poses[begin].size();e++)
//		{
//			anim_pose[e] = poses[begin][e].slerp(t, poses[end][e]);
//		}
//		// Propogate relative rotations via FK to retrieve absolute transformations
//		RotationList vQ;
//		vector<Vector3d> vT;
//		igl::forward_kinematics(C, BE, P, anim_pose, vQ, vT);
//		const int dim = C.cols();
//		MatrixXd T(BE.rows()*(dim + 1), dim);
//		for (int e = 0;e<BE.rows();e++)
//		{
//			Affine3d a = Affine3d::Identity();
//			a.translate(vT[e]);
//			a.rotate(vQ[e]);
//			T.block(e*(dim + 1), 0, dim + 1, dim) =
//				a.matrix().transpose().block(0, 0, dim + 1, dim);
//		}
//		// Compute deformation via LBS as matrix multiplication
//		if (use_dqs)
//		{
//			igl::dqs(V, W, vQ, vT, U);
//		}
//		else
//		{
//			U = M*T;
//		}
//
//		// Also deform skeleton edges
//		MatrixXd CT;
//		MatrixXi BET;
//		igl::deform_skeleton(C, BE, T, CT, BET);
//
//		viewer.data.set_vertices(U);
//		viewer.data.set_edges(CT, BET, sea_green);
//		viewer.data.compute_normals();
//		if (viewer.core.is_animating)
//		{
//			anim_t += anim_t_dir;
//		}
//		else
//		{
//			recompute = false;
//		}
//	}
//	return false;
//}
//
//bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int mods)
//{
//	recompute = true;
//	switch (key)
//	{
//	case 'D':
//	case 'd':
//		use_dqs = !use_dqs;
//		return true;
//	case ' ':
//		viewer.core.is_animating = !viewer.core.is_animating;
//		return true;
//	}
//	return false;
//}
//
//int main(int argc, char *argv[])
//{
//	using namespace Eigen;
//	using namespace std;
//	igl::readOBJ("../data/arm.obj", V, F);
//	U = V;
//	igl::readTGF("../data/arm.tgf", C, BE);
//	// retrieve parents for forward kinematics
//	igl::directed_edge_parents(BE, P);
//	RotationList rest_pose;
//	igl::directed_edge_orientations(C, BE, rest_pose);
//	poses.resize(4, RotationList(4, Quaterniond::Identity()));
//	// poses[1] // twist
//	const Quaterniond twist(AngleAxisd(igl::PI, Vector3d(1, 0, 0)));
//	poses[1][2] = rest_pose[2] * twist*rest_pose[2].conjugate();
//	const Quaterniond bend(AngleAxisd(-igl::PI*0.7, Vector3d(0, 0, 1)));
//	poses[3][2] = rest_pose[2] * bend*rest_pose[2].conjugate();
//
//	igl::readDMAT("../data/arm-weights.dmat", W);
//	igl::lbs_matrix(V, W, M);
//
//	// Plot the mesh with pseudocolors
//	igl::viewer::Viewer viewer;
//	viewer.data.set_mesh(U, F);
//	viewer.data.set_edges(C, BE, sea_green);
//	viewer.core.show_lines = false;
//	viewer.core.show_overlay_depth = false;
//	viewer.core.line_width = 1;
//	viewer.core.trackball_angle.normalize();
//	viewer.callback_pre_draw = &pre_draw;
//	viewer.callback_key_down = &key_down;
//	viewer.core.is_animating = false;
//	viewer.core.camera_zoom = 2.5;
//	viewer.core.animation_max_fps = 30.;
//	cout << "Press [d] to toggle between LBS and DQS" << endl <<
//		"Press [space] to toggle animation" << endl;
//	viewer.launch();
//}

int main(int argc, char *argv[])
{
	// Load in points + normals from .pwn file
	Eigen::MatrixXd P, N;
	{
		Eigen::MatrixXd D;
		std::vector<std::vector<double> > vD;
		std::string line;
		std::fstream in;
		in.open(argc>1 ? argv[1] : "../data/hand.pwn");
		while (in)
		{
			std::getline(in, line);
			std::vector<double> row;
			std::stringstream stream_line(line);
			double value;
			while (stream_line >> value) row.push_back(value);
			if (!row.empty()) vD.push_back(row);
		}
		if (!igl::list_to_matrix(vD, D)) return EXIT_FAILURE;
		assert(D.cols() == 6 && "pwn file should have 6 columns");
		P = D.leftCols(3);
		N = D.rightCols(3);
	}

	//Construct HRBF Approximation to surface
	typedef HRBF_fit<float, 3, Rbf_pow3<float> > HRBF;

	std::vector<HRBF::Vector> points;
	points.reserve(P.rows());
	for (int i = 0; i < P.rows(); i++) {
		if (i % 1000 == 0) {
			auto p = P.row(i);
			points.emplace_back(HRBF::Vector(p[0], p[1], p[2]));
		}
	}
	/*std::cout << points.size() << std::endl;
	std::cout << points[0] << std::endl;*/

	std::vector<HRBF::Vector> normals;
	normals.reserve(N.rows());
	for (int i = 0; i < N.rows(); i++) {
		if (i % 1000 == 0) {
			auto n = N.row(i);
			normals.emplace_back(HRBF::Vector(n[0], n[1], n[2]));
		}
	}

	HRBF fit;
	HRBF::Vector pts[] = { HRBF::Vector(0.f, 0.f, 0.f), HRBF::Vector(1.f, 0.f, 0.f), HRBF::Vector(0.f, 0.f, 2.f) };
	HRBF::Vector nrs[] = { HRBF::Vector(-1.f, 0.f, 0.f), HRBF::Vector(1.f, 0.f, 0.f), HRBF::Vector(0.f, 0.f, 1.f) };
	const int size = sizeof(pts) / sizeof(HRBF::Vector);
	std::cout << pts << std::endl;
	std::vector<HRBF::Vector> p(pts, pts + size);
	std::vector<HRBF::Vector> nr(nrs, nrs + size);

	fit.hermite_fit(points, normals);

	const int n = P.rows();
	// Grid dimensions
	int nx, ny, nz;
	// Maximum extent (side length of bounding box) of points
	double max_extent =
		(P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
	// padding: number of cells beyond bounding box of input points
	const double pad = 8;
	// choose grid spacing (h) so that shortest side gets 30+2*pad samples
	double h = max_extent / double(30 + 2 * pad);
	// Place bottom-left-front corner of grid at minimum of points minus padding
	Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad*h;
	// Grid dimensions should be at least 3 
	nx = std::max((P.col(0).maxCoeff() - P.col(0).minCoeff() + (2.*pad)*h) / h, 3.);
	ny = std::max((P.col(1).maxCoeff() - P.col(1).minCoeff() + (2.*pad)*h) / h, 3.);
	nz = std::max((P.col(2).maxCoeff() - P.col(2).minCoeff() + (2.*pad)*h) / h, 3.);
	// Compute positions of grid nodes
	Eigen::MatrixXd x(nx*ny*nz, 3);
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int k = 0; k < nz; k++)
			{
				// Convert subscript to index
				const auto ind = i + nx*(j + k * ny);
				x.row(ind) = corner + h*Eigen::RowVector3d(i, j, k);
			}
		}
	}
	Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

	for (int i = 0; i < g.rows();i++) {
		auto eval = x.row(i);		
		g[i]= fit.eval(HRBF::Vector(eval[0], eval[1], eval[2]));		
	}

	// Reconstruct mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
	//poisson_surface_reconstruction(P, N, V, F);

	// Create a libigl Viewer object to toggle between point cloud and mesh
	igl::viewer::Viewer viewer;
	std::cout << R"(
  P,p      view point cloud
  M,m      view mesh
)";
	const auto set_points = [&]()
	{
		viewer.data.clear();
		viewer.data.set_points(P, Eigen::RowVector3d(1, 1, 1));
		viewer.data.add_edges(P, (P + 0.01*N).eval(), Eigen::RowVector3d(1, 0, 0));
	};
	set_points();
	viewer.callback_key_pressed = [&](igl::viewer::Viewer&, unsigned int key, int)
	{
		switch (key)
		{
		case 'P':
		case 'p':
			set_points();
			return true;
		case 'M':
		case 'm':
			viewer.data.clear();
			viewer.data.set_mesh(V, F);
			return true;
		}
		return false;
	};
	viewer.core.point_size = 2;
	viewer.launch();

	return EXIT_SUCCESS;
}