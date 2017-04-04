#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/readDMAT.h>
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/readTGF.h>
#include <igl/viewer/Viewer.h>
#include <igl/column_to_quats.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>
#include <set>

#include "hrbf_core.h"
#include "hrbf_phi_funcs.h"

#include <igl/copyleft/marching_cubes.h>
#include <igl/parula.h>
#include <igl/random_points_on_mesh.h>
#include <igl/adjacency_list.h>

#include<float.h>
#include "skin_implicit.h"

typedef HRBF_fit<float, 3, Rbf_pow3<float> > HRBF;
struct Bone {
	Eigen::MatrixXd vertices;
	Eigen::MatrixXi faces;
	Eigen::MatrixXd normals;
	HRBF hrbf;
	Eigen::MatrixXd hrbf_points;
	Eigen::MatrixXd hrbf_normals;

	//For debug
	Eigen::MatrixXd hrbf_recon_verts;
	Eigen::MatrixXi hrbf_recon_faces;
};
struct ImplicitMesh {
	Eigen::MatrixXd vertices;
	Eigen::MatrixXi faces;
	Eigen::MatrixXd weights;
	std::vector<Bone> bones;
	Eigen::MatrixXd C;
	Eigen::MatrixXi BE;
	Eigen::VectorXi parents;
	std::vector<std::vector<int>> neighbors;
	std::vector<std::vector<double>> neighbor_weights;

	//Implicit Mesh optimization stuff
	Eigen::MatrixXd betas;
	Eigen::MatrixXd original_offsets;
	Eigen::MatrixXd offsets;
	Eigen::MatrixXd previous_gradients;
};
float reparam(float d) {
	double threshold = 0.5;
	if (d <= -threshold)
		return 1.0f;
	if (d >= threshold)
		return 0.0f;
	d /= threshold;
	return
		(-3.0f / 16.0f) * d * d * d * d * d +
		(5.0f / 8.0f) * d * d * d +
		(-15.0f / 16.0f) * d +
		0.5f;
}
float reparam_grad(float d) {
	double threshold = 3.0;
	if (d >= -threshold && d <= threshold) {
		return 0.0f;
	}	
	d /= threshold;
	return
		(-15.0f / (16.0f*threshold)) * d * d * d * d +
		(15.0f / (8.0f*threshold)) * d * d +
		(-15.0f / (16.0f*threshold));		
}
float distance(const ImplicitMesh &mesh, int bone, const Eigen::MatrixXd &transforms, const Eigen::Vector3d &point, Eigen::Vector3d &gradient) {
	//get transform
	Eigen::MatrixXd T_3 = transforms.block(bone*(3 + 1), 0, 3 + 1, 3).transpose();
	Eigen::MatrixXd T(4, 4);
	T << T_3, 0, 0, 0, 1;
	Eigen::MatrixXd inverse = T.inverse();

	Eigen::Vector4d new_point = inverse * point.homogeneous();
	double distance = mesh.bones[bone].hrbf.eval(HRBF::Vector(new_point[0], new_point[1], new_point[2]));
	HRBF::Vector grad = mesh.bones[bone].hrbf.grad(HRBF::Vector(new_point[0], new_point[1], new_point[2]));

	Eigen::Vector4d grad_2(grad(0), grad(1), grad(2), 0);
	Eigen::Vector4d grad_3 = -(T * grad_2);
	gradient = grad_3.topRows(3);
	return reparam(distance);
}
float vert_distance(const ImplicitMesh &mesh, const Eigen::MatrixXd &vertices, int vertex_idx, const Eigen::MatrixXd &transforms, Eigen::Vector3d &gradient) {
	Eigen::Vector3d g;
	float d, result = -FLT_MAX;
	for (int j = 0; j < mesh.weights.cols();j++) { //Applied the union operator only
	//for (int j = 0; j < 3;j++) {
		if (mesh.weights(vertex_idx, j)>0.01) {
			d = distance(mesh, j, transforms, vertices.row(vertex_idx), g);
			if (d > result) {
				result = d;
				gradient = g;
			}
		}
	}
	return result;
}



typedef std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;

Eigen::MatrixXd bone_colors(4, 3);


const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
Eigen::MatrixXd V, W, C, U, M;
Eigen::MatrixXi T,F, BE;
Eigen::VectorXi P;
std::vector<RotationList > poses;
RotationList hand_pose;
double anim_t = 0.0;
double anim_t_dir = 0.015;

ImplicitMesh mesh;

//Debug Variables
Eigen::MatrixXd N_vertices;
bool show_normals = false;
bool show_partition = false;
bool show_hrbf_points = false;
bool show_hrbfs = false;
bool show_all_hrbfs = false;

//Booleans
bool use_dqs = false;
bool recompute = true;
bool show_weights = false;
int bone_index = 0;
bool debug_implicits = false;
bool use_arm = false;

using namespace Eigen;
using namespace std;

bool pre_draw(igl::viewer::Viewer & viewer)
{	
	if (recompute)
	{
		RotationList vQ;
		vector<Vector3d> vT;
		if (use_arm) {
			// Find pose interval
			const int begin = (int)floor(anim_t) % poses.size();
			const int end = (int)(floor(anim_t) + 1) % poses.size();
			const double t = anim_t - floor(anim_t);

			// Interpolate pose and identity
			RotationList anim_pose(poses[begin].size());
			for (int e = 0;e < poses[begin].size();e++)
			{
				anim_pose[e] = poses[begin][e].slerp(t, poses[end][e]);
			}
			// Propogate relative rotations via FK to retrieve absolute transformations
			igl::forward_kinematics(C, BE, P, anim_pose, vQ, vT);
		}
		else {
			// Interpolate pose and identity
			RotationList anim_pose(hand_pose.size());
			for (int e = 0;e<hand_pose.size();e++)
			{
				anim_pose[e] = hand_pose[e].slerp(anim_t, Quaterniond::Identity());
			}
			igl::forward_kinematics(C, BE, P, anim_pose, vQ, vT);
		}		
		
		const int dim = C.cols();
		MatrixXd T(BE.rows()*(dim + 1), dim);
		for (int e = 0;e<BE.rows();e++)
		{
			Affine3d a = Affine3d::Identity();
			a.translate(vT[e]);
			a.rotate(vQ[e]);
			T.block(e*(dim + 1), 0, dim + 1, dim) =
				a.matrix().transpose().block(0, 0, dim + 1, dim);
		}
		
		// Compute deformation via LBS as matrix multiplication
		if (use_dqs)
		{
			//igl::dqs(V, W, vQ, vT, U);
			//skin_implicit(V, F, W, T, vQ, vT, U);
			//Do implicit skinning here
			//igl::dqs(V, W, vQ, vT, U);
			U = M*T;

			int v_count = mesh.vertices.rows();

			mesh.betas = MatrixXd::Zero(v_count, 1);
			for (int iter = 0; iter < 10; iter++) {

				//Compute new HRBF distance
				for (int vert = 0; vert < U.rows();vert++) {
					Vector3d gradient;
					mesh.offsets(vert) = vert_distance(mesh, U,vert, T, gradient);

					if (iter > 0 && gradient.dot(mesh.previous_gradients.row(vert))< 0.5736f) {
						mesh.betas(vert) = 1.0f;
						continue;
					}
					mesh.previous_gradients.row(vert) = gradient;

					double delta = mesh.original_offsets(vert) - mesh.offsets(vert);
					if (std::abs(delta) < 0.0001f) {
						mesh.betas(vert) = 0.001;
						continue;
					}
					U.row(vert)+= 0.55*delta*gradient;
					mesh.offsets(vert) += 0.55 *delta;
				}
				for (int vert = 0; vert < U.rows();vert++) {
					if (mesh.betas(vert) > 0.0f) {
						continue;
					}
					float mu = std::abs(mesh.original_offsets(vert) - mesh.offsets(vert)) -1.0f;
					mu = std::max(0.0f, 1.0f - mu*mu*mu*mu);

					Eigen::RowVector3d barycenter(0,0,0);
					for (int j = 0;j < mesh.neighbors[vert].size();j++) {
						barycenter += mesh.neighbor_weights[vert][j] * U.row(mesh.neighbors[vert][j]);
					}
					U.row(vert) = ((1.0f - mu) * U.row(vert)) + (mu * barycenter);

				}
			}
			MatrixXd averages(U.rows(), 1);
			for (unsigned n = 0; n < 3; ++n) {
				for (unsigned i = 0; i < U.rows(); ++i) {
					float sum = mesh.betas(i);
					for (unsigned j = 0; j < mesh.neighbors[i].size(); ++j)
						sum += mesh.betas(mesh.neighbors[i][j]);
					averages(i) = sum / (mesh.neighbors[i].size() + 1);
				}
				mesh.betas = averages;				
			}

			// Compute centroids
			MatrixXd centroids(U.rows(), 3);
			for (unsigned i = 0; i < U.rows(); ++i) {
				Eigen::RowVector3d centroid(0,0,0);
				for (unsigned j = 0; j < mesh.neighbors[i].size(); ++j)
					centroid += U.row(mesh.neighbors[i][j]);
				centroids.row(i) =  (1.0f/ (float)mesh.neighbors[i].size())* centroid;
			}

			// Laplacian smoothing
			for (unsigned i = 0; i < U.rows(); ++i)
				U.row(i) = ((1.0f - mesh.betas(i)) * U.row(i)) + (mesh.betas(i) * centroids.row(i));

		}
		else
		{
			U = M*T;
		}

		// Also deform skeleton edges
		MatrixXd CT;
		MatrixXi BET;
		igl::deform_skeleton(C, BE, T, CT, BET);

		viewer.data.clear();
		if (show_partition) {
			viewer.data.set_mesh(mesh.bones[bone_index].vertices, mesh.bones[bone_index].faces);
		}
		else if (show_hrbfs) {
			viewer.data.set_mesh(mesh.bones[bone_index].hrbf_recon_verts, mesh.bones[bone_index].hrbf_recon_faces);
		}
		else if (show_all_hrbfs) {
			int size_V = 0;
			int size_F = 0;
			for (int i = 0; i < mesh.bones.size();i++) {
				size_V += mesh.bones[i].hrbf_recon_verts.rows();
				size_F += mesh.bones[i].hrbf_recon_faces.rows();
			}
			
			MatrixXd all_vertices(size_V,3);
			MatrixXi all_faces(size_F, 3);
			MatrixXd C(all_faces.rows(), 3);
			int vert_row = 0;
			int face_row = 0;
			int color_row = 0;
			int face_accum = 0;
			for (int i = 0;i < mesh.bones.size();i++) {
				int recon_rows = mesh.bones[i].hrbf_recon_verts.rows();
				all_vertices.block(vert_row, 0, recon_rows, 3) = mesh.bones[i].hrbf_recon_verts;
				vert_row += recon_rows;

				int recon_face_rows = mesh.bones[i].hrbf_recon_faces.rows();
				all_faces.block(face_row, 0, recon_face_rows, 3) = mesh.bones[i].hrbf_recon_faces.array() + face_accum;
				face_row += recon_face_rows;
				face_accum += mesh.bones[i].hrbf_recon_verts.rows();

				C.block(color_row, 0, recon_face_rows,3) = bone_colors.row(i%bone_colors.rows()).replicate(recon_face_rows, 1);
				color_row += recon_face_rows;
					
			}			
			
			viewer.data.set_mesh(all_vertices, all_faces);
			viewer.data.set_colors(C);

		}
		else {
			viewer.data.set_mesh(U, F);
		}
		//viewer.data.set_vertices(U);
		viewer.data.set_edges(CT, BET, sea_green);
		viewer.data.compute_normals();
		if (show_normals) {
			igl::per_vertex_normals(U, F, N_vertices);
			viewer.data.add_edges(U, (U + 0.01*N_vertices).eval(), Eigen::RowVector3d(1, 0, 0));
		}		

		if (show_weights) {
			Eigen::MatrixXd bone_weight = mesh.offsets;//W.col(bone_index);
			Eigen::MatrixXd Color;
			igl::parula(bone_weight, bone_weight.minCoeff(), bone_weight.maxCoeff(), Color);
			viewer.data.set_colors(Color);
			viewer.data.add_edges(U, (U + 0.01*mesh.previous_gradients).eval(), Eigen::RowVector3d(1, 0, 0));
		}
		if (show_hrbf_points) {
			viewer.data.set_points(mesh.bones[bone_index].hrbf_points, Eigen::RowVector3d(1, 0, 0));
			viewer.data.add_edges(mesh.bones[bone_index].hrbf_points, (mesh.bones[bone_index].hrbf_points + 0.05*mesh.bones[bone_index].hrbf_normals).eval(), Eigen::RowVector3d(0, 1, 0));
		}
		

		if (viewer.core.is_animating)
		{
			anim_t += anim_t_dir;
			if (!use_arm) {
				anim_t_dir *= (anim_t >= 1.0 || anim_t <= 0.0 ? -1.0 : 1.0);
			}
		}
		else
		{
			recompute = false;
		}
	}
	return false;
}

bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int mods)
{
	recompute = true;
	switch (key)
	{
	case 'D':
	case 'd':
		use_dqs = !use_dqs;
		return true;
	case 'W':
	case 'w':
		show_weights = !show_weights;
		return true;
	case 'E':
	case 'e':
		bone_index = (bone_index + 1 == W.cols() ? 0 : bone_index + 1);
		return true;
	case 'N':
	case 'n':
		show_normals = !show_normals;
		return true;
	case 'P':
	case 'p':
		show_partition = !show_partition;
		return true;
	case 'M':
	case 'm':
		show_hrbf_points = !show_hrbf_points;
		return true;
	case 'B':
	case 'b':
		show_hrbfs = !show_hrbfs;
		return true;
	case 'V':
	case 'v':
		show_all_hrbfs = !show_all_hrbfs;
		return true;		
	case ' ':
		viewer.core.is_animating = !viewer.core.is_animating;
		return true;
	}
	return false;
}




int main(int argc, char *argv[])
{
	if (use_arm) {
		igl::readOBJ("../data/arm.obj", V, F);
		V = 3 * V;
		U = V;
		cout << "Vertices: " << V.rows() << endl;
		igl::readTGF("../data/arm.tgf", C, BE);
		C = 3 * C;
		// retrieve parents for forward kinematics
		igl::directed_edge_parents(BE, P);

		RotationList rest_pose;
		igl::directed_edge_orientations(C, BE, rest_pose);
		poses.resize(4, RotationList(4, Quaterniond::Identity()));
		// poses[1] // twist
		const Quaterniond twist(AngleAxisd(igl::PI, Vector3d(1, 0, 0)));
		poses[1][2] = rest_pose[2] * twist*rest_pose[2].conjugate();
		const Quaterniond bend(AngleAxisd(-igl::PI*0.7, Vector3d(0, 0, 1)));
		poses[3][2] = rest_pose[2] * bend*rest_pose[2].conjugate();

		igl::readDMAT("../data/arm-weights.dmat", W);

		MatrixXd mask = MatrixXd::Zero(W.rows(), W.cols());
		W = (W.array() < 0.01).select(mask, W);
	}
	else {
		igl::readMESH("../data/hand.mesh", V, T, F);
		U = V;
		igl::readTGF("../data/hand.tgf", C, BE);
		// retrieve parents for forward kinematics
		igl::directed_edge_parents(BE, P);

		// Read pose as matrix of quaternions per row
		MatrixXd Q;
		igl::readDMAT("../data/hand-pose.dmat", Q);
		igl::column_to_quats(Q, hand_pose);
		assert(hand_pose.size() == BE.rows());

		igl::readDMAT("../data/hand-weights.dmat", W);
		cout << W.rows() << " " << W.cols() << endl;
		anim_t = 1.0;
		anim_t_dir = -0.03;
	}
		
	igl::lbs_matrix(V, W, M);
	
	igl::per_vertex_normals(U, F, N_vertices);
	cout << "Total faces: " << F.rows() << endl;
	cout << "Total Vertices: " << V.rows() << endl;
	
	mesh.vertices = V;
	mesh.faces = F;
	mesh.weights = W;
	mesh.bones.resize(W.cols());
	mesh.C = C;
	mesh.BE = BE;
	mesh.parents = P;

	///Segment the Mesh
	std::cout << "Segmenting Mesh-> Bones: " << W.cols() << endl;
	//Find the bone that has the most influence on the vertices
	Eigen::MatrixXd max_index(W.rows(), 1);
	Eigen::MatrixXd vertex_bone_index(W.rows(), 1);
	for (int i = 0; i < W.rows();i++) {
		double max_val = W.row(i).maxCoeff(&max_index(i));
	}	
		
	for (int i = 0; i < W.cols();i++) {

		//Get all the vertices that belong to the bone
		int num_vertices = 0;
		for (int j = 0;j < max_index.rows();j++) {
			if (max_index(j) == i) {
				num_vertices++;
			}
		}
		Eigen::MatrixXd bone_vertices(num_vertices, 3);
		int count = 0;
		for (int j = 0;j < max_index.rows();j++) {
			if (max_index(j) == i) {
				bone_vertices.row(count) = V.row(j);
				vertex_bone_index(j) = count;
				count++;
			}			
		}	

		Eigen::SparseMatrix<int> bone_faces;
		std::vector<Eigen::Triplet<int>> triplets;
		count = 0;
		for (int j = 0; j < F.rows();j++) {
			bool add = true;
			for (int k = 0; k < 3;k++) {
				if (max_index(F(j, k)) != i) {
					add = false;
					break;
				}
			}
			if (add) {
				triplets.emplace_back(count, 0, vertex_bone_index(F(j, 0)));
				triplets.emplace_back(count, 1, vertex_bone_index(F(j, 1)));
				triplets.emplace_back(count, 2, vertex_bone_index(F(j, 2)));
				count++;
			}
		}
		bone_faces.resize(count, 3);
		bone_faces.setFromTriplets(triplets.begin(), triplets.end());
	
		mesh.bones[i].vertices = bone_vertices;
		mesh.bones[i].faces = MatrixXi(bone_faces);
		igl::per_vertex_normals(mesh.bones[i].vertices, mesh.bones[i].faces, mesh.bones[i].normals);

		cout << "\tVertices: " << mesh.bones[i].vertices.rows() << 
				" Faces: " << mesh.bones[i].faces.rows() << endl;
	}
	int total_verts = 0, total_faces = 0;
	for (int i = 0;i < mesh.bones.size();i++) {
		total_verts += mesh.bones[i].vertices.rows();
		total_faces += mesh.bones[i].faces.rows();
	}
	cout << "Total Used Vertices: " << total_verts << endl;
	cout << "Total Used Faces: " << total_faces << endl;
	//cout << "C values" << endl;
	//cout << mesh.C << endl;
	//cout << "BE Values" << endl;
	//cout << mesh.BE << endl;

	///Compute HRBFs for each bone
	int num_points = 100;
	for (int i = 0; i < mesh.bones.size(); i++) {	
		Eigen::MatrixXd B;
		Eigen::MatrixXi FI;	

		//sample unique random points on the mesh for the bone				
		igl::random_points_on_mesh(num_points, mesh.bones[i].vertices, mesh.bones[i].faces, B, FI);
		vector<int> vertex_indices;		
		for (int j = 0;j < FI.rows();j++) {
			vertex_indices.push_back(mesh.bones[i].faces(FI(j), 0));			
		}
		std::set<int> s(vertex_indices.begin(), vertex_indices.end());

		Eigen::MatrixXd P(s.size(), 3);
		Eigen::MatrixXd N(s.size(), 3);
		int v_idx = 0;
		for (auto it = s.begin(); it != s.end(); it++) {
			P.row(v_idx) = mesh.bones[i].vertices.row(*it);
			N.row(v_idx) = mesh.bones[i].normals.row(*it);
			v_idx++;
		}
		
		//Remove extremum point that lie too close to the joint
		double h = 0.05;
		Eigen::SparseMatrix<double> new_points;
		Eigen::SparseMatrix<double> new_normals;
		std::vector<Eigen::Triplet<double>> tripletsP;
		std::vector<Eigen::Triplet<double>> tripletsN;
		Eigen::Vector3d bone_diff = mesh.C.row(mesh.BE(i, 1)) - mesh.C.row(mesh.BE(i, 0));
		double square_norm = bone_diff.squaredNorm();
		int count = 0;
		for (int j = 0; j < P.rows(); j++) {
			Eigen::Vector3d vertex_dist = P.row(j) - mesh.C.row(mesh.BE(i, 0));
			double val = vertex_dist.dot(bone_diff) / square_norm;
			if (val > h && val < (1 - h)) {				
				tripletsP.emplace_back(count, 0, P(j, 0));
				tripletsP.emplace_back(count, 1, P(j, 1));
				tripletsP.emplace_back(count, 2, P(j, 2));
				tripletsN.emplace_back(count, 0, N(j, 0));
				tripletsN.emplace_back(count, 1, N(j, 1));
				tripletsN.emplace_back(count, 2, N(j, 2));
				count++;
			}
		}
		//Add HRBF centers near the joints to have proper closure of the distance field
		MatrixXd max_zero = (mesh.bones[i].vertices.rowwise() - mesh.C.row(mesh.BE(i, 0))).rowwise().norm();
		double s_j0 = max_zero.minCoeff();
		MatrixXd max_one = mesh.bones[i].vertices.rowwise() - mesh.C.row(mesh.BE(i, 1));
		double s_j1 = max_one.rowwise().norm().minCoeff();
		RowVector3d n_j0 = mesh.C.row(mesh.BE(i, 0)) - mesh.C.row(mesh.BE(i, 1));
		RowVector3d n_j1 = mesh.C.row(mesh.BE(i, 1)) - mesh.C.row(mesh.BE(i, 0));
		n_j0.normalize(); n_j1.normalize();		
		MatrixXd p_j0 = mesh.C.row(mesh.BE(i, 0)) + (s_j0*n_j0);
		MatrixXd p_j1 = mesh.C.row(mesh.BE(i, 1)) + (s_j1*n_j1);		
		tripletsP.emplace_back(count, 0, p_j0(0));
		tripletsP.emplace_back(count, 1, p_j0(1));
		tripletsP.emplace_back(count, 2, p_j0(2));
		tripletsP.emplace_back(count+1, 0, p_j1(0));
		tripletsP.emplace_back(count+1, 1, p_j1(1));
		tripletsP.emplace_back(count+1, 2, p_j1(2));
		tripletsN.emplace_back(count, 0, n_j0(0));
		tripletsN.emplace_back(count, 1, n_j0(1));
		tripletsN.emplace_back(count, 2, n_j0(2));
		tripletsN.emplace_back(count + 1, 0, n_j1(0));
		tripletsN.emplace_back(count + 1, 1, n_j1(1));
		tripletsN.emplace_back(count + 1, 2, n_j1(2));
		count += 2;

		new_points.resize(count, 3);
		new_points.setFromTriplets(tripletsP.begin(), tripletsP.end());
		new_normals.resize(count, 3);
		new_normals.setFromTriplets(tripletsN.begin(), tripletsN.end());
		P = (MatrixXd)new_points;
		N = (MatrixXd)new_normals;

		mesh.bones[i].hrbf_points = P;	
		mesh.bones[i].hrbf_normals = N;
		cout << "Points: " << P.rows() << endl;

		std::vector<HRBF::Vector> points;
		points.reserve(P.rows());
		for (int j = 0; j < P.rows(); j++) {			
			auto p = P.row(j);
			points.emplace_back(HRBF::Vector(p[0], p[1], p[2]));			
		}		
		
		std::vector<HRBF::Vector> normals;
		normals.reserve(N.rows());
		for (int j = 0; j < N.rows(); j++) {			
			auto n = N.row(j);
			normals.emplace_back(HRBF::Vector(n[0], n[1], n[2]));			
		}

		mesh.bones[i].hrbf.hermite_fit(points, normals);				
		
		if (debug_implicits) {
			// Grid dimensions
			int nx, ny, nz;
			// Maximum extent (side length of bounding box) of points
			double max_extent =
				(P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
			// padding: number of cells beyond bounding box of input points
			const double pad = 8;
			// choose grid spacing (h) so that shortest side gets 30+2*pad samples
			double h_space = max_extent / double(30 + 2 * pad);
			// Place bottom-left-front corner of grid at minimum of points minus padding
			Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad*h_space;
			// Grid dimensions should be at least 3 
			nx = std::max((P.col(0).maxCoeff() - P.col(0).minCoeff() + (2.*pad)*h_space) / h_space, 3.);
			ny = std::max((P.col(1).maxCoeff() - P.col(1).minCoeff() + (2.*pad)*h_space) / h_space, 3.);
			nz = std::max((P.col(2).maxCoeff() - P.col(2).minCoeff() + (2.*pad)*h_space) / h_space, 3.);
			// Compute positions of grid nodes
			Eigen::MatrixXd x(nx*ny*nz, 3);
			for (int w = 0; w < nx; w++)
			{
				for (int j = 0; j < ny; j++)
				{
					for (int k = 0; k < nz; k++)
					{
						// Convert subscript to index
						const auto ind = w + nx*(j + k * ny);
						x.row(ind) = corner + h_space*Eigen::RowVector3d(w, j, k);
					}
				}
			}
			Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);
		
			for (int j = 0; j < g.rows();j++) {
				auto eval = x.row(j);		
				g[j]= mesh.bones[i].hrbf.eval(HRBF::Vector(eval[0], eval[1], eval[2]));		
			}
		
			// Reconstruct mesh		
			Eigen::MatrixXd vert;
			Eigen::MatrixXi face;
			igl::copyleft::marching_cubes(g, x, nx, ny, nz, vert, face);
			mesh.bones[i].hrbf_recon_faces = face;
			mesh.bones[i].hrbf_recon_verts = vert;
		}
		else {
			mesh.bones[i].hrbf_recon_faces = mesh.bones[i].faces;
			mesh.bones[i].hrbf_recon_verts = mesh.bones[i].vertices;
		}
	}
	cout << "Computing Initial Offsets" << endl;
	//Set identity transform to calculate initial signed distances.
	MatrixXd T(BE.rows()*(3 + 1), 3);
	for (int e = 0;e<BE.rows();e++)
	{
		Affine3d a = Affine3d::Identity();		
		T.block(e*(3 + 1), 0, 3 + 1, 3) =
			a.matrix().transpose().block(0, 0, 3 + 1, 3);
	}
	//Pre-calculate initial offsets in the HRBF field.
	MatrixXd offsets(mesh.vertices.rows(), 1);
	for (int i = 0;i < mesh.vertices.rows();i++) {
		Vector3d dummy;
		offsets(i) = vert_distance(mesh, mesh.vertices, i, T, dummy);
	}
	mesh.original_offsets = offsets;
	mesh.previous_gradients = MatrixXd::Zero(mesh.vertices.rows(), 3);
	mesh.offsets = MatrixXd::Zero(mesh.vertices.rows(), 1);

	cout << "Computing Neighbors and Weights" << endl;
	//Compute neighbors 
	igl::adjacency_list(mesh.faces, mesh.neighbors);
	cout <<"Neighbor Vertex Size: " << mesh.neighbors.size() << endl;
	cout <<"Vertices: "<< mesh.vertices.rows() << endl;
	cout << mesh.faces.maxCoeff() << endl;
	// Compute mean value coordinates
	mesh.neighbor_weights.resize(mesh.vertices.rows());
	for (unsigned i = 0; i < mesh.neighbors.size(); i++) {

		// Project neighbors on tangent plane
		std::vector<Eigen::RowVector3d> projections;
		for (unsigned j = 0; j < mesh.neighbors[i].size(); j++) {
			Eigen::RowVector3d delta = mesh.vertices.row(mesh.neighbors[i][j]) - mesh.vertices.row(i);
			projections.push_back(delta - N_vertices.row(i) * N_vertices.row(i).dot(delta));
		}

		// Compute angles
		std::vector<float> angles;
		for (unsigned j = 0; j < projections.size(); ++j) {
			float cosine = projections[j].normalized().dot(projections[(j + 1) % projections.size()]);			
			angles.push_back(std::acos(cosine));
		}

		// Compute barycentric coordinates
		float sum = 0;
		for (unsigned j = 0; j < projections.size(); ++j) {
			float length = projections[j].norm();
			float tan1 = std::tan(angles[(j + projections.size() - 1) % projections.size()] * 0.5f);
			float tan2 = std::tan(angles[j] * 0.5f);
			float weight = (tan1 + tan2) / length;
			mesh.neighbor_weights[i].push_back(weight);
			sum += weight;
		}

		// Normalize weights
		for (unsigned j = 0; j < projections.size(); ++j)
			mesh.neighbor_weights[i][j] /= sum;
	}


	bone_colors << 0.0, 1.0, 1.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0,
		1.0, 1.0, 0.0;

	igl::viewer::Viewer viewer;
	viewer.data.set_mesh(U, F);
	viewer.data.set_edges(C, BE, sea_green);
	
	/*Eigen::MatrixXd Color;
	igl::parula(mesh.original_offsets, mesh.original_offsets.minCoeff(), mesh.original_offsets.maxCoeff(), Color);
	viewer.data.set_colors(Color);*/

	/*Eigen::MatrixXd neigh_points(mesh.neighbors[0].size(), 3);
	for (int i = 0;i < neigh_points.rows();i++) {
		neigh_points.row(i) = mesh.vertices.row(mesh.neighbors[0][i]);
	}
	viewer.data.set_points(neigh_points, Eigen::RowVector3d(1, 0, 0));*/

	viewer.core.show_lines = false;
	viewer.core.show_overlay_depth = false;
	viewer.core.line_width = 1;
	viewer.core.point_size = 10;
	viewer.core.trackball_angle.normalize();
	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	viewer.core.is_animating = false;
	viewer.core.camera_zoom = 2.5;
	viewer.core.animation_max_fps = 30.;
	viewer.core.background_color = Eigen::Vector4f(0.5, 0.5, 0.5,1);
	cout << "Press [d] to toggle between LBS and DQS" << endl <<
		"Press [w] to show weights and [e] to cycle between bones" <<endl <<
		"Press [space] to toggle animation" << endl;
	viewer.launch();
}

////int main(int argc, char *argv[])
////{
////
////	
////	// Load in points + normals from .pwn file
////	Eigen::MatrixXd P, N;
////	/*{
////		Eigen::MatrixXd D;
////		std::vector<std::vector<double> > vD;
////		std::string line;
////		std::fstream in;
////		in.open(argc>1 ? argv[1] : "../data/sphere.pwn");
////		while (in)
////		{
////			std::getline(in, line);
////			std::vector<double> row;
////			std::stringstream stream_line(line);
////			double value;
////			while (stream_line >> value) row.push_back(value);
////			if (!row.empty()) vD.push_back(row);
////		}
////		if (!igl::list_to_matrix(vD, D)) return EXIT_FAILURE;
////		assert(D.cols() == 6 && "pwn file should have 6 columns");
////		P = D.leftCols(3);
////		N = D.rightCols(3);
////	}*/
////	{
////		igl::readOBJ("../data/arm.obj", V, F);
////		P = 3*V;
////
////		igl::per_vertex_normals(V, F, N);
////	}
////	MatrixXd B;
////	MatrixXi FI;
////	igl::random_points_on_mesh(10, P, F, B, FI);	
////	vector<int> vertex_indices;
////	for (int j = 0;j < FI.rows();j++) {
////		vertex_indices.push_back(F(FI(j), 0));		
////	}
////	std::set<int> s(vertex_indices.begin(), vertex_indices.end());
////
////	Eigen::MatrixXd P_N(s.size(), 3);
////	Eigen::MatrixXd N_N(s.size(), 3);
////
////	cout << "Set size: " << s.size() << endl;
////	int index = 0;
////	for (auto it = s.begin(); it != s.end(); it++) {
////		P_N.row(index) = P.row(*it);
////		N_N.row(index) = N.row(*it);
////		index++;
////	}	
////	P = P_N;
////	N = N_N;
////
////	//Construct HRBF Approximation to surface
////	typedef HRBF_fit<float, 3, Rbf_pow3<float> > HRBF;
////
////	std::vector<HRBF::Vector> points;
////	points.reserve(P.rows());
////	for (int i = 0; i < P.rows(); i++) {
////		//if (i % 300 == 0) {
////			auto p = P.row(i);
////			points.emplace_back(HRBF::Vector(p[0], p[1], p[2]));
////		//}
////	}
////	/*std::cout << points.size() << std::endl;
////	std::cout << points[0] << std::endl;*/
////
////	std::vector<HRBF::Vector> normals;
////	normals.reserve(N.rows());
////	for (int i = 0; i < N.rows(); i++) {
////		//if (i % 300 == 0) {
////			auto n = N.row(i);
////			normals.emplace_back(HRBF::Vector(n[0], n[1], n[2]));
////		//}
////	}
////
////	HRBF fit;
////	HRBF::Vector pts[] = { HRBF::Vector(0.f, 0.f, 0.f), HRBF::Vector(1.f, 0.f, 0.f), HRBF::Vector(0.f, 0.f, 2.f) };
////	HRBF::Vector nrs[] = { HRBF::Vector(-1.f, 0.f, 0.f), HRBF::Vector(1.f, 0.f, 0.f), HRBF::Vector(0.f, 0.f, 1.f) };
////	const int size = sizeof(pts) / sizeof(HRBF::Vector);
////	std::cout << pts << std::endl;
////	std::vector<HRBF::Vector> p(pts, pts + size);
////	std::vector<HRBF::Vector> nr(nrs, nrs + size);
////
////	fit.hermite_fit(points, normals);
////
////	const int n = P.rows();
////	// Grid dimensions
////	int nx, ny, nz;
////	// Maximum extent (side length of bounding box) of points
////	double max_extent =
////		(P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
////	// padding: number of cells beyond bounding box of input points
////	const double pad = 8;
////	// choose grid spacing (h) so that shortest side gets 30+2*pad samples
////	double h = max_extent / double(30 + 2 * pad);
////	// Place bottom-left-front corner of grid at minimum of points minus padding
////	Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad*h;
////	// Grid dimensions should be at least 3 
////	nx = std::max((P.col(0).maxCoeff() - P.col(0).minCoeff() + (2.*pad)*h) / h, 3.);
////	ny = std::max((P.col(1).maxCoeff() - P.col(1).minCoeff() + (2.*pad)*h) / h, 3.);
////	nz = std::max((P.col(2).maxCoeff() - P.col(2).minCoeff() + (2.*pad)*h) / h, 3.);
////	// Compute positions of grid nodes
////	Eigen::MatrixXd x(nx*ny*nz, 3);
////	for (int i = 0; i < nx; i++)
////	{
////		for (int j = 0; j < ny; j++)
////		{
////			for (int k = 0; k < nz; k++)
////			{
////				// Convert subscript to index
////				const auto ind = i + nx*(j + k * ny);
////				x.row(ind) = corner + h*Eigen::RowVector3d(i, j, k);
////			}
////		}
////	}
////	Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);
////
////	for (int i = 0; i < g.rows();i++) {
////		auto eval = x.row(i);		
////		g[i]= fit.eval(HRBF::Vector(eval[0], eval[1], eval[2]));		
////	}
////
////	// Reconstruct mesh
////	Eigen::MatrixXd V;
////	Eigen::MatrixXi F;
////	igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
////	//poisson_surface_reconstruction(P, N, V, F);
////
////	// Create a libigl Viewer object to toggle between point cloud and mesh
////	igl::viewer::Viewer viewer;
////	std::cout << R"(
////  P,p      view point cloud
////  M,m      view mesh
////)";
////	const auto set_points = [&]()
////	{
////		viewer.data.clear();
////		viewer.data.set_points(P, Eigen::RowVector3d(1, 1, 1));
////		viewer.data.add_edges(P, (P + 0.01*N).eval(), Eigen::RowVector3d(1, 0, 0));
////	};
////	set_points();
////	viewer.callback_key_pressed = [&](igl::viewer::Viewer&, unsigned int key, int)
////	{
////		switch (key)
////		{
////		case 'P':
////		case 'p':
////			set_points();
////			return true;
////		case 'M':
////		case 'm':
////			viewer.data.clear();
////			viewer.data.set_mesh(V, F);
////			return true;
////		}
////		return false;
////	};
////	viewer.core.point_size = 2;
////	viewer.launch();
////
////	return EXIT_SUCCESS;
////}

//#include <igl/boundary_conditions.h>
//#include <igl/colon.h>
//#include <igl/column_to_quats.h>
//#include <igl/directed_edge_parents.h>
//#include <igl/forward_kinematics.h>
//#include <igl/jet.h>
//#include <igl/lbs_matrix.h>
//#include <igl/deform_skeleton.h>
//#include <igl/normalize_row_sums.h>
//#include <igl/readDMAT.h>
//#include <igl/readMESH.h>
//#include <igl/readTGF.h>
//#include <igl/viewer/Viewer.h>
//#include <igl/bbw.h>
//#include <igl/writeDMAT.h>
//#include <igl/dqs.h>
//#include <igl/directed_edge_orientations.h>
////#include <igl/embree/bone_heat.h>
//
//#include <Eigen/Geometry>
//#include <Eigen/StdVector>
//#include <vector>
//#include <algorithm>
//#include <iostream>
//
//
//bool use_arm = false;
//
//typedef
//std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> >
//RotationList;
//
//const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
//int selected = 0;
//Eigen::MatrixXd V, W, U, C, M;
//Eigen::MatrixXi T, F, BE;
//Eigen::VectorXi P;
//RotationList hand_pose;
//std::vector<RotationList > poses;
//double anim_t = 1.0;
//double anim_t_dir = -0.03;
//
//bool pre_draw(igl::viewer::Viewer & viewer)
//{
//	using namespace Eigen;
//	using namespace std;
//	if (viewer.core.is_animating)
//	{
//
//		RotationList vQ;
//		vector<Vector3d> vT;
//		if(use_arm){
//			//Find pose interval
//			const int begin = (int)floor(anim_t) % poses.size();
//			const int end = (int)(floor(anim_t) + 1) % poses.size();
//			const double t = anim_t - floor(anim_t);
//
//			// Interpolate pose and identity
//			RotationList anim_pose(poses[begin].size());
//			for (int e = 0;e<poses[begin].size();e++)
//			{
//				anim_pose[e] = poses[begin][e].slerp(t, poses[end][e]);
//			}
//			igl::forward_kinematics(C, BE, P, anim_pose, vQ, vT);
//		}
//		else{
//			// Interpolate pose and identity
//			RotationList anim_pose(hand_pose.size());
//			for (int e = 0;e<hand_pose.size();e++)
//			{
//				anim_pose[e] = hand_pose[e].slerp(anim_t, Quaterniond::Identity());
//			}
//			igl::forward_kinematics(C, BE, P, anim_pose, vQ, vT);
//		}
//		// Propogate relative rotations via FK to retrieve absolute transformations
//				
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
//		//U = M*T;
//		igl::dqs(V, W, vQ, vT, U);
//
//		// Also deform skeleton edges
//		MatrixXd CT;
//		MatrixXi BET;
//		igl::deform_skeleton(C, BE, T, CT, BET);
//
//		viewer.data.set_vertices(U);
//		viewer.data.set_edges(CT, BET, sea_green);
//		viewer.data.compute_normals();
//		anim_t += anim_t_dir;
//		if(!use_arm){
//			anim_t_dir *= (anim_t >= 1.0 || anim_t <= 0.0 ? -1.0 : 1.0);
//		}
//	}
//	return false;
//}
//
//void set_color(igl::viewer::Viewer &viewer)
//{
//	Eigen::MatrixXd C;
//	igl::jet(W.col(selected).eval(), true, C);
//	//viewer.data.set_colors(C);
//}
//
//bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int mods)
//{
//	switch (key)
//	{
//	case ' ':
//		viewer.core.is_animating = !viewer.core.is_animating;
//		break;
//	case '.':
//		selected++;
//		selected = std::min(std::max(selected, 0), (int)W.cols() - 1);
//		set_color(viewer);
//		break;
//	case ',':
//		selected--;
//		selected = std::min(std::max(selected, 0), (int)W.cols() - 1);
//		set_color(viewer);
//		break;
//	}
//	return true;
//}
//
//int main(int argc, char *argv[])
//{
//	using namespace Eigen;
//	using namespace std;
//
//	if(use_arm){
//		igl::readOBJ("../data/arm.obj", V, F);
//		V = 3 * V;
//		U = V;
//		cout << "Vertices: " << V.rows() << endl;
//		igl::readTGF("../data/arm.tgf", C, BE);
//		C = 3 * C;
//		// retrieve parents for forward kinematics
//		igl::directed_edge_parents(BE, P);
//		
//		RotationList rest_pose;
//		igl::directed_edge_orientations(C, BE, rest_pose);
//		poses.resize(4, RotationList(4, Quaterniond::Identity()));
//		// poses[1] // twist
//		const Quaterniond twist(AngleAxisd(igl::PI, Vector3d(1, 0, 0)));
//		poses[1][2] = rest_pose[2] * twist*rest_pose[2].conjugate();
//		const Quaterniond bend(AngleAxisd(-igl::PI*0.7, Vector3d(0, 0, 1)));
//		poses[3][2] = rest_pose[2] * bend*rest_pose[2].conjugate();
//
//		igl::readDMAT("../data/arm-weights.dmat", W);
//
//		MatrixXd mask = MatrixXd::Zero(W.rows(), W.cols());
//		W = (W.array() < 0.01).select(mask, W);
//		anim_t =0.0f;
//		anim_t_dir =0.015;
//	}
//	else{
//
//		igl::readMESH("../data/hand.mesh", V, T, F);
//		U = V;
//		igl::readTGF("../data/hand.tgf", C, BE);
//		// retrieve parents for forward kinematics
//		igl::directed_edge_parents(BE, P);
//
//
//
//
//		// Read pose as matrix of quaternions per row
//		MatrixXd Q;
//		igl::readDMAT("../data/hand-pose.dmat", Q);
//		igl::column_to_quats(Q, hand_pose);
//		assert(hand_pose.size() == BE.rows());
//
//		igl::readDMAT("../data/hand-weights.dmat", W);
//		cout << W.rows() << " " << W.cols() << endl;
//	}
//
//	// precompute linear blend skinning matrix
//	igl::lbs_matrix(V, W, M);
//
//	// Plot the mesh with pseudocolors
//	igl::viewer::Viewer viewer;
//	viewer.data.set_mesh(U, F);
//	set_color(viewer);
//	viewer.data.set_edges(C, BE, sea_green);
//	viewer.core.show_lines = false;
//	viewer.core.show_overlay_depth = false;
//	viewer.core.line_width = 1;
//	viewer.callback_pre_draw = &pre_draw;
//	viewer.callback_key_down = &key_down;
//	viewer.core.is_animating = false;
//	viewer.core.animation_max_fps = 30.;
//	cout <<
//		"Press '.' to show next weight function." << endl <<
//		"Press ',' to show previous weight function." << endl <<
//		"Press [space] to toggle animation." << endl;
//	viewer.launch();
//}