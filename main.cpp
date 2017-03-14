//#include "biharmonic_precompute.h"
//#include "biharmonic_solve.h"
//#include "arap_precompute.h"
//#include "arap_single_iteration.h"
//#include <igl/min_quad_with_fixed.h>
//#include <igl/read_triangle_mesh.h>
//#include <igl/viewer/Viewer.h>
//#include <igl/project.h>
//#include <igl/unproject.h>
//#include <igl/snap_points.h>
//#include <igl/unproject_onto_mesh.h>
//#include <Eigen/Core>
//#include <iostream>
//#include <stack>
//
//// Undoable
//struct State
//{
//  // Rest and transformed control points
//  Eigen::MatrixXd CV, CU;
//  bool placing_handles = true;
//} s;
//
//int main(int argc, char *argv[])
//{
//  // Undo Management
//  std::stack<State> undo_stack,redo_stack;
//  const auto push_undo = [&](State & _s=s)
//  {
//    undo_stack.push(_s);
//    // clear
//    redo_stack = std::stack<State>();
//  };
//  const auto undo = [&]()
//  {
//    if(!undo_stack.empty())
//    {
//      redo_stack.push(s);
//      s = undo_stack.top();
//      undo_stack.pop();
//    }
//  };
//  const auto redo = [&]()
//  {
//    if(!redo_stack.empty())
//    {
//      undo_stack.push(s);
//      s = redo_stack.top();
//      redo_stack.pop();
//    }
//  };
//
//  Eigen::MatrixXd V,U;
//  Eigen::MatrixXi F;
//  long sel = -1;
//  Eigen::RowVector3f last_mouse;
//  igl::min_quad_with_fixed_data<double> biharmonic_data, arap_data;
//  Eigen::SparseMatrix<double> arap_K;
//
//  // Load input meshes
//  igl::read_triangle_mesh(
//    (argc>1?argv[1]:"../shared/data/decimated-knight.off"),V,F);
//  U = V;
//  igl::viewer::Viewer viewer;
//  std::cout<<R"(
//[click]  To place new control point
//[drag]   To move control point
//[space]  Toggle whether placing control points or deforming
//M,m      Switch deformation methods
//U,u      Update deformation (i.e., run another iteration of solver)
//R,r      Reset control points 
//⌘ Z      Undo
//⌘ ⇧ Z    Redo
//)";
//  enum Method
//  {
//    BIHARMONIC = 0,
//    ARAP = 1,
//    NUM_METHODS = 2,
//  } method = BIHARMONIC;
//
//  const auto & update = [&]()
//  {
//    // predefined colors
//    const Eigen::RowVector3d orange(1.0,0.7,0.2);
//    const Eigen::RowVector3d yellow(1.0,0.9,0.2);
//    const Eigen::RowVector3d blue(0.2,0.3,0.8);
//    const Eigen::RowVector3d green(0.2,0.6,0.3);
//    if(s.placing_handles)
//    {
//      viewer.data.set_vertices(V);
//      viewer.data.set_colors(blue);
//      viewer.data.set_points(s.CV,orange);
//    }else
//    {
//      // SOLVE FOR DEFORMATION
//      switch(method)
//      {
//        default:
//        case BIHARMONIC:
//        {
//          Eigen::MatrixXd D;
//          biharmonic_solve(biharmonic_data,s.CU-s.CV,D);
//          U = V+D;
//          break;
//        }
//        case ARAP:
//        {
//          arap_single_iteration(arap_data,arap_K,s.CU,U);
//          break;
//        }
//      }
//      viewer.data.set_vertices(U);
//      viewer.data.set_colors(method==BIHARMONIC?orange:yellow);
//      viewer.data.set_points(s.CU,method==BIHARMONIC?blue:green);
//    }
//    viewer.data.compute_normals();
//  };
//  viewer.callback_mouse_down = 
//    [&](igl::viewer::Viewer&, int, int)->bool
//  {
//    last_mouse = Eigen::RowVector3f(
//      viewer.current_mouse_x,viewer.core.viewport(3)-viewer.current_mouse_y,0);
//    if(s.placing_handles)
//    {
//      // Find closest point on mesh to mouse position
//      int fid;
//      Eigen::Vector3f bary;
//      if(igl::unproject_onto_mesh(
//        last_mouse.head(2),
//        viewer.core.view * viewer.core.model,
//        viewer.core.proj, 
//        viewer.core.viewport, 
//        V, F, 
//        fid, bary))
//      {
//        long c;
//        bary.maxCoeff(&c);
//        Eigen::RowVector3d new_c = V.row(F(fid,c));
//        if(s.CV.size()==0 || (s.CV.rowwise()-new_c).rowwise().norm().minCoeff() > 0)
//        {
//          push_undo();
//          s.CV.conservativeResize(s.CV.rows()+1,3);
//          // Snap to closest vertex on hit face
//          s.CV.row(s.CV.rows()-1) = new_c;
//          update();
//          return true;
//        }
//      }
//    }else
//    {
//      // Move closest control point
//      Eigen::MatrixXf CP;
//      igl::project(
//        Eigen::MatrixXf(s.CU.cast<float>()),
//        viewer.core.view * viewer.core.model, 
//        viewer.core.proj, viewer.core.viewport, CP);
//      Eigen::VectorXf D = (CP.rowwise()-last_mouse).rowwise().norm();
//      sel = (D.minCoeff(&sel) < 30)?sel:-1;
//      if(sel != -1)
//      {
//        last_mouse(2) = CP(sel,2);
//        push_undo();
//        update();
//        return true;
//      }
//    }
//    return false;
//  };
//
//  viewer.callback_mouse_move = [&](igl::viewer::Viewer &, int,int)->bool
//  {
//    if(sel!=-1)
//    {
//      Eigen::RowVector3f drag_mouse(
//        viewer.current_mouse_x,
//        viewer.core.viewport(3) - viewer.current_mouse_y,
//        last_mouse(2));
//      Eigen::RowVector3f drag_scene,last_scene;
//      igl::unproject(
//        drag_mouse,
//        viewer.core.view*viewer.core.model,
//        viewer.core.proj,
//        viewer.core.viewport,
//        drag_scene);
//      igl::unproject(
//        last_mouse,
//        viewer.core.view*viewer.core.model,
//        viewer.core.proj,
//        viewer.core.viewport,
//        last_scene);
//      s.CU.row(sel) += (drag_scene-last_scene).cast<double>();
//      last_mouse = drag_mouse;
//      update();
//      return true;
//    }
//    return false;
//  };
//  viewer.callback_mouse_up = [&](igl::viewer::Viewer&, int, int)->bool
//  {
//    sel = -1;
//    return false;
//  };
//  viewer.callback_key_pressed = 
//    [&](igl::viewer::Viewer &, unsigned int key, int mod)
//  {
//    switch(key)
//    {
//      case 'M':
//      case 'm':
//      {
//        method = (Method)(((int)(method)+1)%((int)(NUM_METHODS)));
//        break;
//      }
//      case 'R':
//      case 'r':
//      {
//        push_undo();
//        s.CU = s.CV;
//        break;
//      }
//      case 'U':
//      case 'u':
//      {
//        // Just trigger an update
//        break;
//      }
//      case ' ':
//        push_undo();
//        s.placing_handles ^= 1;
//        if(!s.placing_handles && s.CV.rows()>0)
//        {
//          // Switching to deformation mode
//          s.CU = s.CV;
//          Eigen::VectorXi b;
//          igl::snap_points(s.CV,V,b);
//          // PRECOMPUTATION FOR DEFORMATION
//          biharmonic_precompute(V,F,b,biharmonic_data);
//          arap_precompute(V,F,b,arap_data,arap_K);
//        }
//        break;
//      default:
//        return false;
//    }
//    update();
//    return true;
//  };
//
//  // Special callback for handling undo
//  viewer.callback_key_down = 
//    [&](igl::viewer::Viewer &, unsigned char key, int mod)->bool
//  {
//    if(key == 'Z' && (mod & GLFW_MOD_SUPER))
//    {
//      (mod & GLFW_MOD_SHIFT) ? redo() : undo();
//      update();
//      return true;
//    }
//    return false;
//  };
//  viewer.callback_pre_draw = 
//    [&](igl::viewer::Viewer &)->bool
//  {
//    if(viewer.core.is_animating && !s.placing_handles && method == ARAP)
//    {
//      arap_single_iteration(arap_data,arap_K,s.CU,U);
//      update();
//    }
//    return false;
//  };
//  viewer.data.set_mesh(V,F);
//  viewer.core.show_lines = false;
//  viewer.core.is_animating = true;
//  viewer.data.face_based = true;
//  update();
//  viewer.launch();
//  return EXIT_SUCCESS;
//}

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

typedef
std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> >
RotationList;

const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
Eigen::MatrixXd V, W, C, U, M;
Eigen::MatrixXi F, BE;
Eigen::VectorXi P;
std::vector<RotationList > poses;
double anim_t = 0.0;
double anim_t_dir = 0.015;
bool use_dqs = false;
bool recompute = true;

bool pre_draw(igl::viewer::Viewer & viewer)
{
	using namespace Eigen;
	using namespace std;
	if (recompute)
	{
		// Find pose interval
		const int begin = (int)floor(anim_t) % poses.size();
		const int end = (int)(floor(anim_t) + 1) % poses.size();
		const double t = anim_t - floor(anim_t);

		// Interpolate pose and identity
		RotationList anim_pose(poses[begin].size());
		for (int e = 0;e<poses[begin].size();e++)
		{
			anim_pose[e] = poses[begin][e].slerp(t, poses[end][e]);
		}
		// Propogate relative rotations via FK to retrieve absolute transformations
		RotationList vQ;
		vector<Vector3d> vT;
		igl::forward_kinematics(C, BE, P, anim_pose, vQ, vT);
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
			igl::dqs(V, W, vQ, vT, U);
		}
		else
		{
			U = M*T;
		}

		// Also deform skeleton edges
		MatrixXd CT;
		MatrixXi BET;
		igl::deform_skeleton(C, BE, T, CT, BET);

		viewer.data.set_vertices(U);
		viewer.data.set_edges(CT, BET, sea_green);
		viewer.data.compute_normals();
		if (viewer.core.is_animating)
		{
			anim_t += anim_t_dir;
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
	case ' ':
		viewer.core.is_animating = !viewer.core.is_animating;
		return true;
	}
	return false;
}

int main(int argc, char *argv[])
{
	using namespace Eigen;
	using namespace std;
	igl::readOBJ("../shared/data/arm.obj", V, F);
	U = V;
	igl::readTGF("../shared/data/arm.tgf", C, BE);
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

	igl::readDMAT("../shared/data/arm-weights.dmat", W);
	igl::lbs_matrix(V, W, M);

	// Plot the mesh with pseudocolors
	igl::viewer::Viewer viewer;
	viewer.data.set_mesh(U, F);
	viewer.data.set_edges(C, BE, sea_green);
	viewer.core.show_lines = false;
	viewer.core.show_overlay_depth = false;
	viewer.core.line_width = 1;
	viewer.core.trackball_angle.normalize();
	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	viewer.core.is_animating = false;
	viewer.core.camera_zoom = 2.5;
	viewer.core.animation_max_fps = 30.;
	cout << "Press [d] to toggle between LBS and DQS" << endl <<
		"Press [space] to toggle animation" << endl;
	viewer.launch();
}

