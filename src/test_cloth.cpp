#include "test_cloth.h"
#include "hrbf_core.h"
#include "hrbf_phi_funcs.h"

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
#include <iostream>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/doublearea.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>


namespace cloth_sim {
	/* Some physics constants */
	#define DAMPING 0.01 // how much to damp the cloth simulation each frame
	#define TIME_STEPSIZE2 0.5*0.5 // how large time step each particle takes each frame
	#define CONSTRAINT_ITERATIONS 5 // how many iterations of constraint satisfaction each frame (more is rigid, less is soft)

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
	


	class Vec3 // a minimal vector class of 3 floats and overloaded math operators
	{
	public:
		float f[3];

		Vec3(float x, float y, float z)
		{
			f[0] = x;
			f[1] = y;
			f[2] = z;
		}

		Vec3() {}

		float length()
		{
			return sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
		}

		Vec3 normalized()
		{
			float l = length();
			return Vec3(f[0] / l, f[1] / l, f[2] / l);
		}

		void operator+= (const Vec3 &v)
		{
			f[0] += v.f[0];
			f[1] += v.f[1];
			f[2] += v.f[2];
		}

		Vec3 operator/ (const float &a)
		{
			return Vec3(f[0] / a, f[1] / a, f[2] / a);
		}

		Vec3 operator- (const Vec3 &v)
		{
			return Vec3(f[0] - v.f[0], f[1] - v.f[1], f[2] - v.f[2]);
		}

		Vec3 operator+ (const Vec3 &v)
		{
			return Vec3(f[0] + v.f[0], f[1] + v.f[1], f[2] + v.f[2]);
		}

		Vec3 operator* (const float &a)
		{
			return Vec3(f[0] * a, f[1] * a, f[2] * a);
		}

		Vec3 operator-()
		{
			return Vec3(-f[0], -f[1], -f[2]);
		}

		Vec3 cross(const Vec3 &v)
		{
			return Vec3(f[1] * v.f[2] - f[2] * v.f[1], f[2] * v.f[0] - f[0] * v.f[2], f[0] * v.f[1] - f[1] * v.f[0]);
		}

		float dot(const Vec3 &v)
		{
			return f[0] * v.f[0] + f[1] * v.f[1] + f[2] * v.f[2];
		}
	};

	/* The particle class represents a particle of mass that can move around in 3D space*/
	class Particle
	{
	private:
		bool movable; // can the particle move or not ? used to pin parts of the cloth

		float mass; // the mass of the particle (is always 1 in this example)
		Vec3 pos; // the current position of the particle in 3D space
		Vec3 old_pos; // the position of the particle in the previous time step, used as part of the verlet numerical integration scheme
		Vec3 acceleration; // a vector representing the current acceleration of the particle
		Vec3 accumulated_normal; // an accumulated normal (i.e. non normalized), used for OpenGL soft shading

	public:
		Particle(Vec3 pos) : pos(pos), old_pos(pos), acceleration(Vec3(0, 0, 0)), mass(1), movable(true), accumulated_normal(Vec3(0, 0, 0)) {}
		Particle() {}

		void addForce(Vec3 f)
		{
			acceleration += f / mass;
		}

		/* This is one of the important methods, where the time is progressed a single step size (TIME_STEPSIZE)
		The method is called by Cloth.time_step()
		Given the equation "force = mass * acceleration" the next position is found through verlet integration*/
		void timeStep()
		{
			if (movable)
			{
				Vec3 temp = pos;
				pos = pos + (pos - old_pos)*(1.0 - DAMPING) + acceleration*TIME_STEPSIZE2;
				old_pos = temp;
				acceleration = Vec3(0, 0, 0); // acceleration is reset since it HAS been translated into a change in position (and implicitely into velocity)	
			}
		}

		Vec3& getPos() { return pos; }

		void resetAcceleration() { acceleration = Vec3(0, 0, 0); }

		void offsetPos(const Vec3 v) { if (movable) pos += v; }

		void makeUnmovable() { movable = false; }

		void addToNormal(Vec3 normal)
		{
			accumulated_normal += normal.normalized();
		}

		Vec3& getNormal() { return accumulated_normal; } // notice, the normal is not unit length

		void resetNormal() { accumulated_normal = Vec3(0, 0, 0); }

	};

	class Constraint
	{
	private:
		float rest_distance; // the length between particle p1 and p2 in rest configuration

	public:
		Particle *p1, *p2; // the two particles that are connected through this constraint

		Constraint(Particle *p1, Particle *p2) : p1(p1), p2(p2)
		{
			Vec3 vec = p1->getPos() - p2->getPos();
			rest_distance = vec.length();
		}

		/* This is one of the important methods, where a single constraint between two particles p1 and p2 is solved
		the method is called by Cloth.time_step() many times per frame*/
		void satisfyConstraint()
		{
			Vec3 p1_to_p2 = p2->getPos() - p1->getPos(); // vector from p1 to p2
			float current_distance = p1_to_p2.length(); // current distance between p1 and p2
			Vec3 correctionVector = p1_to_p2*(1 - rest_distance / current_distance); // The offset vector that could moves p1 into a distance of rest_distance to p2
			Vec3 correctionVectorHalf = correctionVector*0.5; // Lets make it half that length, so that we can move BOTH p1 and p2.
			p1->offsetPos(correctionVectorHalf); // correctionVectorHalf is pointing from p1 to p2, so the length should move p1 half the length needed to satisfy the constraint.
			p2->offsetPos(-correctionVectorHalf); // we must move p2 the negative direction of correctionVectorHalf since it points from p2 to p1, and not p1 to p2.	
		}

	};

	class Cloth
	{
	private:

		int num_particles_width; // number of particles in "width" direction
		int num_particles_height; // number of particles in "height" direction
								  // total number of particles is num_particles_width*num_particles_height

		std::vector<Particle> particles; // all particles that are part of this cloth
		std::vector<Constraint> constraints; // alle constraints between particles as part of this cloth

		Particle* getParticle(int x, int y) { return &particles[y*num_particles_width + x]; }
		void makeConstraint(Particle *p1, Particle *p2) { constraints.push_back(Constraint(p1, p2)); }


		/* A private method used by drawShaded() and addWindForcesForTriangle() to retrieve the
		normal vector of the triangle defined by the position of the particles p1, p2, and p3.
		The magnitude of the normal vector is equal to the area of the parallelogram defined by p1, p2 and p3
		*/
		Vec3 calcTriangleNormal(Particle *p1, Particle *p2, Particle *p3)
		{
			Vec3 pos1 = p1->getPos();
			Vec3 pos2 = p2->getPos();
			Vec3 pos3 = p3->getPos();

			Vec3 v1 = pos2 - pos1;
			Vec3 v2 = pos3 - pos1;

			return v1.cross(v2);
		}

		/* A private method used by windForce() to calcualte the wind force for a single triangle
		defined by p1,p2,p3*/
		void addWindForcesForTriangle(Particle *p1, Particle *p2, Particle *p3, const Vec3 direction)
		{
			Vec3 normal = calcTriangleNormal(p1, p2, p3);
			Vec3 d = normal.normalized();
			Vec3 force = normal*(d.dot(direction));
			p1->addForce(force);
			p2->addForce(force);
			p3->addForce(force);
		}

		/* A private method used by drawShaded(), that draws a single triangle p1,p2,p3 with a color*/
		void drawTriangle(Particle *p1, Particle *p2, Particle *p3, const Vec3 color)
		{
			glColor3fv((GLfloat*)&color);

			glNormal3fv((GLfloat *)&(p1->getNormal().normalized()));
			glVertex3fv((GLfloat *)&(p1->getPos()));

			glNormal3fv((GLfloat *)&(p2->getNormal().normalized()));
			glVertex3fv((GLfloat *)&(p2->getPos()));

			glNormal3fv((GLfloat *)&(p3->getNormal().normalized()));
			glVertex3fv((GLfloat *)&(p3->getPos()));
		}

	public:

		/* This is a important constructor for the entire system of particles and constraints*/
		Cloth(float posx, float posy, float width, float height, int num_particles_width, int num_particles_height) : num_particles_width(num_particles_width), num_particles_height(num_particles_height)
		{
			particles.resize(num_particles_width*num_particles_height); //I am essentially using this vector as an array with room for num_particles_width*num_particles_height particles

																		// creating particles in a grid of particles from (0,0,0) to (width,-height,0)
			for (int x = 0; x<num_particles_width; x++)
			{
				for (int y = 0; y<num_particles_height; y++)
				{
					Vec3 pos = Vec3(posx + (width * (x/(float)num_particles_width)),
					-height * (y/(float)num_particles_height) + posy,
					0);
					/*Vec3 pos = Vec3(width * (x / (float)num_particles_width),
						0, -height * (y / (float)num_particles_height));*/
					particles[y*num_particles_width + x] = Particle(pos); // insert particle in column x at y'th row
				}
			}

			// Connecting immediate neighbor particles with constraints (distance 1 and sqrt(2) in the grid)
			for (int x = 0; x<num_particles_width; x++)
			{
				for (int y = 0; y<num_particles_height; y++)
				{
					if (x<num_particles_width - 1) makeConstraint(getParticle(x, y), getParticle(x + 1, y));
					if (y<num_particles_height - 1) makeConstraint(getParticle(x, y), getParticle(x, y + 1));
					if (x<num_particles_width - 1 && y<num_particles_height - 1) makeConstraint(getParticle(x, y), getParticle(x + 1, y + 1));
					if (x<num_particles_width - 1 && y<num_particles_height - 1) makeConstraint(getParticle(x + 1, y), getParticle(x, y + 1));
				}
			}


			// Connecting secondary neighbors with constraints (distance 2 and sqrt(4) in the grid)
			for (int x = 0; x<num_particles_width; x++)
			{
				for (int y = 0; y<num_particles_height; y++)
				{
					if (x<num_particles_width - 2) makeConstraint(getParticle(x, y), getParticle(x + 2, y));
					if (y<num_particles_height - 2) makeConstraint(getParticle(x, y), getParticle(x, y + 2));
					if (x<num_particles_width - 2 && y<num_particles_height - 2) makeConstraint(getParticle(x, y), getParticle(x + 2, y + 2));
					if (x<num_particles_width - 2 && y<num_particles_height - 2) makeConstraint(getParticle(x + 2, y), getParticle(x, y + 2));
				}
			}
			
			for (int i = 0;i<num_particles_height; i++)
			{
				//getParticle(0 ,i)->offsetPos(Vec3(0.5,0.0,0.0)); // moving the particle a bit towards the center, to make it hang more natural - because I like it ;)
				getParticle(0 ,i)->makeUnmovable(); 			
			}
		}
	
		
		void get_draw(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
			Eigen::SparseMatrix<float> vertices;
			Eigen::SparseMatrix<int> faces;
			std::vector<Eigen::Triplet<float>> v_triplets;
			std::vector<Eigen::Triplet<int>> f_triplets;

			v_triplets.resize(num_particles_width*num_particles_height * 3);
			int count = 0;
			for (int x = 0; x<num_particles_width; x++)
			{
				for (int y = 0; y<num_particles_height; y++)
				{
					Particle * p1 = getParticle(x, y);					
					v_triplets.emplace_back(count, 0, p1->getPos().f[0]);
					v_triplets.emplace_back(count, 1, p1->getPos().f[1]);
					v_triplets.emplace_back(count, 2, p1->getPos().f[2]);
				
					count++;
				}
			}
			vertices.resize(num_particles_width*num_particles_height, 3);
			vertices.setFromTriplets(v_triplets.begin(), v_triplets.end());
			V = Eigen::MatrixXd(vertices);

			count = 0;
			for (int x = 0; x<num_particles_width-1; x++)
			{
				for (int y = 0; y<num_particles_height-1; y++)
				{
					int t1 = x * num_particles_height + y;
					int t2 = x * num_particles_height + (y + 1);
					int t3 = (x + 1) *num_particles_height + y;
					int t4 = (x + 1)*num_particles_height + (y + 1);

					f_triplets.emplace_back(count, 0, t1);
					f_triplets.emplace_back(count, 1, t2);
					f_triplets.emplace_back(count, 2, t3);

					f_triplets.emplace_back(count+1, 0, t2);
					f_triplets.emplace_back(count+1, 1, t4);
					f_triplets.emplace_back(count+1, 2, t3);

					count+=2;
				}
			}
			faces.resize(count, 3);
			faces.setFromTriplets(f_triplets.begin(), f_triplets.end());
			F = Eigen::MatrixXi(faces);
		}

		/* this is an important methods where the time is progressed one time step for the entire cloth.
		This includes calling satisfyConstraint() for every constraint, and calling timeStep() for all particles
		*/
		void timeStep()
		{
			std::vector<Constraint>::iterator constraint;
			for (int i = 0; i<CONSTRAINT_ITERATIONS; i++) // iterate over all constraints several times
			{
				for (constraint = constraints.begin(); constraint != constraints.end(); constraint++)
				{
					(*constraint).satisfyConstraint(); // satisfy constraint.
				}
			}

			std::vector<Particle>::iterator particle;
			for (particle = particles.begin(); particle != particles.end(); particle++)
			{
				(*particle).timeStep(); // calculate the position of each particle at the next time step.
			}
		}

		/* used to add gravity (or any other arbitrary vector) to all particles*/
		void addForce(const Vec3 direction)
		{
			std::vector<Particle>::iterator particle;
			for (particle = particles.begin(); particle != particles.end(); particle++)
			{
				(*particle).addForce(direction); // add the forces to each particle
			}

		}

		/* used to add wind forces to all particles, is added for each triangle since the final force is proportional to the triangle area as seen from the wind direction*/
		void windForce(const Vec3 direction)
		{
			for (int x = 0; x<num_particles_width - 1; x++)
			{
				for (int y = 0; y<num_particles_height - 1; y++)
				{
					addWindForcesForTriangle(getParticle(x + 1, y), getParticle(x, y), getParticle(x, y + 1), direction);
					addWindForcesForTriangle(getParticle(x + 1, y + 1), getParticle(x + 1, y), getParticle(x, y + 1), direction);
				}
			}
		}

		/* used to detect and resolve the collision of the cloth with the ball.
		This is based on a very simples scheme where the position of each particle is simply compared to the sphere and corrected.
		This also means that the sphere can "slip through" if the ball is small enough compared to the distance in the grid bewteen particles
		*/
		void ballCollision(const Vec3 center, const float radius)
		{
			std::vector<Particle>::iterator particle;
			for (particle = particles.begin(); particle != particles.end(); particle++)
			{
				Vec3 v = (*particle).getPos() - center;
				float l = v.length();
				if (v.length() < radius) // if the particle is inside the ball
				{
					(*particle).offsetPos(v.normalized()*(radius - l)); // project the particle to the surface of the ball
				}
			}
		}		
	};




	Cloth cloth1(-10,5,7, 5, 25, 20); // one Cloth object of the Cloth class

	typedef std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;
	Eigen::MatrixXd bone_colors(4, 3);


	Eigen::MatrixXd V, U, W, C, M;
	Eigen::MatrixXi T, F, F1, BE;
	Eigen::VectorXi P;
	std::vector<RotationList > poses;
	RotationList hand_pose;
	double anim_t = 0.0;
	double anim_t_dir = 0.015;
	ImplicitMesh mesh;

	bool step = false;
	bool smoothing = false;
	

	bool pre_draw(igl::viewer::Viewer & viewer)
	{
		using namespace std;	
		using namespace Eigen;
		//Cloth simulation updates
		cloth1.get_draw(V, F);		

		Eigen::MatrixXd N;
		igl::per_vertex_normals(V, F, N);		
		/*viewer.data.set_mesh(V, F);
		viewer.data.set_normals(N);*/
		if (step) {
			cloth1.addForce(Vec3(0, -0.2, 0)*TIME_STEPSIZE2); // add gravity each frame, pointing down	
			cloth1.windForce(Vec3(0.2, 0.1, 0.01)*TIME_STEPSIZE2); // generate some wind each frame
			cloth1.timeStep(); // calculate the particle positions of the next frame
		}
		

		//Animation updates
		MatrixXd all_vertices(V.rows() + U.rows(), 3);
		MatrixXi all_faces(F.rows() + F1.rows(),3);
		all_vertices.block(0, 0, V.rows(), V.cols()) = V;
		all_vertices.block(V.rows(), 0, U.rows(), U.cols()) = U;
		all_faces.block(0, 0, F.rows(), F.cols()) = F;
		all_faces.block(F.rows(), 0, F1.rows(), F1.cols()) = F1.array() + V.rows();

		MatrixXd C(all_faces.rows(), 3);
		C.block(0, 0, F.rows(), 3) = bone_colors.row(0).replicate(F.rows(), 1);
		C.block(F.rows(), 0, F1.rows(), 3) = bone_colors.row(1).replicate(F1.rows(), 1);

		Eigen::MatrixXd N2;
		igl::per_vertex_normals(U, F1, N2);
		MatrixXd normals(all_vertices.rows(),3);
		normals.block(0, 0, N.rows(), N.cols()) = N;
		normals.block(N.rows(), 0, N2.rows(), N2.cols()) = N2;

		viewer.data.set_mesh(all_vertices, all_faces);
		viewer.data.set_normals(normals);
		viewer.data.set_colors(C);

		return false;
	}
	bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int mods)
	{		
		using namespace std;
		switch (key)
		{
			case ' ':
				viewer.core.is_animating = !viewer.core.is_animating;
				return true;
			case 'D':
			case 'd':
				step = !step;
				return true;
			case 'S':
			case 's':
				smoothing = !smoothing;
				cout << "Smothing val is: " << smoothing << endl;
				return true;
			case 'A':
			case 'a':
				cloth1.addForce(Vec3(0, -0.2, 0)*TIME_STEPSIZE2); // add gravity each frame, pointing down	
				cloth1.windForce(Vec3(0.2, 0.1, 0.01)*TIME_STEPSIZE2); // generate some wind each frame
				cloth1.timeStep(); // calculate the particle positions of the next frame	
				return true;
		}
		return false;
	}
	int test_cloth() {

		using namespace std;
		using namespace Eigen;

		cout << "Running Cloth Simulation" << endl;
		igl::readOBJ("../data/arm.obj", U, F1);
		U = 3 * U;				
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
		igl::lbs_matrix(V, W, M);

		mesh.vertices = V;
		mesh.faces = F1;
		mesh.weights = W;
		mesh.bones.resize(W.cols());
		mesh.C = C;
		mesh.BE = BE;
		mesh.parents = P;
	

		bone_colors << 1.0, 1.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0,
			1.0, 1.0, 0.0;

		igl::viewer::Viewer viewer;
		viewer.core.show_lines = false;
		viewer.core.show_overlay_depth = false;
		viewer.core.line_width = 1;
		viewer.core.point_size = 10;
		viewer.core.trackball_angle.normalize();
		viewer.core.is_animating = false;
		viewer.callback_pre_draw = &pre_draw;
		viewer.callback_key_down = &key_down;
		viewer.core.camera_zoom = 2.5;
		viewer.core.animation_max_fps = 60.;
		viewer.core.background_color = Eigen::Vector4f(1.0, 1.0, 1.0, 1);
		cout << "Press [d] to turn on automatic update" << endl <<
			"Press [a] to increment individual time step" << endl <<
			"Press [space] to toggle animation" << endl;
		viewer.launch();

		return 0;
	}
}