// Frickin 5th Project File Lets Make this Work.cpp : This file contains the 'main' function. Program execution begins and ends there.
// ALRIGHT WE IN THE LAST STRETCH BROSKI LES MAKE PRIYANKA FRIICKIN PROUD SONNNNNNNNNNNN


//probably dont need all these libraries idk lmao
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <Eigen/Dense>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <chrono>



using namespace Eigen;

class Disk {


public:
	float x, y;
	//Disks will x and y values to find a and b values for the unit vectors?
	float theta, Mnet, Fnet, arclength;
	Matrix<float, 3, 3> G_mat;
	Vector3d coords;
	//Create a function to find Fi-1. According to equation 9, i to i-1 unit vectors will point towards i-1 (origin). Use x,y to find unit vectors
	//From that, can find Moments of the driving forces
	Disk(float th, float arclength, float Fnet, float Mnet, Matrix<float, 3, 3> G_mat) :theta(th), arclength(arclength), Fnet(Fnet), Mnet(Mnet) {}

	Vector3d coord_get(float arclength, float th, Matrix<float, 3, 3> G_mat) {

		float x = (arclength / th) * sin(th);
		float y = (arclength / th) * (1 - cos(th));

		Vector3d coords(x, y, 0);

	}



};

class Robot {

public:
	//Figure out all variables going in. Radius? Need it for inertia. Ask priyanka

	float FI, FII, d_radius, arclength, theta, g, E; // E is youngs modulus
	int disk_num;
	Matrix<Disk, 1, Dynamic> disks;
	Vector2d weights;
	Matrix<float, 2, 1> f_tend;
	Robot(float l, float r, float g, int d, float E, Vector2d weights, Matrix<float, 2, 1> f_tend) :d_radius(r), arclength(l), disk_num(d), g(g), E(E), weights(weights), f_tend(f_tend) {

		
	}

	Matrix<float, 3, 3> get_transmat(int q, int p, float l, float theta) {

		if (q < p) {
			std::cout << "bad input";
		}
		else {

			Matrix<float, 3, 3> T = Matrix<float, 3, 3>::Identity(); //Set a "global T" probably a better way to do this, i did it the python way i guess

			for (int i = p; i < q; i++) {

				Matrix<float, 3, 3> m0;  //m0 is Ri, m1 is Ri-1 
							//placed 0.6 as a test till i figure out how to make unknowns carry over for the solver
				m0 << cos(theta), -sin(theta), (l / theta)* sin(theta),
					sin(theta), cos(theta), (l / theta)* (1 - cos(theta)),
					0, 0, 1;

				T = T * m0;

			}
			return T;
		}
	}

	Matrix<Vector3f, 2, 1> get_segmentfm(Matrix<Vector3d, 2, 1> p_ten, Matrix<float, 2, 1> f_tend, int i, int n, Matrix<float, 2, 1> weights, Matrix<Vector3d, 2, 1> prev_MF, float l, float th) {
		//gets current force and moment on each tendon 

		Matrix<float, 3, 3> T_mat = get_transmat(i, i - 1, l, th); //figure out how to get rid of the need for an l
		Vector3f prev_coord(T_mat(2, 0), T_mat(2, 1), 0);

		Vector3f ai = T_mat * (p_ten[0]);
		Vector3f bi = T_mat * (p_ten[1]); //p_ten houses tendon locations

			//figure out how to combine ai and bi 

		ai = (-1 * ai) + p_ten(0);  //ai to ai-1
		bi = (-1 * bi) + p_ten(1);

		Vector3f F1_rel; //combine into one matrix to allow for one operation 
		Vector3f F2_rel;

		if (i == n) { //i = segment #, n = total disks


			F1_rel = (f_tend[0] * (ai)) / ai.norm();
			F2_rel = (f_tend[1] * (bi)) / (bi).norm();

		}

		else if (i < n) {
			Matrix<float, 3, 3> T_mat1 = get_transmat(i + 1, i, l, th); // i+1 to i-1
			T_mat1 = T_mat * T_mat1;

			Vector3f ai_1 = T_mat1 * (p_ten[0]); // instead of i+1 to i-1, reuse i to i-1 for claculations
			Vector3f bi_1 = T_mat1 * (p_ten[1]); //p_ten houses tendon locations

			ai_1 = ai_1 - ai;  //i to i+1 vector
			bi_1 = bi_1 + bi;


			F1_rel = f_tend[0] * ((ai / ai.norm()) + (ai_1 / ai_1.norm()));
			F2_rel = f_tend[1] * ((bi / bi.norm()) + (bi_1 / bi_1.norm()));

		}

		Matrix<float, 3, 3> G_mat = get_transmat(i - 1, 0, l, th);

		Matrix<float, 3, 1> Fg_d = G_mat * weights[0]; //force g on disk i in reference of i-1 MAKE SURE TO DEFINE MASSES[0] AS [0, mass*g, 0]
		Matrix<float, 3, 1> Mg_d = prev_coord.cross(Fg_d); //moment of Fg_d

		Matrix<float, 3, 1>	MFi = prev_coord.cross(prev_MF[1]);

		Vector3f M1_rel = ai.cross(F1_rel);
		Vector3f M2_rel = bi.cross(F2_rel);

		prev_MF[1] = T_mat * prev_MF[1]; //puts previous force in terms of i-1 reference frame
		prev_MF[0] = T_mat * prev_MF[0]; //previous moment in terms of i-1 reference frame

		Vector3f F_net = F1_rel + F2_rel + prev_MF[1] + Fg_d;
		Vector3f M_net = M1_rel + M2_rel + prev_MF[0] + Mg_d + MFi;

		Matrix<Vector3f, 2, 1> nets(F_net, M_net);

		return nets;
		

	}

	bool BuildRobot(Vector2d masses, float d_radius, float arclength, float g, float FI, float FII, float E, int disk_num, float th) {

		Matrix<float, 2, 1> f_tend;
		Matrix<Vector3d, 2, 1> p_tend;
		Vector3d r1(d_radius, 0, 0), r2(-d_radius, 0, 0);

		p_tend << r1, r2; // Creates matrix housing positions of tendons in x-y coordinates
		f_tend << FI, FII; // houses forces of each tendon

		masses = g * masses;
		Matrix<Vector3d, 2, 1> prev_MF;
		Vector3d prev_M(0, 0, 0), prev_F(0, 0, 0);
		prev_MF << prev_M, prev_F;
		//MatrixXd theta = get_theta(); //yeah.... change this to what works
		for (int i = disk_num; i == 0; i--) {

			//float th = theta[i][0];
			VectorXd temp;
			Vector2d currFM = get_segmentfm(p_tend, f_tend, i, disk_num, masses, prev_MF, arclength, th);
			Matrix<float, 3, 3> G_mati = get_transmat(i, 0, arclength, th);
			Disk di(th, arclength, currFM[0], currFM[1], G_mati);
			
			

		}

	}
	
	//sets up parameters for theta solver
	struct theta_params {
		double arclength;
		double young_mod;
		float inertia;
		Matrix<float, 2, 1> p_ten, f_tend;
		int cur_disk;
		Vector2d weights;
		Matrix<Vector3f, 2, 1> prev_MF;
		
	};

	int thetas(const gsl_vector* x, void* p, gsl_vector* f) {

		// I thonk I initizialized this right?? Maybe? I see no errors?
		struct theta_params* params = (struct theta_params*)p;
		const double arclength = (params->arclength);
		const double young_mod = (params->young_mod);
		const float inertia = (params->inertia);
		const Matrix<float, 2, 1> p_ten = (params->p_ten);
		const Matrix<float, 2, 1> f_tend = (params->f_tend);
		const int i = (params->cur_disk);
		const Vector2d weights = (params->weights);
		const Matrix<Vector3f, 2, 1> prev_MF = (params->prev_MF);
		const double theta0 = gsl_vector_get(x, 0);

		const double y0 = theta0 / arclength;
		const double y1 = ((get_segmentfm(p_ten, f_tend, i-1, i, weights, prev_MF, arclength, theta0))[1]).norm();

		gsl_vector_set(f, 0, y0);
		gsl_vector_set(f, 1, y1); //May to actually use "magnitude" thing
		//^alternatively I can actually convert it into a GSL matrix like a normal human being

		//CHECK IF IM SETTING THETA AS A VECTOR OR IF ITS JUST 

		return GSL_SUCCESS;
	};


	int get_theta(float inertia, Matrix<float, 2, 1> p_ten, int cur_disk, Matrix<Vector3f, 2, 1> prev_MF) {
		//cur_disk is the segment that is being worked on. First ieration is cur_disk = last disk
		const gsl_multiroot_fsolver_type* T;
		gsl_multiroot_fsolver* s;

		int status;
		size_t i, iter = 0;


		//number of components for the vector, will solve for each segment then append to an array at the end of each iteration
		//^may affect runtime to append?
		//could expand to n = number of disks, make a n_disk dimension theta vector, must change get_FM to accomodate 
		const size_t n = 1;
		struct theta_params p = { arclength,  E, inertia, p_ten, f_tend, cur_disk, weights, prev_MF};
		
		//Attempted to implement type change solution. I am doing something wrong
		//I swear to Allah that I will figure GSL out or so help me 
		Robot* ptr2 = this;
		auto ptr = [=](const gsl_vector* x, void* p, gsl_vector* f)->int {return ptr2->thetas(x, p, f); };
		gsl_function_pp<decltype(ptr)> Fp(ptr);
		gsl_multiroot_function* f = static_cast<gsl_multiroot_function*>(&Fp);
		
		
		double x_init =  -10.0;
		gsl_vector* x = gsl_vector_alloc(n);

		gsl_vector_set(x, 0, x_init);

		T = gsl_multiroot_fsolver_hybrids;
		s = gsl_multiroot_fsolver_alloc(T, 1);
		//In the example they have a dereferenece operator before f, however f is already a pointer and i got rid of it
		//It works apparently?
		gsl_multiroot_fsolver_set(s, f, x);


		//below is ripped from the GNU website so i can see the thing happen
		print_state (iter, s);

		do {
			iter++;
			status = gsl_multiroot_fsolver_iterate(s);

			print_state(iter, s);

			if (status)   /* check if solver is stuck */
				break;

			status =
				gsl_multiroot_test_residual(s->f, 1e-7);
		} while (status == GSL_CONTINUE && iter < 1000);

		printf("status = %s\n", gsl_strerror(status));

		gsl_multiroot_fsolver_free(s);
		gsl_vector_free(x);
		return 0;
	
	}

	//below is also ripped from the GNU to see solver do its thing
	int print_state(size_t iter, gsl_multiroot_fsolver* s)
	{
		printf("iter = %3u x = % .3f % .3f "
			"f(x) = % .3e % .3e\n",
			iter,
			gsl_vector_get(s->x, 0),
			gsl_vector_get(s->x, 1),
			gsl_vector_get(s->f, 0),
			gsl_vector_get(s->f, 1));
	}
	

};

//template?
template< typename F >
class gsl_function_pp : public gsl_multiroot_function {
public:
	gsl_function_pp(const F& func) : _func(func) {
		function = &gsl_function_pp::invoke;
		params = this;
	}
private:
	const F& _func;
	static double invoke(double x, void* params) {
		return static_cast<gsl_function_pp*>(params)->_func(x);
	}
};

int main() {
	
	//Initialize values
	int ft1 = 2, ft2 = 2, num_d = 1, roh = 6450; // tendon forces and disk number
	float radius_d = 0.5, g = 9.81, E = 54 * 10 ^ 9, l = 0.2;  // radius of disk
	float d_mass = 0.005, t_mass = 0.0005; // mass of each disk and tube segment

	Matrix<Vector3f, 2, 1> prev(0,0);
	//Organizing Matrices
	Matrix<float, 2, 1> r1, r2, p_tend, f_tend;
	Vector2d weights(d_mass*9.81, t_mass*9.81); // matrix intialization
	r1 << radius_d, 0;
	r2 << -radius_d, 0;
	p_tend << r1, r2; // Creates matrix housing positions of tendons in x-y coordinates
	f_tend << ft1, ft2; // houses forces of each tendon


	Robot x(l, radius_d, 9.81, num_d, E, weights, f_tend);

	auto y = x.get_theta(2, p_tend, 1, prev);
	

	


	
	
	






	//timing functions
	auto t1 = std::chrono::high_resolution_clock::now();
	//Place function i want to test here:
	

	auto t2 = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

	std::cout << duration;
	return 0;
	
}



