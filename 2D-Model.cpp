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


	Vector3d coord_get(float arclength, float th, Matrix<float, 3, 3> G_mat) {

		float x = (arclength / th) * sin(th);
		float y = (arclength / th) * (1 - cos(th));

		Vector3d coords(x, y, 0);

	}



};


//template?
//template< typename F >
//class gsl_multiroot_function_pp : public gsl_multiroot_function {
//public:
//	gsl_multiroot_function_pp(const F& func) : _func(func) {
//		func = &gsl_multiroot_function_pp::invoke;
//		params = this;
//	}
//private:
//	const F& _func;
//	static double invoke(double x, void* params) {
//		return static_cast<gsl_multiroot_function_pp*>(params)->_func(x);
//	}
//};

class Robot {

public:

	float d_radius, arclength, theta, g, E; // E is youngs modulus
	int disk_num;
	Matrix<Disk, 1, Dynamic> disks;
	Matrix<Vector3f, 2, 1> weights;
	Matrix<float, 2, 1> f_tend;


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

	Matrix<Vector3f, 2, 1> get_segmentfm(Matrix<Vector3f, 2, 1> p_ten, Matrix<float, 2, 1> f_tend, int i, int n, Matrix<Vector3f, 2, 1> weights, Matrix<Vector3f, 2, 1> prev_MF, float l, float th) {
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

		Vector3f Fg_d = G_mat * weights[0]; //force g on disk i in reference of i-1 MAKE SURE TO DEFINE MASSES[0] AS [0, mass*g, 0]
		Vector3f Mg_d = prev_coord.cross(Fg_d); //moment of Fg_d

		Vector3f e = prev_MF[1];

		Vector3f MFi = prev_coord.cross(e);

		Vector3f M1_rel = ai.cross(F1_rel);
		Vector3f M2_rel = bi.cross(F2_rel);

		Vector3f prevF = T_mat * prev_MF[1]; //puts previous force in terms of i-1 reference frame
		Vector3f prevM = T_mat * prev_MF[0]; //previous moment in terms of i-1 reference frame

		Vector3f F_net = F1_rel + F2_rel + prevF + Fg_d;
		Vector3f M_net = M1_rel + M2_rel + prevM + Mg_d; // + MFi  add this when you fix the cross product thing


		Matrix<Vector3f, 2, 1> nets;

		nets << F_net, M_net;

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


	}




};

int print_state(size_t iter, gsl_multiroot_fsolver* s)
{
	printf("iter = %3u x = % .3f % .3f "
		"f(x) = % .3e % .3e\n",
		iter,
		gsl_vector_get(s->x, 0),
		gsl_vector_get(s->x, 1),
		gsl_vector_get(s->f, 0),
		gsl_vector_get(s->f, 1));
	//gsl_vector_get(s->x, 1),

	return 0;
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


MatrixXd get_segmentfm(Matrix<Vector3f, 2, 1> p_ten, Matrix<float, 2, 1> f_tend, int n, Matrix<Vector3f, 2, 1> weights, float l, MatrixXd th, float E, float I) {
	//gets current force and moment on each tendon 
	MatrixXd thetas;
	thetas.resize(n, 1);

	for (int i = n; i == 1; i - 1) {
		Matrix<float, 3, 3> T_mat = get_transmat(i, i - 1, l, th(i,0)); //figure out how to get rid of the need for an l
		Vector3f prev_coord(T_mat(2, 0), T_mat(2, 1), 0);

		Vector3f ai = T_mat * (p_ten[0]);
		Vector3f bi = T_mat * (p_ten[1]); //p_ten houses tendon locations

		//figure out how to combine ai and bi 

		ai = (-1 * ai) + p_ten(0);  //ai to ai-1
		bi = (-1 * bi) + p_ten(1);

		Vector3f F1_rel; //combine into one matrix to allow for one operation 
		Vector3f F2_rel;
		Matrix<Vector3f, 2, 1> prev_MF;

		if (i == n) { //i = segment #, n = total disks



			F1_rel = (f_tend[0] * (ai)) / ai.norm();
			F2_rel = (f_tend[1] * (bi)) / (bi).norm();
			Vector3f M(0, 0, 0), F(0, 0, 0);
			prev_MF << M, F;

		}

		else if (i < n) {
			Matrix<float, 3, 3> T_mat1 = get_transmat(i + 1, i, l, th(i,0)); // i+1 to i-1
			T_mat1 = T_mat * T_mat1;

			Vector3f ai_1 = T_mat1 * (p_ten[0]); // instead of i+1 to i-1, reuse i to i-1 for claculations
			Vector3f bi_1 = T_mat1 * (p_ten[1]); //p_ten houses tendon locations

			ai_1 = ai_1 - ai;  //i to i+1 vector
			bi_1 = bi_1 + bi;


			F1_rel = f_tend[0] * ((ai / ai.norm()) + (ai_1 / ai_1.norm()));
			F2_rel = f_tend[1] * ((bi / bi.norm()) + (bi_1 / bi_1.norm()));


		}

		Matrix<float, 3, 3> G_mat = get_transmat(i - 1, 0, l, th(i,0));

		Vector3f Fg_d = G_mat * weights[0]; //force g on disk i in reference of i-1 MAKE SURE TO DEFINE MASSES[0] AS [0, mass*g, 0]
		Vector3f Mg_d = prev_coord.cross(Fg_d); //moment of Fg_d

		Vector3f e = prev_MF[1];

		Vector3f MFi = prev_coord.cross(e);

		Vector3f M1_rel = ai.cross(F1_rel);
		Vector3f M2_rel = bi.cross(F2_rel);

		Vector3f prevF = T_mat * prev_MF[1]; //puts previous force in terms of i-1 reference frame
		Vector3f prevM = T_mat * prev_MF[0]; //previous moment in terms of i-1 reference frame

		Vector3f F_net = F1_rel + F2_rel + prevF + Fg_d;
		Vector3f M_net = M1_rel + M2_rel + prevM + Mg_d; // + MFi  add this when you fix the cross product thing

		thetas(i, 0) = M_net[1] * l / (E * I) - th(i,0);

		prev_MF << M_net, F_net;

		
	}
	return thetas;
}

struct theta_params {
	double arclength;
	double young_mod;
	float inertia;
	Matrix<float, 2, 1> f_tend;
	int n;
	Matrix<Vector3f, 2, 1> prev_MF, p_ten, weights;

};

int thetas(const gsl_vector* x, void* p, gsl_vector* f) {

	// I thonk I initizialized this right?? Maybe? I see no errors?
	struct theta_params* params = (struct theta_params*)p;
	const double arclength = (params->arclength);
	const double young_mod = (params->young_mod);
	const float inertia = (params->inertia);
	const Matrix<Vector3f, 2, 1> p_ten = (params->p_ten);
	const Matrix<float, 2, 1> f_tend = (params->f_tend);
	const int n = (params->n);
	const Matrix<Vector3f, 2, 1> weights = (params->weights);
	const Matrix<Vector3f, 2, 1> prev_MF = (params->prev_MF);

	MatrixXd th;
	for (int i = 0; i < n; i++)
	{
		th(i, 0) = gsl_vector_get(x, i);
	}

	MatrixXd expressions = get_segmentfm(p_ten, f_tend, n, weights, arclength, th, young_mod, inertia);
	double exp = 0;
	for (int i = 0; i < n; i++)
	{
		exp = exp + expressions(i, 0);
	}
	gsl_vector_set(f, 0, exp);


	return GSL_SUCCESS;
};

int main() {
	//Initialize values

	Matrix<Vector3f, 2, 1> p_ten, weights, prev_MF;
	Matrix<float, 2, 1> f_tend;
	f_tend << 5, 0;
	Vector3f pM(0, 0, 0), pF(0, 0, 0);

	prev_MF << pM, pF;
	Vector3f d_weight(0, -9.81 * 0.05, 0), t_weight(0, 0, 0);
	weights << d_weight, t_weight;

	int disks = 1, r = 0.001, l = 0.02, E = 54 * 10 ^ 9, I = 0.6;



	//Matrix<Vector3f, 2,1> x = get_segmentfm(p_ten, f_tend, 1, 1, weights, prev_MF, l, 0.0981);

	//std::cout << x[0], x[1];


	//timing functions
	auto t1 = std::chrono::high_resolution_clock::now();
	//Place function i want to test here:

	const gsl_multiroot_fsolver_type* T;
	gsl_multiroot_fsolver* s;

	int status;
	size_t i, iter = 0;

	const size_t n = disks;
	struct theta_params p = { l, E, I, f_tend, disks, prev_MF, p_ten, weights };
	gsl_multiroot_function f = { &thetas, n, &p };

	double x_init[1] = { 0.1 };
	gsl_vector* x = gsl_vector_alloc(n);

	
	gsl_vector_set(x, 0, x_init[0]);
	
	

	
	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc(T, n);
	gsl_multiroot_fsolver_set(s, &f, x); //here

	print_state(iter, s);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);

		print_state(iter, s);

		if (status)   //check if solver is stuck 
			break;

		status =
			gsl_multiroot_test_residual(s->f, 1e-7);
	} while (status == GSL_CONTINUE && iter < 1000);

	printf("status = %s\n", gsl_strerror(status));

	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);


	auto t2 = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

	std::cout << duration;

	return 0;

}