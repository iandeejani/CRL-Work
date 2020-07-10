// Frickin 5th Project File Lets Make this Work.cpp : This file contains the 'main' function. Program execution begins and ends there.
//



#include <iostream>
#include <cmath>
#include <Eigen/Dense>


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
	Disk(float th, float alarclength, float Mnet, float Fnet, Matrix<float, 3, 3> G_mat) :theta(th), arclength(arclength), Mnet(Mnet), Fnet(Fnet) {}
	
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
	Matrix<Disk, Dynamic, 1> disks;
	Vector2d masses;
	Robot(float FI, float FII, float l, float r, float g, int d, float E, Vector2d masses) :FI(FI), FII(FII), d_radius(r), arclength(l), disk_num(d), g(g), E(E), masses(masses) {
	
		
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

	Vector2d get_segmentfm(Matrix<Vector3d, 2, 1> p_ten, Matrix<float, 2, 1> f_tend, int i, int n, Matrix<float, 2, 1> weights, Matrix<Vector3d, 2, 1> prev_MF, float l, float th) {
		//gets current force and moment on each tendon 

		Matrix<float, 3, 3> T_mat = get_transmat(i, i-1, l, th); //figure out how to get rid of the need for an l
		Vector3d prev_coord(T_mat(2,0), T_mat(2,1), 0);

		Vector3d ai = T_mat * (p_ten[0]);
		Vector3d bi = T_mat * (p_ten[1]); //p_ten houses tendon locations

		ai = (-1 * ai) + p_ten(0);  //ai to ai-1
		bi = (-1 *bi) + p_ten(1);
		
		Vector3d F1_rel;
		Vector3d F2_rel;

		if (i == n) { //i = segment #, n = total disks


			Vector3d F1_rel = (f_tend[0] * (ai)) / ai.norm;
			Vector3d F2_rel = (f_tend[1] * (bi)) / (bi).norm;

		}

		else if (i < n) {
			Matrix<float, 3, 3> T_mat1 = get_transmat(i + 1, i - 1, l, th); // i+1 to i-1
			
			Matrix<float, 3, 1> ai_1 = T_mat1 * (p_ten[0]);
			Matrix<float, 3, 1> bi_1 = T_mat1 * (p_ten[1]); //p_ten houses tendon locations
			
			Vector3d ai_1 = ai_1 - ai;  //i to i+1 vector
			Vector3d bi_1 = bi_1 + bi;


			F1_rel = f_tend[0] * ((ai / ai.norm) + (ai_1 / ai_1.norm));  
			F2_rel = f_tend[1] * ((bi / bi.norm) + (bi_1 / bi_1.norm));

		}

		Matrix<float, 3, 3> G_mat = get_transmat(i-1, 0, l, th);
		
		Matrix<float, 3, 1> Fg_d = G_mat * weights[0]; //force g on disk i in reference of i-1 MAKE SURE TO DEFINE MASSES[0] AS [0, mass*g, 0]
		Matrix<float, 3, 1> Mg_d = prev_coord.cross(Fg_d); //moment of Fg_d

		Matrix<float, 3, 1>	MFi = prev_coord.cross(prev_MF[1]);

		Vector3d M1_rel = ai.cross(F1_rel);
		Vector3d M2_rel = bi.cross(F2_rel);

		prev_MF[1] = T_mat * prev_MF[1]; //puts previous force in terms of i-1 reference frame
		prev_MF[0] = T_mat * prev_MF[0]; //previous moment in terms of i-1 reference frame

		Vector3d F_net = F1_rel + F2_rel + prev_MF[1] + Fg_d;
		Vector3d M_net = M1_rel + M2_rel + prev_MF[0] + Mg_d + MFi;
		
		Vector2d nets(M_net, F_net);

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
		for (int i = disk_num; i = 0; i--) {

			Vector2d currFM = get_segmentfm(p_tend, f_tend, i, disk_num, masses, prev_MF, arclength, th);
			
			

		}

	}



};




int main() {
	/*
	//Initialize values
	int ft1 = 2, ft2 = 2, num_d = 1, roh = 6450; // tendon forces and disk number
	float radius_d = 0.5, g = 9.81, E = 54 * 10 ^ 9;  // radius of disk
	float d_mass = 0.005, t_mass = 0.0005; // mass of each disk and tube segment

	//Organizing Matrices
	Matrix<float, 2, 1> r1, r2, p_tend, f_tend;
	Vector2d masses(d_mass, t_mass); // matrix intialization
	r1 << radius_d, 0;
	r2 << -radius_d, 0;
	p_tend << r1, r2; // Creates matrix housing positions of tendons in x-y coordinates
	f_tend << ft1, ft2; // houses forces of each tendon





	Matrix<float, 3, 3> x;
	

	std::cout << x;


	Vector3d gd(0, -9.81 * masses(0), 3), gt(0, -9.81 * masses(1), 3);
	std::cout << gd;
	std::cout << gt;
	*/

	Matrix<float, 3, 3> asf;
	asf << 1, 2, 1, 1, 2, 1, 0, 0, 1;


	Disk bob(2, 2, 2, 2, asf);

	bob.coord_get(2, 2, asf);
}



// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
