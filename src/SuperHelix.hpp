#ifndef __SUPERHELIX_HPP__
#define __SUPERHELIX_HPP__

#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <unistd.h>


#ifdef __APPLE__
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
#endif

using namespace Eigen;


class SuperHelix {
public:
	// Number of segments
	int n;				
	float hair_radius;
	float rho;
	float total_len;
	float integration_step;
	float E;
	float a1,a2;
	float I1, I2;
	float mu;
	float J;
	float EI[3];
	float sigma;
	float epsilon;
	float helix_radius;
	float helix_step;
	float constantofmult;
	Vector3f g;

	// Position of the left end of the segment
	std::vector< Vector3f> L;

	std::vector< std::vector <Vector3f> > nL; 
	std::vector< Vector3f > rL; // INitialise size in a constrctie
	std::vector< float > len;
	std::vector< float > cum_len;

	VectorXf q;
	VectorXf qprev;
	VectorXf qrest;
	VectorXf q1;
	VectorXf q2;
	Vector3f initial_position;
	Vector3f n_initial[3];
	float seglen;
	int nrounds;

	MatrixXf M;
	MatrixXf K;
	VectorXf A;
	SuperHelix();

	Vector3f getOmega(Vector3f n0, Vector3f n1, Vector3f n2, int segment);
	float getOmegaNorm(Vector3f n0, Vector3f n1, Vector3f n2, int segment);
	// Vector3f getn(int index, float s, int segment);
	void calculateLnormals ();
	void calculateLrs();
	Vector3f getparallel(Vector3f n, Vector3f w);
	Vector3f getperpendicular(Vector3f n, Vector3f w);

	Vector3f getr(float s, int segment);

	Vector3f getJacobian(float s, int segment, int index);
	MatrixXf getM();
	MatrixXf getK();
	VectorXf getA();
	VectorXf getQint();

	void update();
	void renderstrand();

};

SuperHelix hairstrand;

#endif