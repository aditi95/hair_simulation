#ifndef _SUPERHELIX_HPP_
#define _SUPERHELIX_HPP_

#include <bits/stdc++.h>
#include <Eigen/Dense>

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
	// Position of the left end of the segment
	std::vector< Vector3f> L;

	std::vector< std::vector <Vector3f> > nL; 
	std::vector< Vector3f > rL; // INitialise size in a constrctie
	std::vector< float > len;
	std::vector< float > cum_len;

	VectorXf q;
	Vector3f initial_position;
	Vector3f n_initial[3];

	MatrixXf jacobian;

	Vector3f getOmega(Vector3f n0, Vector3f n1, Vector3f n2, int segment);
	float getOmegaNorm(Vector3f n0, Vector3f n1, Vector3f n2, int segment);
	// Vector3f getn(int index, float s, int segment);
	void calculateLnormals ();
	void calculateLrs();
	Vector3f getparallel(Vector3f n, Vector3f w);
	Vector3f getperpendicular(Vector3f n, Vector3f w);

	Vector3f getr(float s, int segment);

	SuperHelix();

};

SuperHelix hairstrand;

#endif
