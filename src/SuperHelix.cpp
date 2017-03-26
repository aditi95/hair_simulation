#include "SuperHelix.hpp"

float deg2rad = 0.0174532925;
Vector3f SuperHelix::getOmega(Vector3f n0, Vector3f n1, Vector3f n2, int segment) {
	return	(q[3*segment] * n0 + 
			 q[3*segment+1] * n1 +
			 q[3*segment+2] * n2).normalized();
}

float SuperHelix::getOmegaNorm (Vector3f n0, Vector3f n1, Vector3f n2, int segment) {
	return	(q[3*segment] * n0 + 
			 q[3*segment+1] * n1 +
			 q[3*segment+2] * n2).norm();
}

// Vector3f SuperHelix::getn(int index, float s, int segment) {

// }

void SuperHelix::calculateLnormals () {

	nL[0][0] = n_initial[0];
	nL[1][0] = n_initial[1];
	nL[2][0] = n_initial[2];

	for(int j=1; j<n; j++) {

		for(int i=0; i<3; i++) {
			
			Vector3f w = getOmega(nL[0][j-1], nL[1][j-1], nL[2][j-1], j-1);
			float omega = getOmegaNorm(nL[0][j-1], nL[1][j-1], nL[2][j-1], j-1);
			Vector3f n_prev_parallel = getparallel(nL[i][j-1], w);
			Vector3f n_prev_perp = getperpendicular(nL[i][j-1], w);

			nL[i][j] =  n_prev_parallel + n_prev_perp * cos(omega * len[j-1] * deg2rad) +
						w.cross(getperpendicular(nL[0][j-1], w)) * sin(omega * len[j-1] * deg2rad);
		
		}
	}
}

Vector3f SuperHelix::getparallel(Vector3f n, Vector3f w) {
	return n.dot(w) * w;
}

Vector3f SuperHelix::getperpendicular(Vector3f n, Vector3f w) {
	return n - (n.dot(w) * w);
}

void SuperHelix::calculateLrs() {

	rL[0] = initial_position;

	for(int j=1; j<n; j++) {
			
		Vector3f w = getOmega(nL[0][j-1], nL[1][j-1], nL[2][j-1], j-1);
		float omega = getOmegaNorm(nL[0][j-1], nL[1][j-1], nL[2][j-1], j-1);

		Vector3f n_prev_parallel = getparallel(nL[0][j-1], w);
		Vector3f n_prev_perp = getperpendicular(nL[0][j-1], w);

		rL[j] = rL[j-1] + n_prev_parallel * len[j-1] + n_prev_perp * sin(omega * len[j-1] * deg2rad) / omega +
					w.cross(getperpendicular(nL[0][j-1], w)) * (1 - cos(omega * len[j-1] * deg2rad)) / omega;
		
	}
}

Vector3f SuperHelix::getr(float s, int segment) {

	Vector3f w = getOmega(nL[0][segment], nL[1][segment], nL[2][segment], segment);
	float omega = getOmegaNorm(nL[0][segment], nL[1][segment], nL[2][segment], segment);

	Vector3f n_parallel = getparallel(nL[0][segment], w);
	Vector3f n_perp = getperpendicular(nL[0][segment], w);

	return rL[segment] + n_parallel * (s - cum_len[segment]) + n_perp * sin(omega * (s - cum_len[segment] * deg2rad)) / omega +
				w.cross(getperpendicular(nL[0][segment], w)) * (1 - cos(omega * (s - cum_len[segment] * deg2rad))) / omega;

}

Vector3f SuperHelix::getJacobian(float s, int segment, int index) {
	Vector3f jacobian;
	float epsilon = 0.01;

	Vector3f oldr = getr(s, segment);

	q[3*segment + index] += epsilon;

	calculateLnormals();
	calculateLrs();

	Vector3f newr = getr(s, segment);

	q[3*segment + index] -= epsilon;

	calculateLnormals();
	calculateLrs();

	return (newr - oldr)/epsilon;
}

MatrixXf SuperHelix::getM() {
	MatrixXf M(3*n, 3*n);

	for(int i=0; i<3*n; i++) {
		for(int j=i; j<3*n; j++) {
			float sum = 0;
			for(int k=0; k<n; k++) {
				float s = cum_len[k];
				while(s <= (cum_len[k]+len[k])) {
					sum += getJacobian(s, (i/3), (i%3)).dot(getJacobian(s, (j/3), (j%3)));
					s += integration_step;
				}
			}
			M(i,j) = M(j,i) = sum;
		}
	}
	return M;
}

MatrixXf SuperHelix::getK() {
	MatrixXf K(3*n, 3*n);

	for(int i=0; i<3*n; i++) {
		for(int j=i+1; j<3*n; j++) {
			K(i,j) = K(j,i) = 0;
		}
		K(i,i) = len[(i/3)] * EI[(i%3)];
	}

	return K;
}

VectorXf SuperHelix::getA() {
	VectorXf A(3*n);
	Vector3f g;
	g[0] = g[2] = 0;
	g[1] = -9.8;
	for(int i=0; i<3*n; i++) {
		float sum = 0;
		for(int k=0; k<n; k++) {
			float s = cum_len[k];
			while(s <= (cum_len[k]+len[k])) {
				sum += getJacobian(s, (i/3), (i%3)).dot(g) * rho;
				s += integration_step;
			}
		}
		A[i] = sum;
	}
	return A;
}

VectorXf SuperHelix::getQint() {
	VectorXf Qint(3*n);
	for(int i=0; i<3*n; i++) {
		Qint[i] = 0;
	}
	return Qint;
}


SuperHelix::SuperHelix() {
	float pi = 3.14;
	n = 3;
	hair_radius = 35e-6;
	rho = 1.3 ;
	integration_step = 5;
	epsilon = 1;

	a1 = hair_radius;
	a2 = hair_radius;
	sigma = 0.48;

	E = 4e9;
	I1 = pi * a1 * a2 * a2 * a2 / 4;
	I2 = pi * a1 * a1 * a1 * a2 / 4;

	mu = E / (2 * (1 + sigma));
	J = pi * a1 * a1 * a1 * a2 * a2 * a2 / (a1 * a1 + a2 * a2);
	EI[0] = mu * J;
	EI[1] = E * I1;
	EI[2] = E * I2;

	M.resize(3*n, 3*n);
	A.resize(3*n);
	K.resize(3*n, 3*n);

	nL.resize(3, std::vector <Vector3f>(n));
	rL.resize(n);

	len.resize(n,0);
	cum_len.resize(n,0);
	len[0] = 50;

	for(int i=1; i<n; i++) {
		len[i] = 50;
		cum_len[i] = len[i] + cum_len[i-1];
	}

	total_len = len[n-1] + cum_len[n-1];

	q.resize(3*n);
	q1.resize(3*n);
	q2.resize(3*n);
	for(int i=0; i<3*n; i++) 
		q1[i] = q2[i] = 0;

	qprev.resize(3*n);
	qprevprev.resize(3*n);

	helix_radius = 1e-2;
	helix_step = 0.5e-2;

	for(int i=0; i<3*n; i++) {

		if(i%3 == 0)
			q[i] = helix_step / (2 * 3.14 * helix_radius * helix_radius);
		else if (i%3 == 1)
			q[i] = 1.0/helix_radius;
		else 
			q[i] = 0;
	}

	initial_position[0] = 0;
	initial_position[1] = 0;
	initial_position[2] = 0;

	n_initial[0] = -1 * Vector3f::UnitY();
	n_initial[1] = Vector3f::UnitX();
	n_initial[2] = Vector3f::UnitZ();
}

void SuperHelix::update() {
	MatrixXf M = getM();
	MatrixXf K = getK();
	MatrixXf A = getA();
	MatrixXf Qint = getQint();

	q1 = (q - qprev)/epsilon;
	q2 = q - 2* qprev + qprevprev;

	MatrixXf AA = M + epsilon * mu * K + epsilon * epsilon * K;
	VectorXf qn = q;
	for(int i=1; i<n; i++) 
		qn[i] *= q[i];
	VectorXf B = epsilon * (M * q1 - epsilon * (A + Qint) + K * (qn - q));
	VectorXf qd = AA.jacobiSvd(ComputeThinU | ComputeThinV).solve(B);
	qprevprev = qprev;
	qprev = q;
	q = q + qd;

	calculateLnormals();
	calculateLrs();

}

void SuperHelix::renderstrand() {
	glPushMatrix(); 
	glBegin(GL_LINE_STRIP);


	for(int i=0; i<n; i++) {
		float start = cum_len[i];
		float end = start + len[i];
		// std::cout << "Segment : " << i << std::endl;
		while (true) {
			if(start > end)
				break;
			glVertex3f((getr(start, i))[0], 2*(getr(start, i))[1], (getr(start, i))[2]);
			// std::cout << (getr(start, i))[0] << " " << (getr(start, i))[1] << " " << (getr(start, i))[2] << std::endl;
			start += 1;
		}

	}

	// glVertex3f(100, 100, 0);
	// glVertex3f(-100, 100, 0);
	// glVertex3f(-100, -100, 0);
	// glVertex3f(100, -100, 0);
	// glVertex3f(100, -100, 100);
	// glVertex3f(-100, -100, 100);
	// glVertex3f(100, 100, 0);
	// glVertex3f(100, 100, 0);
	// glVertex3f(100, 100, 0);

	// glVertex3f(v1[0],v1[1],v1[2]);
	// glVertex3f(v2[0],v2[1],v2[2]);
	glEnd();
	// glBegin(GL_LINE_STRIP);
	// for(int i = 0; i <360; i+=10)
	// {
	// 	glVertex3f(100*cos(i*0.0174532925), 100*sin(i*0.0174532925), 1);
	// }
	// glEnd();
	glPopMatrix(); 

}

void display_func(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0,0,300,0,0,-1,0,1,0);

	// glTranslatef(0.0,-50.0,0.0);
	glPushMatrix();

	hairstrand.update();
	hairstrand.renderstrand(); 
	glutSwapBuffers();

	glutPostRedisplay();
}

void reshape_func(int w, int h)
{
		glViewport(0, 0, (GLsizei) w, (GLsizei) h);
		glMatrixMode(GL_PROJECTION);
		gluPerspective(60.0, 1.0, 1.0, 1000.0);
		glutPostRedisplay();
}

void keyboard_func (unsigned char key, int x, int y)
{
  switch (key)
    {
    case 27:
      exit(0);
      break;
    default:
      break;
    }
}



int main(int argc, char **argv) {
	// std::cout << "Hello" << std::endl;
	hairstrand.calculateLnormals();
	// std::cout << "Hello" << std::endl;
	hairstrand.calculateLrs();
	hairstrand.qprevprev = hairstrand.qprev = hairstrand.q;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(512,512);
	glutInitWindowPosition(300, 300);
	glutCreateWindow("Hair Simulation");
	glutDisplayFunc(display_func);
	glutReshapeFunc(reshape_func);
	glutKeyboardFunc(keyboard_func);
	glutMainLoop();


	return 0;
}
