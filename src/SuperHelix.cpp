#include "SuperHelix.hpp"

float deg2rad = 1;
int nguides = 2;
int nint = 1;

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

	return rL[segment] + n_parallel * (s - cum_len[segment]) + n_perp * sin(omega * (s - cum_len[segment]) * deg2rad) / omega +
				w.cross(getperpendicular(nL[0][segment], w)) * (1 - cos(omega * (s - cum_len[segment]) * deg2rad)) / omega;

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
			float s = cum_len[0];
			while(s <= (cum_len[n-1]+len[n-1])) {
				sum += getJacobian(s, (i/3), (i%3)).dot(getJacobian(s, (j/3), (j%3))) * integration_step;
				s += integration_step;
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
	for(int i=0; i<3*n; i++) {
		float sum = 0;
		float s = cum_len[0];
		while(s <= (cum_len[n-1]+len[n-1])) {
			sum += (getJacobian(s, (i/3), (i%3)).dot(g) * rho * 3.14 * hair_radius * hair_radius) * integration_step;
			s += integration_step;
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

	constantofmult = 100000;
	g[0] = g[2] = 0;
	g[1] = 9.8 * constantofmult;
	n = 3;
	hair_radius = 35e-6 * constantofmult;
	rho = 1.3e-3 ;
	integration_step = 1;
	epsilon = 0.01;
	// epsilon = 1.0 / (rho * pi * hair_radius * hair_radius);

	seglen = 20;
	
	a1 = hair_radius;
	a2 = hair_radius;
	sigma = 0.48;

	E = 2e09 / constantofmult;
	I1 = pi * a1 * a2 * a2 * a2 / 4;
	I2 = pi * a1 * a1 * a1 * a2 / 4;

	mu = E / (2 * (1 + sigma));
	J = pi * a1 * a1 * a1 * a2 * a2 * a2 / (a1 * a1 + a2 * a2);
	EI[0] = mu * J;
	EI[1] = E * I1;
	EI[2] = E * I2;

	helix_radius = 1e-2 * constantofmult;
	helix_step = 0.5e-2 * constantofmult;

	M.resize(3*n, 3*n);
	A.resize(3*n);
	K.resize(3*n, 3*n);

	nL.resize(3, std::vector <Vector3f>(n));
	rL.resize(n);

	len.resize(n,0);
	cum_len.resize(n,0);
	len[0] = seglen;

	for(int i=1; i<n; i++) {
		len[i] = seglen;
		cum_len[i] = len[i] + cum_len[i-1];
	}

	total_len = len[n-1] + cum_len[n-1];

	q.resize(3*n);
	q1.resize(3*n);
	q2.resize(3*n);
	for(int i=0; i<3*n; i++) 
		q1[i] = q2[i] = 0;

	qprev.resize(3*n);
	qrest.resize(3*n);


	// for(int i=0; i<3*n; i++) {
	// 	if(i%3 == 0)
	// 		q[i] = 0.1; //helix_step / (2 * 3.14 * helix_radius * helix_radius);
	// 	else if (i%3 == 1)
	// 		q[i] = 0.1; //1.0/helix_radius;
	// 	else
	// 		q[i] = 0;
	// }

	for(int i=0; i<3*n; i++) {
		if(i%3 == 0)
			qrest[i] = helix_step / (2 * 3.14 * helix_radius * helix_radius);
		else if (i%3 == 1)
			qrest[i] = 1.0/helix_radius;
		else
			qrest[i] = 0;
	}

	q = qrest;
	// qrest = q;

	initial_position[0] = 0;
	initial_position[1] = 0;
	initial_position[2] = 0;

	// n_initial[0] = Vector3f(+1.0/std::sqrt(2), +1.0/std::sqrt(2), 0);
	// n_initial[1] = Vector3f(-1.0/std::sqrt(2), +1.0/std::sqrt(2), 0);

	n_initial[0] = -1 * Vector3f::UnitX();
	n_initial[1] = Vector3f::UnitY();
	n_initial[2] = -Vector3f::UnitZ();

	nrounds = 1;
}

void SuperHelix::update() {
	M = getM();
	A = getA();
	MatrixXf Qint = getQint();

	q1 = (q - qprev)/epsilon;

	MatrixXf AA = M + epsilon * mu * K + epsilon * epsilon * K;

	VectorXf B = epsilon * (M * q1 - epsilon * (A + Qint) + K * (qrest - q));
	VectorXf qd = AA.jacobiSvd(ComputeThinU | ComputeThinV).solve(B);

	qprev = q;
	q = q + qd;

	for(int i=0; i<3*n; i++)
		std::cout << q[i] << " " << qd[i] << std::endl;

	calculateLnormals();
	calculateLrs();

}

void SuperHelix::renderstrand() {
	glPushMatrix(); 
	glBegin(GL_LINE_STRIP);


	for(int i=0; i<pos.size(); i++) {
		glVertex3f(pos[i][0], pos[i][1], pos[i][2]);
	}

	glEnd();
	glPopMatrix(); 

}

void SuperHelix::setpositions() {
	pos.clear();
	for(int i=0; i<n; i++) {
		float start = cum_len[i];
		float end = start + len[i];
		while (true) {
			if(start > end)
				break;
			pos.push_back(getr(start,i));
			start += 1;
		}
	}
}

void renderinterpolatedstrands() {
	glPushMatrix(); 
	for(int j=0; j<nint; j++) {
		for(int i=0; i<(nguides-1); i++) {
			glBegin(GL_LINE_STRIP);
			float x = (j+1)*1.0 /(nint+1);
			Vector3f myPos;
			for(int k=0; k<guidestrands[i].pos.size(); k++) {
				myPos = guidestrands[i].pos[k]*x + guidestrands[i+1].pos[k]*(1-x);
				glVertex3f(myPos[0], myPos[1], myPos[2]);
			}
			glEnd();
		}
	}
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

	for(int i=0; i<nguides; i++) {
		for(int j=0; j<guidestrands[i].nrounds; j++) {
			guidestrands[i].update();
			guidestrands[i].setpositions();
		}
		guidestrands[i].renderstrand(); 
	}

	renderinterpolatedstrands();

	glutSwapBuffers();

	// sleep(1);
	// std::cout << "Hello" << std::endl;
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

// SuperHelix* guidestrands;
// SuperHelix* interpolatedstrands;

int main(int argc, char **argv) {

	guidestrands = new SuperHelix[nguides];
	// interpolatedstrands = new SuperHelix[nint];

	// if(strcmp(argv[1], "straight") == 0) {
	// 	hairstrand.n_initial[0] = -1 * Vector3f::UnitX();
	// 	hairstrand.n_initial[1] = Vector3f::UnitY();
	// 	hairstrand.n_initial[2] = -Vector3f::UnitZ();

	// 	hairstrand.nrounds = 1;
	// }
	// else if(strcmp(argv[1], "curly") == 0) {

	// 	hairstrand.epsilon = 0.01;

	// 	hairstrand.n_initial[0] = -1 * Vector3f::UnitX();
	// 	hairstrand.n_initial[1] = Vector3f::UnitY();
	// 	hairstrand.n_initial[2] = -Vector3f::UnitZ();

	// 	for(int i=0; i<3*hairstrand.n; i++) {
	// 		if(i%3 == 0)
	// 			hairstrand.qrest[i] = 0.01;
	// 		else if (i%3 == 1)
	// 			hairstrand.qrest[i] = 0.1;
	// 		else
	// 			hairstrand.qrest[i] = 0;
	// 	}

	// 	hairstrand.q = hairstrand.qrest;
	// 	hairstrand.seglen = 40;

	// 	hairstrand.len[0] = hairstrand.seglen;

	// 	for(int i=1; i<hairstrand.n; i++) {
	// 		hairstrand.len[i] = hairstrand.seglen;
	// 		hairstrand.cum_len[i] = hairstrand.len[i] + hairstrand.cum_len[i-1];
	// 	}

	// 	hairstrand.constantofmult = 100000;

	// 	hairstrand.total_len = hairstrand.len[hairstrand.n-1] + hairstrand.cum_len[hairstrand.n-1];

	// 	hairstrand.nrounds = 10;


	// }

	guidestrands[1].initial_position = Vector3f(20,20,0);

	for(int i=0; i<nguides; i++) {
		guidestrands[i].calculateLnormals();
		guidestrands[i].calculateLrs();
		guidestrands[i].qprev = guidestrands[i].q;
		guidestrands[i].K = guidestrands[i].getK();
	}


	// for(int i=0; i<nint; i++) {
	// 	interpolatedstrands[i].calculateLnormals();
	// 	interpolatedstrands[i].calculateLrs();
	// 	interpolatedstrands[i].qprev = interpolatedstrands[i].q;
	// 	interpolatedstrands[i].K = interpolatedstrands[i].getK();
	// }
	// std::cout << "Hello" << std::endl;
	hairstrand.calculateLnormals();
	// std::cout << "Hello" << std::endl;
	hairstrand.calculateLrs();
	hairstrand.qprev = hairstrand.q;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(800,800);
	glutInitWindowPosition(800, 800);
	glutCreateWindow("Hair Simulation");
	glutDisplayFunc(display_func);
	glutReshapeFunc(reshape_func);
	glutKeyboardFunc(keyboard_func);
	glutMainLoop();


	return 0;
}
