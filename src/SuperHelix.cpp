#include "SuperHelix.hpp"

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

			nL[i][j] =  n_prev_parallel + n_prev_perp * cos(omega * len[j-1]) + 
						w.cross(getperpendicular(nL[0][j-1], w)) * sin(omega * len[j-1]);
		
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

		rL[j] = rL[j-1] + n_prev_parallel * len[j-1] + n_prev_perp * sin(omega * len[j-1]) / omega + 
					w.cross(getperpendicular(nL[0][j-1], w)) * (1 - cos(omega * len[j-1])) / omega;
		
	}
}

Vector3f SuperHelix::getr(float s, int segment) {

	Vector3f w = getOmega(nL[0][segment], nL[1][segment], nL[2][segment], segment);
	float omega = getOmegaNorm(nL[0][segment], nL[1][segment], nL[2][segment], segment);

	Vector3f n_parallel = getparallel(nL[0][segment], w);
	Vector3f n_perp = getperpendicular(nL[0][segment], w);

	return rL[segment] + n_parallel * (s - cum_len[segment]) + n_perp * sin(omega * (s - cum_len[segment])) / omega + 
				w.cross(getperpendicular(nL[0][segment], w)) * (1 - cos(omega * (s - cum_len[segment]))) / omega;

}

SuperHelix::SuperHelix() {
	n = 3;

	nL.resize(3, std::vector <Vector3f>(n));
	rL.resize(n);

	len.resize(n,0);
	cum_len.resize(n,0);
	len[0] = 50;

	for(int i=1; i<n; i++) {
		len[i] = 50;
		cum_len[i] = len[i] + cum_len[i-1];
	}

	q.resize(3*n);

	for(int i=0; i<3*n; i++) {
		if(i%3 == 0)
			q[i] = 10;
		else if (i%3 == 1)
			q[i] = 0;
		else 
			q[i] = 0;
	}

	initial_position[0] = 0;
	initial_position[1] = 0;
	initial_position[2] = 0;

	n_initial[0] = Vector3f::UnitX();
	n_initial[1] = Vector3f::UnitY();
	n_initial[2] = Vector3f::UnitZ();
}

void renderstrand() {
	glPushMatrix(); 
	glBegin(GL_LINE_STRIP);


	for(int i=0; i<hairstrand.n; i++) {
		float start = hairstrand.cum_len[i];
		float end = start + hairstrand.len[i];
		std::cout << "Segment : " << i << std::endl;
		while (true) {
			if(start > end)
				break;
			glVertex3f((hairstrand.getr(start, i))[0], 2*(hairstrand.getr(start, i))[1], (hairstrand.getr(start, i))[2]);
			std::cout << (hairstrand.getr(start, i))[0] << " " << (hairstrand.getr(start, i))[1] << " " << (hairstrand.getr(start, i))[2] << std::endl;
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

  renderstrand(); 
  
  glutSwapBuffers();
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