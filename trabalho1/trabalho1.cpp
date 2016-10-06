/*
	Name: test01.2.Two.Windows.cpp
	Author: Rodrigo Luis de Souza da Silva
	Author: Bruno José Dembogurski
	Date: 27/01/2011
	Last Update: 27/01/2016
	Description: Managing two windows using PixelLab.
*/

#include "pixelLab.h"
#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include <stdio.h>
#include <stdlib.h>
#include "spline.h"

#include <iostream>
#include <vector>

#define CWIDTH 256
#define CHEIGHT 256

#define BOX_SIZE 232

#define DASHED_AMOUNT 4

#define BOX_LEFT (CWIDTH - BOX_SIZE)
#define BOX_RIGHT CWIDTH
#define BOX_TOP 0
#define BOX_BOTTOM BOX_SIZE

#define DASHED_SPACING (BOX_SIZE / DASHED_AMOUNT)

typedef struct {
	double x;
	double y;
} Point;

// Image Objects
PixelLab *img = NULL, *imgOriginal = NULL;
int window1 = 0;
int window2 = 0;
int movingIndex = -1;

std::vector<double> ptsx;
std::vector<double> ptsy;
tk::spline f;

pixel **m;

void idle() {
	glutPostRedisplay();
}

void printPts() {
	for (int i = 0; i < (int) ptsx.size(); i++) {
		std::cout << "(" << ptsx[i] /*<< "," << ptsy[i] */<< ")#";
	}

	std::cout << std::endl;
}

Point screenToBoxCoords(const double x, const double y) {
	return Point {(256 * (x - BOX_LEFT)) / BOX_SIZE,
		BOX_SIZE - (256 * (y - BOX_LEFT)) / BOX_SIZE};
}

Point boxToScreenCoords(const double x, const double y) {
	return Point {BOX_LEFT + ((x / 256) * BOX_SIZE),
		BOX_BOTTOM - ((y / 256) * BOX_SIZE)};
}

int addPoint(double x, double y) {
	Point p = screenToBoxCoords(x, y);

	int i;
	for (i = 0; x > ptsx[i]; i++);

	if (abs(p.x - ptsx[i-1]) > 2) {
		ptsx.insert(ptsx.begin() + i, p.x);
		ptsy.insert(ptsy.begin() + i, p.y);
	}

	return i;
}

void removePoint(int index) {
	ptsx.erase(ptsx.begin() + index);
	ptsy.erase(ptsy.begin() + index);
}

void updateCurve() {
	if (ptsx.size() > 2) f.set_points(ptsx, ptsy);
}

static void display1(void) {
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// reverse 'Y' to allow correct mouse click (0 value is in the upper left corner for mouse functions)
	gluOrtho2D(0, CWIDTH, CHEIGHT, 0); // You could use "glOrtho(0, 64, 256, 0, -5, 5);" instead if you want 3D features

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// draw left vertical color box
	glBegin(GL_QUADS);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex2d(0, 0);
		glVertex2d(16, 0);
		glVertex2d(16, 232);
		glVertex2d(0, 232);
	glEnd();

	glBegin(GL_QUADS);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex2d(1, 1);
		glVertex2d(15, 1);
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex2d(15, 231);
		glVertex2d(1, 231);
	glEnd();

	// draw bottom horizontal color box
		glBegin(GL_QUADS);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex2d(24, 240);
		glVertex2d(24, 256);
		glVertex2d(256, 256);
		glVertex2d(256, 240);
	glEnd();

	glBegin(GL_QUADS);
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex2d(25, 241);
		glVertex2d(25, 255);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex2d(255, 255);
		glVertex2d(255, 241);
	glEnd();

	// draw center control box
	glBegin(GL_QUADS);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex2d(24, 0);
		glVertex2d(256, 0);
		glVertex2d(256, 232);
		glVertex2d(24, 232);
	glEnd();

	glBegin(GL_QUADS);
		glColor3f(0.94f, 0.94f, 0.94f);
		glVertex2d(25, 1);
		glVertex2d(255, 1);
		glVertex2d(255, 231);
		glVertex2d(25, 231);
	glEnd();

	// draw dashed lines
	glEnable(GL_LINE_STIPPLE);
	glColor3f(0.0f, 0.0f, 0.0f);
	glLineStipple(2, 0xAAAA);
	glLineWidth(1);

	// vertical dashed lines
	for (int i = 1; i < DASHED_AMOUNT; ++i) {
		glBegin(GL_LINES);
			glVertex2f(BOX_LEFT + (i * DASHED_SPACING), BOX_TOP);
			glVertex2f(BOX_LEFT + (i * DASHED_SPACING), BOX_BOTTOM);
		glEnd();
	}

	// horizontal dashed lines
	for (int i = 1; i < DASHED_AMOUNT; ++i) {
		glBegin(GL_LINES);
			glVertex2f(BOX_LEFT, BOX_TOP + (i * DASHED_SPACING));
			glVertex2f(BOX_RIGHT, BOX_TOP + (i * DASHED_SPACING));
		glEnd();
	}

	glDisable(GL_LINE_STIPPLE);

	// draw curve
	glBegin(GL_LINE_STRIP);
		glColor3f(0.0f, 0.0f, 0.0f);

		if (ptsx.size() > 2) {
			for (int x = 0; x < 256; ++x) {
				Point p = boxToScreenCoords(x, f(x));
				glVertex2f(p.x, p.y);
			}
		} else {
			for (int x = 0; x < 256; ++x) {
				Point p;
				if (x <= ptsx[0]) {
					p = boxToScreenCoords(x, 0.0);
				} else if (x >= ptsx[ptsx.size() - 1]) {
					p = boxToScreenCoords(x, 256.0);
				} else {
					p = boxToScreenCoords(x,
						ptsy[0] + ((x - ptsx[0]) / ptsx[1]) * (ptsy[1] - ptsy[0]));
				}

				glVertex2f(p.x, p.y);
			}
		}
	glEnd();

	// draw points
	for (int i = 0; i < (int) ptsx.size(); i++) {
		Point p = boxToScreenCoords(ptsx[i], ptsy[i]);

		glBegin(GL_LINE_STRIP);
			if (movingIndex == i) glColor3f(1.0f, 0.0f, 0.0f);
			else glColor3f(0.0f, 0.0f, 0.0f);

			glVertex2f(p.x - 2, p.y - 2);
			glVertex2f(p.x + 2, p.y - 2);
			glVertex2f(p.x + 2, p.y + 2);
			glVertex2f(p.x - 2, p.y + 2);
			glVertex2f(p.x - 2, p.y - 2);
		glEnd();
	}

	glutSwapBuffers();
}

static void display2(void) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	img->ViewImage();
	glutSwapBuffers();
}

void modifyImage() {
	// Get Image Data
	imgOriginal->GetDataAsMatrix(m);

	int aux;
	if (img) {
		for (int y = 0; y < img->GetHeight(); y++) {
			for (int x = 0; x < img->GetWidth(); x++) {
				 aux = f(m[y][x].value);

				 if (aux < 0) m[y][x].value = 0;
				 else if (aux > 255) m[y][x].value = 255;
				 else m[y][x].value = aux;
			}
		}
	}
	// Update PixelLab Object
	img->SetDataAsMatrix(m);

	// Update both windows
	glutPostWindowRedisplay(window1);
	glutPostWindowRedisplay(window2);
}

// Keyboard
static void key(unsigned char key, int x, int y) {
	 switch (key) {
		case 'q':
			exit(0);
		break;
	 }

	 glutPostRedisplay();
}

int isWithinBounds(int x, int y, int cx, int cy, int radio) {
	if (x < cx + radio && x > cx - radio &&
			y < cy + radio && y > cy - radio) {
		return 1;
	}

	return 0;
}

// Mouse callback - Capture mouse click in the brightness window
void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON) {
		if (state == GLUT_DOWN) {
			if (isWithinBounds(x, y, 140, 116, 116)) {
				for (int i = 0; i < (int) ptsx.size(); i++) {
					if (isWithinBounds(x, y, ptsx[i], ptsy[i], 20)) {
						movingIndex = i;

						return;
					}
				}

				addPoint(x, y);

				updateCurve();
				modifyImage();
			}
		} else if (state == GLUT_UP) {
			movingIndex = -1;
		}
	} /*else if (button == GLUT_RIGHT_BUTTON) {
		if (state == GLUT_DOWN) {
			if (ptsx.size() > 2) {
				for (int i = 0; i < (int) ptsx.size(); i++) {
					if (isWithinBounds(x, y, ptsx[i], ptsy[i], 20)) {
						ptsx.erase(ptsx.begin() + i);
						ptsy.erase(ptsy.begin() + i);

						updateCurve();
						modifyImage();
					}
				}
			}
		}
	}*/
}

// Motion callback - Capture mouse motion when left button is clicked
void motion(int x, int y ) {
	if (movingIndex > 0 && movingIndex < (int) ptsx.size() - 1) {
		if (x > ptsx[movingIndex - 1] &&
				x < ptsx[movingIndex + 1]) {
			removePoint(movingIndex);
			movingIndex = addPoint(x, y);

			updateCurve();
			modifyImage();
		} else {
			movingIndex = -1;
		}
	}
}

void init() {
	imgOriginal = new PixelLab();
	imgOriginal->Read("figs/lenaGray.png");

	ptsx.push_back(0.0);
	ptsy.push_back(0.0);

	ptsx.push_back(128.0);
	ptsy.push_back(128.0);

	ptsx.push_back(256.0);
	ptsy.push_back(256.0);

	updateCurve();

	img = new PixelLab();
	img->Copy(imgOriginal);

	// Allocate Matrix
	img->AllocatePixelMatrix(&m, img->GetHeight(), img->GetWidth());
}

int main(int argc, char *argv[]) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

	// Init image objects
	init();

	// Create two windows
	window2 = glutCreateWindow("Display Window");
	window1 = glutCreateWindow("Control Window");

	///////////////////////////////////////////////////////////
	// Control Window
	glutSetWindow(window1);	// Change current window to 1 (all callbacks will be set to this window)
	glutDisplayFunc(display1);
	glutPositionWindow(20, 30);
	glutReshapeWindow(CWIDTH, CHEIGHT);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutKeyboardFunc(key);

	///////////////////////////////////////////////////////////
	// Display Window
	glutSetWindow(window2); // Change current window to 2
	glutDisplayFunc(display2);
	glutReshapeWindow(imgOriginal->GetWidth(),imgOriginal->GetHeight());
	glutPositionWindow(25 + CWIDTH, 30);
	glutKeyboardFunc(key);

	glutMainLoop();

	// Deallocate Matrix
	img->DeallocatePixelMatrix(&m, img->GetHeight(), img->GetWidth() );

	return 0;
}
