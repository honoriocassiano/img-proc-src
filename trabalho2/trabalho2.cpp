#include "pixelLab.h"
#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include <stdio.h>
#include <stdlib.h>

#include "Debug.h"

#include <iostream>
#include <vector>

#define ZOOM_STEP 0.005
#define MOVE_STEP 1

int width1 = 0;
int height1 = 0;
double posx1 = 0.0;
double posy1 = 0.0;
double factor1 = 1.0;

int width2 = 0;
int height2 = 0;
double posx2 = 0.0;
double posy2 = 0.0;
double factor2 = 1.0;

int width3 = 0;
int height3 = 0;

typedef struct {
	int x;
	int y;
} Point;

// Image Objects
PixelLab *img1 = NULL, *img2 = NULL, *img3 = NULL;

pixel **mat1, **mat2, **mat3;

int window1 = 0;
int window2 = 0;
int window3 = 0;

std::string currentImg1, currentImg2;

Point lastClick;
bool moving = false;

void idle() {
	glutPostRedisplay();
}

static void display1() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0, 0, 0, 0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-100.0, 100.0, -100.0, 100.0, -1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glPushMatrix();

	// View first image
	//float factor = std::min(width1 / ((float) img1->GetWidth()),
	//		height1 / ((float) img1->GetHeight()));

	//float x = -1 + ((width1  - (img1->GetWidth()) * factor)) / (height1 * 2.0f);
	//float y = ((width1 - (img1->GetHeight()) * factor)) / (height1 * 2.0f);

	//std::cout << "img1: x=" << x << ",y=" << y << std::endl;

	//glRasterPos2f(x, y);
	glTranslatef(posx1, posy1, 0.0);
	glRasterPos2f(-100, -100);
	glPixelZoom(factor1, factor1);

	img1->ViewImage();

	glPopMatrix();

	glutSwapBuffers();
}

static void display2() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(1, 1, 1, 1);
	glLoadIdentity();

	glPushMatrix();

	// View second image
	//float factor = std::min(width2 / ((float) img2->GetWidth()),
	//		height2 / ((float) img2->GetHeight()));

	//float x = ((width2  - (img2->GetWidth()) * factor * 2) / 2.0f) / (width2 * 2.0f);
	//float y = ((height2  - (img2->GetWidth()) * factor * 2) / 2.0f) / (height2 * 2.0f);

	//std::cout << "img2: x=" << x << ",y=" << y << std::endl;

	//glRasterPos2f(x, y);
	//glRasterPos2f(posx2, posy2);
	glPixelZoom(factor2, factor2);

	img2->ViewImage();

	glPopMatrix();

	glutSwapBuffers();
}

static void display3() {
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);



	glutSwapBuffers();
}

void modifyImage() {
	/*int value = brightnessSlider - 128; // To make the middle of the slider, zero
	int pixelVal;

	for (int y = 0; y < imgOriginal->GetHeight(); y++) {
		for (int x = 0; x < imgOriginal->GetWidth(); x++) {
			pixelVal = imgOriginal->GetGrayValue(x, y);
			pixelVal += value;
			if (pixelVal >= 255)
				pixelVal = 255;
			if (pixelVal < 0)
				pixelVal = 0;
			img->SetGrayValue(x, y, (uByte) pixelVal);
		}
	}*/

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

// Special Keys callback
void specialKeys1(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_UP:
		factor1 += ZOOM_STEP;
		modifyImage();
		break;

	case GLUT_KEY_DOWN:
		factor1 -= ZOOM_STEP;
		modifyImage();
		break;

	case 'q':
		exit(0);
	}

	glutPostRedisplay();
}

// Special Keys callback
void specialKeys3(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_UP:
		//factor1 += ZOOM_STEP;
		modifyImage();
		break;

	case GLUT_KEY_DOWN:
		//factor1 -= ZOOM_STEP;
		modifyImage();
		break;

	case 'q':
		exit(0);
	}

	glutPostRedisplay();
}

void specialKeys2(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_UP:
		factor2 += ZOOM_STEP;
		modifyImage();
		break;

	case GLUT_KEY_DOWN:
		factor2 -= ZOOM_STEP;
		modifyImage();
		break;

	case 'q':
		exit(0);
	}

	glutPostRedisplay();
}

// Mouse callback - Capture mouse click
void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON) {
		if (state == GLUT_DOWN)
			lastClick = {x, y};
		else
			moving = false;
	}

	modifyImage();
}

// Motion callback - Capture mouse motion
void motion1(int x, int y) {
	if (moving) {
		posx1 += (x - lastClick.x) * MOVE_STEP;
		posy1 += (lastClick.y - y) * MOVE_STEP;

		lastClick = {x, y};
	} else {
		lastClick = {x, y};
		moving = true;
	}

	modifyImage();
}

void init() {
	img1 = new PixelLab();
	img1->Read("figs/lenaGray.png");
	currentImg1 = "figs/lenaGray.png";

	width1 = img1->GetWidth();
	height1 = img1->GetHeight();

	img2 = new PixelLab();
	img2->Read("figs/gearsGray.png");
	currentImg2 = "figs/gearsGray.png";

	width2 = img2->GetWidth();
	height2 = img2->GetHeight();

	img3 = new PixelLab();

	width3 = std::max(img1->GetWidth(), img2->GetWidth());
	height3 = std::max(img2->GetWidth(), img2->GetWidth());

	img3->CreateImage(width3, height3, 3);
	img3->AllocatePixelMatrix(&mat3, width3, height3);

	// Creates a gradient
	for(int y = 0; y < img3->GetHeight(); y++) {
		for(int x = 0; x < img3->GetWidth(); x++) {
			mat3[y][x].R = 255;
			mat3[y][x].G = (y > 255) ? 255 : y;
			mat3[y][x].B = 0;
		}
	}

	// Set Data as a Matrix
	img3->SetDataAsMatrix(mat3);

	img3->ViewImage();

	img3->DeallocatePixelMatrix(&mat3, width3, height3);
}

int main(int argc, char *argv[]) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

	// Init image objects
	init();

	// Create three windows
	std::string title1 = "Image 1: " + currentImg1;
	std::string title2 = "Image 2: " + currentImg2;
	std::string title3 = "Resulting image";

	window1 = glutCreateWindow(title1.c_str());
	window2 = glutCreateWindow(title2.c_str());
	window3 = glutCreateWindow(title3.c_str());

	// Display Window 1
	glutSetWindow(window1); // Change current window to 1
	glutDisplayFunc(display1);
	glutReshapeWindow(img1->GetWidth(), img1->GetHeight());
	glutPositionWindow(50, 50);
	glutKeyboardFunc(key);
	glutSpecialFunc(specialKeys1);
	glutMotionFunc(motion1);

	// Display Window 2
	glutSetWindow(window2); // Change current window to 2
	glutDisplayFunc(display2);
	glutReshapeWindow(img2->GetWidth(), img2->GetHeight());
	glutPositionWindow(100 + img1->GetWidth(), 50);
	glutKeyboardFunc(key);
	glutSpecialFunc(specialKeys2);

	// Display Window 2
	glutSetWindow(window3); // Change current window to 3
	glutDisplayFunc(display3);
	glutReshapeWindow(
			std::max(img1->GetWidth(), img2->GetWidth()),
			std::max(img1->GetHeight(), img2->GetHeight())
	);
	glutPositionWindow(150 + img1->GetWidth() + img2->GetWidth(), 50);
	glutKeyboardFunc(key);
	glutSpecialFunc(specialKeys3);

	glutMainLoop();

	return 0;
}
