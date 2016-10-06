#include "pixelLab.h"
#ifdef WIN32
#include <windows.h>
#endif
#include <iostream>
#include <cmath>
#include "GFFT.h"
#include "GL/glui.h"
#include "GL/glut.h"

// Image Objects
PixelLab *img = NULL;
PixelLab *imgMod = NULL;

GLUI_RadioGroup *radio, *radio2;
GLUI_RadioButton *buttons[3];
int main_window;
int UI_width = 0;
int option = 0;
int option2 = 0;
int brightness = 0;
int prevBrightness = 0;

char input[512];
char output[512];

void idle() {
	if (glutGetWindow() != main_window)
		glutSetWindow(main_window);
	glutPostRedisplay();
}

pixel** slice(pixel **m, int sx, int sy, int width, int height) {
	pixel **s = new pixel*[height];

	for (int i = 0; i < height; ++i) {
		s[i] = new pixel[width];
	}

	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			s[x][y].value = m[x + sx][y + sy].value;
		}
	}

	return s;
}

void computeUIWidth() {
	int aux, vw;
	GLUI_Master.get_viewport_area(&aux, &aux, &vw, &aux);
	UI_width = img->GetWidth() - vw;
}

void reshape(int x, int y) {
	glutPostRedisplay();
}

static void display(void) {
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glRasterPos2i(-1, -1); // Fix raster position
	imgMod->ViewImage();
	GLUI_Master.auto_set_viewport();

	glutSwapBuffers();
}

void applyFFT() {
	if (img) {
		double *data = new double[imgMod->GetHeight() * imgMod->GetWidth()];

		for (int y = 0; y < imgMod->GetHeight(); y++) {
			for (int x = 0; x < imgMod->GetWidth(); x++) {
				data[x + y * imgMod->GetWidth()] = imgMod->GetGrayValue(x, y);
			}
		}

		GFFT<17, double> fourier;

		fourier.fft(data);

		for (int y = 0; y < imgMod->GetHeight(); y++) {
			for (int x = 0; x < imgMod->GetWidth(); x++) {
				imgMod->SetGrayValue(x, y, data[x + y * imgMod->GetWidth()]);
			}
		}

		free(data);
	}
}

void applyIFFT() {
	if (img) {

	}
}

void shiftFFT() {
  	PixelLab tmp, q0, q1, q2, q3;

  	pixel **m;

	// first crop the image, if it has an odd number of rows or columns

  	/*PixelLab oi;
	oi.CreateImage(imgMod->GetWidth() & -2, imgMod->GetHeight() & -2, 1);*/
	imgMod->AllocatePixelMatrix(&m, imgMod->GetWidth(), imgMod->GetHeight());
	imgMod->GetDataAsMatrix(m);
	//image = image(Rect(0, 0, image.GetWidth() & -2, image.GetHeight() & -2));

	int cx = imgMod->GetWidth() / 2;
	int cy = imgMod->GetHeight() / 2;

	// rearrange the quadrants of Fourier image
	// so that the origin is at the image center

	// q0 = image(Rect(0, 0, cx, cy));
	// q1 = image(Rect(cx, 0, cx, cy));
	// q2 = image(Rect(0, cy, cx, cy));
	// q3 = image(Rect(cx, cy, cx, cy));

	q0.CreateImage(cx, cy);
	q1.CreateImage(cx, cy);
	q2.CreateImage(cx, cy);
	q3.CreateImage(cx, cy);
	tmp.CreateImage(cx, cy);


	q0.SetDataAsMatrix(slice(m, 0, 0, cx, cy));
	q1.SetDataAsMatrix(slice(m, cx, 0, cx, cy));
	q2.SetDataAsMatrix(slice(m, 0, cy, cx, cy));
	q3.SetDataAsMatrix(slice(m, cx, cy, cx, cy));

	/*q0.Copy(&tmp);
	q3.Copy(&q0);
	tmp.Copy(&q3);

	q1.Copy(&tmp);
	q2.Copy(&q1);
	tmp.Copy(&q2);*/

	for (int y = 0; y < cy; y++) {
		for (int x = 0; x < cx; x++) {
			imgMod->SetGrayValue(x, y, q3.GetGrayValue(x, y));
			imgMod->SetGrayValue(x + cx, y, q2.GetGrayValue(x, y));
			imgMod->SetGrayValue(x, y + cy, q1.GetGrayValue(x, y));
			imgMod->SetGrayValue(x + cx, y + cy, q0.GetGrayValue(x, y));
		}
	}

	// imgMod->SetDataAsMatrix(m);
	imgMod->DeallocatePixelMatrix(&m, imgMod->GetHeight(), imgMod->GetWidth());
}

void control(int value) {
	switch (value) {
	case 1:
		img->Read(input);
		imgMod->Copy(img);
		glutReshapeWindow(imgMod->GetWidth() + UI_width, imgMod->GetHeight());
		break;
	case 2:
		imgMod->Save(output);
		break;
	case 3:
		if (option == 0) {
			shiftFFT();
			//applyFFT();
		} else if (option == 1) {
			applyIFFT();
		}
		break;
	case 4:
		imgMod->Copy(img);
	}
}

static void key(unsigned char key, int x, int y) {
	switch (key) {
	case 27:
	case 'q':
		exit(0);
		break;
	}
	glutPostRedisplay();
}

int main(int argc, char *argv[]) {
	glutInit(&argc, argv);
	strcpy(input, "figs/lenaGray.png");
	strcpy(output, "figs/output.png");

	img = new PixelLab(input);
	imgMod = new PixelLab();
	imgMod->Copy(img);

	glutInitWindowSize(img->GetWidth(), img->GetHeight());
	glutInitWindowPosition(100, 100);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

	main_window = glutCreateWindow("Image Filtering in the Spatial Domain");

	glutKeyboardFunc(key);
	glutIdleFunc(idle);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);

	// GLUI
	// Use the first line above to get the interface in the same window as the graphics
	// Use the second line (comment the first) to get the interface in a separeted window
	GLUI *glui = GLUI_Master.create_glui_subwindow(main_window,
	GLUI_SUBWINDOW_RIGHT);
	//GLUI *glui = GLUI_Master.create_glui( "" );

	GLUI_Master.set_glutReshapeFunc(reshape);

	GLUI_Panel *io_panel = glui->add_panel((char *) "Input/Output");
	GLUI_EditText *edit1, *edit2;
	edit1 = glui->add_edittext_to_panel(io_panel, (char *) "",
	GLUI_EDITTEXT_TEXT, input);
	edit2 = glui->add_edittext_to_panel(io_panel, (char *) "",
	GLUI_EDITTEXT_TEXT, output);
	edit1->set_w(185);
	edit2->set_w(185);
	glui->add_column_to_panel(io_panel, false);
	GLUI_Button *b1 = glui->add_button_to_panel(io_panel, (char*) "Load", 1,
			control);
	GLUI_Button *b2 = glui->add_button_to_panel(io_panel, (char*) "Save", 2,
			control);
	b1->set_w(50);
	b2->set_w(50);

	GLUI_Panel *channel_panel = glui->add_panel((char *) "Transform");
	//glui->add_statictext_to_panel( channel_panel, (char *) "Low-pass");
	radio = glui->add_radiogroup_to_panel(channel_panel, &option, 3, (GLUI_Update_CB) NULL);
	glui->add_radiobutton_to_group(radio, (char *) "Fast Fourier Transform");
	glui->add_radiobutton_to_group(radio, (char *) "Inverse Fast Fourier Transform");

	glui->add_button((char *) "Apply", 3, control);
	glui->add_button((char *) "Reset", 4, control);
	glui->add_button((char *) "Quit", 0, (GLUI_Update_CB) exit);

	glui->set_main_gfx_window(main_window);
	GLUI_Master.set_glutIdleFunc(idle);
	computeUIWidth(); // Compute the size of the user interface
	glutReshapeWindow(imgMod->GetWidth() + UI_width, imgMod->GetHeight());

	glutMainLoop();

	return 1;
}
