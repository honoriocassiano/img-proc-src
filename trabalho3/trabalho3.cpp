#include "pixelLab.h"
#ifdef WIN32
#include <windows.h>
#endif
#include <stdlib.h>
#include <math.h>
#include "GL/glui.h"
#include "GL/glut.h"

#include "filter_matrix.h"

// Image Objects
PixelLab *img    = NULL;
PixelLab *imgMod = NULL;

GLUI_Panel *size_panel;
GLUI_RadioGroup *radio, *radio2;
GLUI_RadioButton *buttons[3];
int main_window;
int UI_width  = 0;
int option = 0;
int option2 = 0;
int brightness = 0;
int prevBrightness = 0;

char input[512];
char output[512];

int getLowPassFilter(int filter, int size, int x, int y) {
	return low[size][filter][x][y];
}

int getHighPassFilter(int filter, int opt, int x, int y) {
	return high[opt][filter][x][y];
}

void idle() {
   if ( glutGetWindow() != main_window)
      glutSetWindow(main_window);
   glutPostRedisplay();
}

void computeUIWidth() {
   int aux, vw;
   GLUI_Master.get_viewport_area(&aux, &aux, &vw, &aux);
   UI_width = img->GetWidth()-vw;
}

void reshape( int x, int y ) {
  glutPostRedisplay();
}

static void display(void) {
   glClearColor (0.0, 0.0, 0.0, 0.0);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glRasterPos2i(-1, -1); // Fix raster position
   imgMod->ViewImage();
   GLUI_Master.auto_set_viewport();

   glutSwapBuffers();
}

void modifyImageLowPass(int filter, int size) {
   if (img) {
	   imgMod->Copy(img);

	   if (size < 0 || size > 2) return;

      for (int y = 0; y < img->GetHeight(); y++) {
         for (int x = 0; x < img->GetWidth(); x++) {
			int sum[3] = {0, 0, 0};
			int tempX, tempY;
			int amount = 0;

			for (int i = 0; i < sizeMap[size]; i++) {
				for (int j = 0; j < sizeMap[size]; j++) {
					tempX = x + j - (sizeMap[size] / 2);
					tempY = y + i - (sizeMap[size] / 2);

					if (tempX >= 0 && tempX < img->GetWidth()
							&& tempY >= 0 && tempY < img->GetHeight()) {
						sum[0] += img->GetR(tempX, tempY) * getLowPassFilter(filter, size, i, j);
						sum[1] += img->GetG(tempX, tempY) * getLowPassFilter(filter, size, i, j);
						sum[2] += img->GetB(tempX, tempY) * getLowPassFilter(filter, size, i, j);

						amount += getLowPassFilter(filter, size, i, j);
					}
				}
			}

			sum[0] /= amount;
			sum[1] /= amount;
			sum[2] /= amount;

			imgMod->SetRGB( x, y, sum[0], sum[1], sum[2]);
         }
      }
   }
}

void modifyImageHighPass(int filter, int opt) {
   if (img) {
	   imgMod->Copy(img);

	   if (opt == 0 || opt == 1) {
			for (int y = 0; y < img->GetHeight(); y++) {
				for (int x = 0; x < img->GetWidth(); x++) {
					int sum[3] = {0, 0, 0};
					int tempX, tempY;

					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 3; j++) {
							tempX = x + j - (3 / 2);
							tempY = y + i - (3 / 2);

							if (tempX >= 0 && tempX < img->GetWidth()
									&& tempY >= 0 && tempY < img->GetHeight()) {
								sum[0] += img->GetR(tempX, tempY) * getHighPassFilter(filter, opt, i, j);
								sum[1] += img->GetG(tempX, tempY) * getHighPassFilter(filter, opt, i, j);
								sum[2] += img->GetB(tempX, tempY) * getHighPassFilter(filter, opt, i, j);
							}
						}
					}

					if (sum[0] < 0) sum[0] = 0;
					if (sum[1] < 0) sum[1] = 0;
					if (sum[2] < 0) sum[2] = 0;

					if (sum[0] > 255) sum[0] = 255;
					if (sum[1] > 255) sum[1] = 255;
					if (sum[2] > 255) sum[2] = 255;

					imgMod->SetRGB( x, y, sum[0], sum[1], sum[2]);
				}
			}
	   } else if (opt == 2) {
		   for (int y = 0; y < img->GetHeight(); y++) {
				for (int x = 0; x < img->GetWidth(); x++) {
					int sum[2][3] = {{0, 0, 0}, {0, 0, 0}};
					int tempX, tempY;

					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 3; j++) {
							tempX = x + j - (3 / 2);
							tempY = y + i - (3 / 2);

							if (tempX >= 0 && tempX < img->GetWidth()
									&& tempY >= 0 && tempY < img->GetHeight()) {
								sum[0][0] += img->GetR(tempX, tempY)
										* getHighPassFilter(filter, HORIZ, i, j);
								sum[1][0] += img->GetR(tempX, tempY)
										* getHighPassFilter(filter, VERT, i, j);

								sum[0][1] += img->GetR(tempX, tempY)
										* getHighPassFilter(filter, HORIZ, i, j);
								sum[1][1] += img->GetR(tempX, tempY)
										* getHighPassFilter(filter, VERT, i, j);

								sum[0][2] += img->GetR(tempX, tempY)
										* getHighPassFilter(filter, HORIZ, i, j);
								sum[1][2] += img->GetR(tempX, tempY)
										* getHighPassFilter(filter, VERT, i, j);
							}
						}
					}

					if (sum[0][0] < 0) sum[0][0] = 0;
					if (sum[1][0] < 0) sum[1][0] = 0;
					if (sum[0][1] < 0) sum[0][1] = 0;
					if (sum[1][1] < 0) sum[1][1] = 0;
					if (sum[0][2] < 0) sum[0][2] = 0;
					if (sum[1][2] < 0) sum[1][2] = 0;

					if (sum[0][0] > 255) sum[0][0] = 255;
					if (sum[1][0] > 255) sum[1][0] = 255;
					if (sum[0][1] > 255) sum[0][1] = 255;
					if (sum[1][1] > 255) sum[1][1] = 255;
					if (sum[0][2] > 255) sum[0][2] = 255;
					if (sum[1][2] > 255) sum[1][2] = 255;

					int sumR = sum[0][0] + sum[1][0];
					int sumG = sum[0][1] + sum[1][1];
					int sumB = sum[0][2] + sum[1][2];

					if (sumR < 0) sumR = 0;
					if (sumG < 0) sumG = 0;
					if (sumB < 0) sumB = 0;

					if (sumR > 255) sumR = 255;
					if (sumG > 255) sumG = 255;
					if (sumB > 255) sumB = 255;

					imgMod->SetRGB( x, y, sumR, sumG, sumB);
				}
			}
	   }
   }
}

void control(int value)
{
   if(value == 1) {
      img->Read(input);
      imgMod->Copy(img);
      glutReshapeWindow(imgMod->GetWidth() + UI_width, imgMod->GetHeight());
   }
   if(value == 2) {
      imgMod->Save(output);
   }

   if (value == 3 || value == 4) { // Radio button
	   if (option == 0) {
		   size_panel->disable();

		   imgMod->Copy(img);
	   }
	   if (option == 1 || option == 2) {
		   size_panel->enable();

		   buttons[0]->set_name("3x3");
		   buttons[1]->set_name("5x5");
		   buttons[2]->set_name("7x7");

		   modifyImageLowPass(option - 1, option2);
	   }
	   if (option == 3 || option == 4) {
		   size_panel->enable();

		   buttons[0]->set_name("Horizontal");
		   buttons[1]->set_name("Vertical");
		   buttons[2]->set_name("Both");

		   modifyImageHighPass(option - 3, option2);
	   }
   }
}

static void key(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 27 :
        case 'q':
            exit(0);
        break;
    }
    glutPostRedisplay();
}

int main(int argc, char *argv[])
{
   glutInit(&argc, argv);
   strcpy(input, "figs/woman.png");
   strcpy(output, "figs/output.png");

   img    = new PixelLab(input);
   imgMod = new PixelLab();
   imgMod->Copy(img);

   glutInitWindowSize(img->GetWidth(),img->GetHeight());
   glutInitWindowPosition(100,100);
   glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

   main_window = glutCreateWindow("Image Filtering in the Spatial Domain");

   glutKeyboardFunc(key);
   glutIdleFunc(idle);
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);

   // GLUI
   // Use the first line above to get the interface in the same window as the graphics
   // Use the second line (comment the first) to get the interface in a separeted window
   GLUI *glui = GLUI_Master.create_glui_subwindow( main_window, GLUI_SUBWINDOW_RIGHT );
   //GLUI *glui = GLUI_Master.create_glui( "" );

   GLUI_Master.set_glutReshapeFunc(reshape);

   GLUI_Panel *io_panel = glui->add_panel( (char *) "Input/Output" );
   GLUI_EditText *edit1, *edit2;
   edit1 = glui->add_edittext_to_panel (io_panel, (char *) "", GLUI_EDITTEXT_TEXT, input);
   edit2 = glui->add_edittext_to_panel (io_panel, (char *) "", GLUI_EDITTEXT_TEXT, output);
   edit1->set_w(185);
   edit2->set_w(185);
   glui->add_column_to_panel( io_panel, false);
   GLUI_Button *b1 = glui->add_button_to_panel( io_panel, (char*) "Load", 1, control );
   GLUI_Button *b2 = glui->add_button_to_panel( io_panel, (char*) "Save", 2, control );
   b1->set_w(50);
   b2->set_w(50);

   GLUI_Panel *channel_panel = glui->add_panel((char *) "Filters" );
   //glui->add_statictext_to_panel( channel_panel, (char *) "Low-pass");
   radio = glui->add_radiogroup_to_panel(channel_panel, &option, 3, control);
   glui->add_radiobutton_to_group( radio, (char *) "No filter" );
   glui->add_radiobutton_to_group( radio, (char *) "Average" );
   glui->add_radiobutton_to_group( radio, (char *) "Gaussian" );
   glui->add_radiobutton_to_group( radio, (char *) "Sobel" );
   glui->add_radiobutton_to_group( radio, (char *) "Prewitt" );
   /*GLUI_Spinner* spinner;
   spinner = glui->add_spinner_to_panel( channel_panel, (char*)"Increase/Decrease brightness", GLUI_SPINNER_INT, &brightness, 4, control);*/

   size_panel = glui->add_panel((char *) "Option" );
   radio2 = glui->add_radiogroup_to_panel(size_panel, &option2, 4, control);
	buttons[0] = glui->add_radiobutton_to_group( radio2, (char *) "3x3" );
	buttons[1] = glui->add_radiobutton_to_group( radio2, (char *) "5x5" );
	buttons[2] = glui->add_radiobutton_to_group( radio2, (char *) "7x7" );

	size_panel->disable();

   glui->add_button((char *) "Quit", 0,(GLUI_Update_CB)exit );

   glui->set_main_gfx_window( main_window );
   GLUI_Master.set_glutIdleFunc( idle );
   computeUIWidth(); // Compute the size of the user interface
   glutReshapeWindow(imgMod->GetWidth() + UI_width, imgMod->GetHeight());

   glutMainLoop();

   return 1;
}
