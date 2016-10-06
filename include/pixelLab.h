/*
  Name:        PixelLab.h
  Version      0.6
  Author:      Rodrigo Luis de Souza da Silva
  Date:        04/07/2009
  Last Update: 28/06/2013
  Description: Simple class to load PNG/BMP images and modify its data.
               You can create an image manually and convert it to a
               pixel matrix and vice-versa.

  Updates:
  28/06/2013 - "ifndef" problem fixed
  21/11/2012 - Pixel lab is in a separate project and creates a dynamic library
*/


//////////////////////////////////////////////////////////////////
// Preprocessor directives
#define GLUT_NO_LIB_PRAGMA
#define GLUT_NO_WARNING_DISABLE
#define GLUT_DISABLE_ATEXIT_HACK

#ifndef PIXEL_LAB_H
#define PIXEL_LAB_H

//////////////////////////////////////////////////////////////////
// New types
//typedef char Byte; // conflicts with zconf.h
typedef unsigned char uByte;

//////////////////////////////////////////////////////////////////
// Classes
// Simple pixel Structure
struct pixel
{
   uByte value;      // Grayscale value
   uByte R, G, B, A; // Color RGBA values
};

// Pixel map class
class PixelLab
{
   public:
      // Constructors and Destructors
      PixelLab(const char *fname);
      PixelLab();
      ~PixelLab();

      // Read image from a bitmap file
      void Read(const char *fname);

      // Save image in a bitmap file
      void Save(const char *fname);

      // Create image from a 'uByte' array
      void CreateImage(int w, int h, uByte *d, int channels = 3);

      // Create an empty image with width x height pixels
      void CreateImage(int width, int height, int channels = 3);

      // Check if the input file name is allowed (BMP or PNG)
      int CheckFileType(char *name);

      // view stored image
      void ViewImage();

      // Set/Get Image as matrix.  The "set" method will allocate the matrix
      void SetDataAsMatrix(pixel **m);
      void GetDataAsMatrix(pixel **m);

      // Convertion methods (Deprecated)
      void ConvertPixelsToData(pixel **m);
      void ConvertDataToPixels(pixel **m);

      // Modify and get pixel information
      uByte  pixel_elem(int x, int y, int elem);
      uByte *pixel_pos(int x, int y);

      // Set/Get image's width and height
      int  GetWidth();
      void SetWidth( int);
      int  GetHeight();
      void SetHeight( int);

      // Set/Get image data
      uByte* GetData();
      void   SetData(uByte*);

      // Set/Get individual RGB values
      void  SetRGB(int x, int y, uByte R, uByte G, uByte B);
      void  GetRGB(int x, int y, uByte &R, uByte &G, uByte &B);
      uByte GetR(int x, int y);
      uByte GetG(int x, int y);
      uByte GetB(int x, int y);
      void  SetR(int x, int y, uByte);
      void  SetG(int x, int y, uByte);
      void  SetB(int x, int y, uByte);

      // Set/Get number of channels
      void SetNumberOfChannels(int);
      int  GetNumberOfChannels();

      // Set/Get Gray values
      void  SetGrayValue(int x, int y, uByte value);
      uByte GetGrayValue(int x, int y);

      // Convert Colored Images to Gray Scales
      void ConvertToGrayScale();

      // Copy PixelData Object
      void Copy( PixelLab *origem);

      // Print image information
      void PrintImageInfo();

      // Add value to all the channels. Values will be truncated if they are bigger than 255 or smaller than 0.
      void AddValueToChannels(int value);

      // Auxiliary functions
      void AllocatePixelMatrix  (pixel ***m, int lin, int col);
      void DeallocatePixelMatrix(pixel ***m, int lin, int col);

   private:
      // Internal saving methods
      void SaveBMP(const char *fname);
      void SavePNG(const char *fname);

      // Internal reading methods
      void ReadBMP(const char *fname);
      void ReadPNG(const char *fname);

   protected:
      uByte *data;
      int   width;
      int   height;
      bool  grayScale;
      int   channels;       // number of channels of the image
      char  imageName[512];
      char  imageType[4];
};

#endif // PIXEL_LAB_H
