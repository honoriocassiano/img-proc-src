#include "pixelLab.h"
#include <fstream>
#include <iostream>
#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include "png.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BMP_TYPE 1
#define PNG_TYPE 2

using namespace std;

//////////////////////////////////////////////////////////////////
PixelLab::PixelLab()
{
   this->grayScale = false;
   this->data = NULL;
   this->width = this->height = 0;
}

//////////////////////////////////////////////////////////////////
PixelLab::PixelLab(const char *fname)
{
   this->grayScale = false;
   this->data = NULL;
   this->width = this->height = 0;
   this->Read(fname);
}

//////////////////////////////////////////////////////////////////
PixelLab::~PixelLab()
{
	if( data )
		delete[] data;
}

//////////////////////////////////////////////////////////////////
void PixelLab::CreateImage(int w, int h, uByte *d, int channels )
{
   this->width    = w;
   this->height   = h;
   this->data     = d;

   this->channels = channels;
   if(this->channels == 1) this->grayScale = true;
}

//////////////////////////////////////////////////////////////////
void PixelLab::CreateImage(int w, int h, int channels )
{
   this->width    = w;
   this->height   = h;
   unsigned size(w * h * 3); // size of the image in bytes (3 is to RGB component).
   this->data = new uByte[size];

   this->channels = channels;
   if(this->channels == 1) this->grayScale = true;
}

//////////////////////////////////////////////////////////////////
int PixelLab::CheckFileType(char *name)
{
   int s = strlen(name);
   if( toupper (name[ s-3 ]) == 'B' &&
       toupper (name[ s-2 ]) == 'M' &&
       toupper (name[ s-1 ]) == 'P')
       return 1;
   if( toupper (name[ s-3 ]) == 'P' &&
       toupper (name[ s-2 ]) == 'N' &&
       toupper (name[ s-1 ]) == 'G')
       return 2;
   return 0;
}

//////////////////////////////////////////////////////////////////
void PixelLab::Read(const char *fname)
{
   // Erase previous data
   if( this->data ) delete[] data;

   int type = this->CheckFileType( (char *) fname);
   switch(type)
   {
      case BMP_TYPE: this->ReadBMP(fname); break;
      case PNG_TYPE: this->ReadPNG(fname); break;
      default:
      {
         cout << "Error: Only PNG or BMP files are allowed." << endl;
         exit(1);
      }
   }
}

//////////////////////////////////////////////////////////////////
void PixelLab::Save(const char *fname)
{
   int type = this->CheckFileType( (char *) fname);
   switch(type)
   {
      case BMP_TYPE: this->SaveBMP(fname); break;
      case PNG_TYPE: this->SavePNG(fname); break;
      default:
      {
         cout << "Error: Only PNG or BMP files are allowed." << endl;
         exit(1);
      }
   }
}

//////////////////////////////////////////////////////////////////
void PixelLab::PrintImageInfo()
{
   cout << "------------------------------------------------" << endl;
   cout << "Image Information" << endl;
   cout << "  Image Name:  " << this->imageName << endl;
   cout << "  Image Type:  " << this->imageType << endl;
   cout << "  Dimension:   " << this->width << " x " << this->height << " pixels" << endl;
   cout << "  Channels:    " << this->channels << endl;
   if( this->channels == 4 )
      cout << "\n Warning!!! RBGA images may not appear correctly\n" << endl;
   cout << "------------------------------------------------" << endl;
}

//////////////////////////////////////////////////////////////////
void PixelLab::ReadPNG(const char *file_name)
{
	char header[8];	// 8 is the maximum size that can be checked
   int k;
   //int number_of_passes;
   png_structp png_ptr;
   png_infop info_ptr;
   //png_byte color_type;
//   png_byte bit_depth;
   png_bytep * row_pointers;

	/* open file and test for it being a png */
	FILE *fp = fopen(file_name, "rb");
	if (!fp)
	{
		cout << "Error [PixelLab::ReadPNG]:\nFile \"" << file_name << "\" could not be opened for reading." << endl;
		exit(1);
	}
	fread(header, 1, 8, fp);

	if (png_sig_cmp((png_byte*)header, 0, 8))
   {
		cout << "Error in [PixelLab::ReadPNG]:\nFile \"" << file_name << "\" is not recognized as a PNG file." << endl;
		exit(1);
   }
	// initialize stuff
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (!png_ptr)
   {
		cout << "Error in [PixelLab::ReadPNG]: png_create_read_struct failed." << endl;
		exit(1);
   }

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
   {
		cout << "Error in [PixelLab::ReadPNG]: png_create_read_struct failed." << endl;
		exit(1);
   }

	if (setjmp(png_jmpbuf(png_ptr)))
   {
		cout << "Error in [PixelLab::ReadPNG]: Error during init_io." << endl;
		exit(1);
   }

	png_init_io(png_ptr, fp);
	png_set_sig_bytes(png_ptr, 8);

	png_read_info(png_ptr, info_ptr);

	//number_of_passes = png_set_interlace_handling(png_ptr);
	png_read_update_info(png_ptr, info_ptr);

	// read file
	if (setjmp(png_jmpbuf(png_ptr)))
   {
		cout << "Error in [PixelLab::ReadPNG]: Error during read_image." << endl;
		exit(1);
   }

	//color_type = info_ptr->color_type;
//	bit_depth = info_ptr->bit_depth;
	this->width = info_ptr->width;
	this->height = info_ptr->height;
   this->channels = info_ptr->channels;
   strcpy(this->imageName, file_name);
   strcpy(this->imageType, "PNG");

   // If the code get here, print image info
   this->PrintImageInfo();

	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
	for (int y=0; y<height; y++)
		row_pointers[y] = (png_byte*) malloc(info_ptr->rowbytes);

	png_read_image(png_ptr, row_pointers);

	int size(this->width * this->height * 3);
	this->data = new uByte[size];

   // copy data to pixel lab struct
   k = 0;
	for (int y=height-1; y>=0; y--)
	{
		png_byte* row = row_pointers[y];
		for (int x=0; x<width; x++)
		{
			png_byte* ptr = &(row[x * this->channels ]);
         if(this->channels == 1)
            data[k  ] = data[k+1] = data[k+2] = ptr[0];
			if(this->channels == 3)
         {
         data[k  ] = ptr[0];
         data[k+1] = ptr[1];
         data[k+2] = ptr[2];
         }
         k+=3;
		}
	}
   fclose(fp);
}

//////////////////////////////////////////////////////////////////
// Save a PNG File
// * In this implementation, bit depth is always 8
// * Color type can be RGB or Grayscale.
void PixelLab::SavePNG(const char *file_name)
{
   png_bytep * row_pointers;

   cout << "Saving File \"" << file_name << "\"";

	// create file
	FILE *fp = fopen(file_name, "wb");
	if (!fp)
	{
		cout << "\nError in [PixelLab::SavePNG]: File \"" << file_name << "\" could not be opened for writing" << endl;
		exit(1);
	}

	// initialize stuff
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (!png_ptr)
	{
		cout << "\nError in [PixelLab::SavePNG]: png_create_write_struct failed" << endl;
		exit(1);
	}

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
	{
		cout << "\nError in [PixelLab::SavePNG]: png_create_info_struct failed" << endl;
		exit(1);
	}

	if (setjmp(png_jmpbuf(png_ptr)))
	{
		cout << "\nError in [PixelLab::SavePNG]: Error during init_io" << endl;
		exit(1);
	}

	png_init_io(png_ptr, fp);


	// write header
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		cout << "\nError in [PixelLab::SavePNG]: Error during writing header" << endl;
		exit(1);
	}

   // Setting color type
   png_byte color_type = PNG_COLOR_TYPE_RGB;
   if( this->channels == 1) color_type = PNG_COLOR_TYPE_GRAY;

   // Setting bit_depth to 8
   png_byte bit_depth = 8; // Default value

	png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, color_type,
	             PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);

	// write bytes
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		cout << "\nError in [PixelLab::SavePNG]: Error during writing bytes" << endl;
		exit(1);
	}

	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * this->height);
	for (int y=0; y<this->height; y++)
      row_pointers[y] = (png_byte*) malloc( this->width * this->channels );

   // copy data from pixelLab struct to row_pointers
   int k = 0;
   int aux;
	for (int y=height-1; y>=0; y--)
	{
	   aux = 0;
		png_byte* row = row_pointers[y];
		for (int x=0; x < width; x++)
		{
         if(this->channels == 1)
         {
            row[aux++] = data[k];
            k+=3; // Even grayscale image are stored as 3 channel image
         }
			if(this->channels == 3)
         {
            row[aux++] = data[k++];
            row[aux++] = data[k++];
            row[aux++] = data[k++];
         }
		}
	}

   // Write data
	png_write_image(png_ptr, row_pointers);

	// end write
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		cout << "\nError in [PixelLab::SavePNG]: Error during end of write" << endl;
		exit(1);
	}

	png_write_end(png_ptr, NULL);

   // cleanup heap allocation
	for (int y=0; y<height; y++)
		free(row_pointers[y]);
	free(row_pointers);

   fclose(fp);
   cout << "      [DONE]"<< endl;
}

//////////////////////////////////////////////////////////////////
void PixelLab::ReadBMP(const char *fname)
{
	unsigned short planes;	// number of planes in image (must be 1)
	unsigned short bpp;			// number of bits per pixel (must be 24)

   strcpy(this->imageName, fname);
   strcpy(this->imageType, "BMP");

	ifstream fin(fname, ios::in | ios::binary);
	if( !fin )
	{
		cerr << "File not found " << fname << '\n';
		exit(1);
	}

	fin.seekg(18, ios::cur);

	fin.read((char *)&width, sizeof(int));
	fin.read((char *)&height, sizeof(int));

	fin.read((char *)&planes, sizeof(unsigned short));
	if( planes != 1 )
	{
		cout << "Planes from " << fname << " is not 1: " << planes << "\n";
		exit(1);
	}

	fin.read((char *)&bpp, sizeof(unsigned short));
	if( bpp != 24 )
	{
		cout << "Error: Bpp from " << fname << " is not 24: " << bpp << "\n";
		cout << "If you want to use grayscale images, use PNG format." << endl;
		exit(1);
	}

	fin.seekg(24, ios::cur);
	this->channels = 3;

	int size(width * height * 3); // size of the image in bytes (3 is to RGB component).
	data = new uByte[size];
	fin.read((char *)data, size);

	uByte tmp;	// temporary color storage for bgr-rgb conversion.

	for( int i(0); i <  size; i += 3 )
	{
		tmp = data[i];
		data[i] = data[i+2];
		data[i+2] = tmp;
	}

	this->PrintImageInfo();
}

//////////////////////////////////////////////////////////////////
void PixelLab::SaveBMP(const char *fname)
{
#ifdef WIN32
   cout << endl << "Saving file \"" << fname << "\"";
	using namespace std;

   BITMAPINFOHEADER bmpInfoHeader = {0};
   bmpInfoHeader.biSize = sizeof(BITMAPINFOHEADER);            // Set the size
   bmpInfoHeader.biBitCount = 24;                              // Bit count
   bmpInfoHeader.biClrImportant = 0;                           // Use all colors
   bmpInfoHeader.biClrUsed = 0;                                // Use as many colors according to bits per pixel
   bmpInfoHeader.biCompression = BI_RGB;                       // Store as un Compressed
   bmpInfoHeader.biHeight = this->height;                      // Set the height in pixels
   bmpInfoHeader.biWidth = this->width;                        // Width of the Image in pixels
   bmpInfoHeader.biPlanes = 1;                                 // Default number of planes
   bmpInfoHeader.biSizeImage = this->width * this->height * 3; // Calculate the image size in bytes
   BITMAPFILEHEADER bfh = {0};

   // This value should be values of BM letters i.e 0×4D42
   // 0×4D = M 0×42 = B storing in reverse order to match with endian
   bfh.bfType=0x4D42;

   // Offset to the RGBQUAD
   bfh.bfOffBits = sizeof(BITMAPINFOHEADER) + sizeof(BITMAPFILEHEADER);
   // Total size of image including size of headers
   bfh.bfSize = bfh.bfOffBits + bmpInfoHeader.biSizeImage;
   // Create the file in disk to write
   void* hFile = CreateFile( fname,GENERIC_WRITE, 0,NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL,NULL);

   if( !hFile ) return; // return if error opening file

   DWORD dwWritten = 0;
   // Write the File header
   WriteFile( hFile, &bfh, sizeof(bfh), &dwWritten , NULL );
   // Write the bitmap info header
   WriteFile( hFile, &bmpInfoHeader, sizeof(bmpInfoHeader), &dwWritten, NULL );

   // Convert from RGB to BGR
   int size(width * height * 3);
   uByte *tempData = new uByte[size];
	for( int i(0); i <  size; i += 3 )
	{
		tempData[i]   = data[i+2];
      tempData[i+1] = data[i+1];
      tempData[i+2] = data[i];
	}

   // Write the RGB Data
   WriteFile( hFile, tempData, bmpInfoHeader.biSizeImage, &dwWritten, NULL );

   // Delete temporary array
   delete[] tempData;

   // Close the file handle
   CloseHandle( hFile );
   cout << "...Done!" << endl;
#else
   cout << endl << "Warning: \"Save\" method in BMP file type is available only for windows (for now)." << endl;
   cout << "          Change the extension to PNG to save your file." << endl;
#endif
}

//////////////////////////////////////////////////////////////////
void PixelLab::ViewImage()
{
   glDrawPixels(this->GetWidth(), this->GetHeight(), GL_RGB, GL_UNSIGNED_BYTE, (uByte*) this->GetData());
}

//////////////////////////////////////////////////////////////////
uByte PixelLab::pixel_elem(int x, int y, int elem)
{
	int pos = (y*width+x) * 3 + elem;
	return data[pos];
}

//////////////////////////////////////////////////////////////////
uByte *PixelLab::pixel_pos(int x, int y)
{
	int pos = (y * width + x) * 3;
	return &data[pos];
}

//////////////////////////////////////////////////////////////////
int PixelLab::GetWidth()
{
   return this->width;
}

//////////////////////////////////////////////////////////////////
void PixelLab::SetWidth(int w)
{
   this->width = w;
}

//////////////////////////////////////////////////////////////////
int PixelLab::GetHeight()
{
   return this->height;
}

//////////////////////////////////////////////////////////////////
void PixelLab::SetHeight(int h)
{
   this->height = h;
}

//////////////////////////////////////////////////////////////////
uByte* PixelLab::GetData()
{
   return this->data;
}

//////////////////////////////////////////////////////////////////
void PixelLab::SetData(uByte *d)
{
   this->data = d;
}

//////////////////////////////////////////////////////////////////
void PixelLab::ConvertToGrayScale()
{
   this->grayScale = true;
   double value;
	int k = 0;
   for(int lin = 0; lin < this->GetHeight(); lin++)
      for(int col = 0; col < this->GetWidth(); col++)
      {
         value =  0.3 * (double) data[k] + 0.59 * (double) data[k+1] + 0.11 * (double) data[k+2];
         data[k] = data[k+1] = data[k+2] = (uByte) value;
         k+=3;
      }
}

//////////////////////////////////////////////////////////////////
void PixelLab::AddValueToChannels(int value)
{
	int size = this->height * this->width * 3;

   for(int i = 0; i < size; i++)
   {
      int intval = data[i]; // because 'data' is unsigned char
      if     ( intval + value > 255) data[i] = 255;
      else if( intval + value < 0)   data[i] = 0;   // Value can be negative
      else                           data[i]+=value;
   }
}

//////////////////////////////////////////////////////////////////
void PixelLab::SetRGB(int x, int y, uByte R, uByte G, uByte B)
{
	int pos = (y*width+x) * 3;
	data[pos] = R;
	data[pos+1] = G;
	data[pos+2] = B;
}

//////////////////////////////////////////////////////////////////
void PixelLab::GetRGB(int x, int y, uByte &R, uByte &G, uByte &B)
{
	int pos = (y*width+x) * 3;
	R = data[pos];
	G = data[pos+1];
	B = data[pos+2];
}

//////////////////////////////////////////////////////////////////
uByte PixelLab::GetR(int x, int y)
{
	int pos = (y*width+x) * 3;
	return data[pos];
}

//////////////////////////////////////////////////////////////////
uByte PixelLab::GetG(int x, int y)
{
	int pos = (y*width+x) * 3;
	return data[pos+1];
}

//////////////////////////////////////////////////////////////////
uByte PixelLab::GetB(int x, int y)
{
	int pos = (y*width+x) * 3;
	return data[pos+2];
}

//////////////////////////////////////////////////////////////////
void PixelLab::SetR(int x, int y, uByte R)
{
	int pos = (y*width+x) * 3;
	data[pos] = R;
}

//////////////////////////////////////////////////////////////////
void PixelLab::SetG(int x, int y, uByte G)
{
	int pos = (y*width+x) * 3;
	data[pos+1] = G;
}

//////////////////////////////////////////////////////////////////
void PixelLab::SetB(int x, int y, uByte B)
{
	int pos = (y*width+x) * 3;
	data[pos+2] = B;
}

//////////////////////////////////////////////////////////////////
void PixelLab::SetGrayValue(int x, int y, uByte value)
{
	int pos = (y*width+x) * 3;
	data[pos]   = value;
	data[pos+1] = value;
	data[pos+2] = value;
}

//////////////////////////////////////////////////////////////////
uByte PixelLab::GetGrayValue(int x, int y)
{
	int pos = (y*width+x) * 3;
	return (data[pos] + data[pos+1] + data[pos+2])/3;
}

//////////////////////////////////////////////////////////////////
void PixelLab::ConvertPixelsToData(pixel **m)
{
   int lin,col,k = 0;
   for(lin = 0; lin < this->GetHeight(); lin++)
      for(col = 0; col < this->GetWidth(); col++)
      {
         // If grayscale, set all components to value
         if(this->grayScale)
            m[lin][col].R = m[lin][col].G = m[lin][col].B = m[lin][col].value;
         data[k]   = m[lin][col].R;
         data[k+1] = m[lin][col].G;
         data[k+2] = m[lin][col].B;
         k+=3;
      }
}

//////////////////////////////////////////////////////////////////
void PixelLab::ConvertDataToPixels(pixel **m)
{
   int lin, col;
   uByte *dataAux = this->data;
   for(lin = 0; lin < this->GetHeight(); lin++)
      for(col = 0; col < this->GetWidth(); col++)
      {
         if(this->grayScale)
            m[lin][col].value = *dataAux; // All components are the same in grayscale mode
         m[lin][col].R = *dataAux++;
         m[lin][col].G = *dataAux++;
         m[lin][col].B = *dataAux++;
      }
}

//////////////////////////////////////////////////////////////////
void PixelLab::SetDataAsMatrix(pixel **m)
{
   int lin,col,k = 0;
   for(lin = 0; lin < this->GetHeight(); lin++)
      for(col = 0; col < this->GetWidth(); col++)
      {
         // If grayscale, set all components to value
         if(this->channels == 1)
            m[lin][col].R = m[lin][col].G = m[lin][col].B = m[lin][col].value;
         data[k]   = m[lin][col].R;
         data[k+1] = m[lin][col].G;
         data[k+2] = m[lin][col].B;
         k+=3;
      }
}

//////////////////////////////////////////////////////////////////
void PixelLab::GetDataAsMatrix(pixel **m)
{
   int lin, col;
   uByte *dataAux = this->data;
   for(lin = 0; lin < this->GetHeight(); lin++)
      for(col = 0; col < this->GetWidth(); col++)
      {
         if(this->channels == 1)
            m[lin][col].value = *dataAux; // All components are the same in grayscale mode
         m[lin][col].R = *dataAux++;
         m[lin][col].G = *dataAux++;
         m[lin][col].B = *dataAux++;
      }

}

//////////////////////////////////////////////////////////////////
void PixelLab::Copy( PixelLab *origem)
{
   // Erase previous data
   if( this->data ) delete[] data;

   // Copy parameters
   this->width    = origem->GetWidth();
   this->height   = origem->GetHeight();
   this->channels = origem->channels;

   if(this->channels == 1) this->grayScale = true;

   unsigned size(this->width * this->height * 3); // size of the image in bytes (3 is to RGB component).
   this->data = new uByte[size];

   // Copy Data
   for( int i = 0; i < (int) size; i++)
      this->data[i] = origem->data[i];
}

//////////////////////////////////////////////////////////////////
void PixelLab::SetNumberOfChannels(int n)
{
   this->channels = n;
   if(this->channels == 1) this->grayScale = true;
}

//////////////////////////////////////////////////////////////////
int PixelLab::GetNumberOfChannels()
{
   return this->channels;
}

//////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS
//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
void PixelLab::AllocatePixelMatrix(pixel ***m, int lin, int col)
{
   pixel **maux; // Auxiliary matrix

   // Dynamic allocation of lines
   maux = new pixel* [lin];
   // Dynamic allocation of columns
   for(int i = 0; i < lin; i++)
      maux[i] = new pixel[col];

   *m = maux; // Update de 'm' pointer
}

//////////////////////////////////////////////////////////////////
void PixelLab::DeallocatePixelMatrix(pixel ***m, int lin, int col)
{
   pixel **maux; // Auxiliary matrix
   maux = *m;
   // Dynamic deallocation of lines
   for(int i = 0; i < lin; i++)
      delete maux[i];
   // Dynamic deallocation of columns
   delete maux;
   maux = NULL;
}
