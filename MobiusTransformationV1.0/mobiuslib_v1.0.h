/*****************************************************************************

			DECLARATIONS of MOBIUS TRANSFORMATION
        
			Copyright (C), Januari 2019, Jan de Nijs
		
*****************************************************************************/

#ifndef __MOBIUSLIB_V10__
#define __MOBIUSLIB_V10__

#include <CImg.h>
#ifdef Success
#undef Success
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <thread>
#include <vector>
#include <chrono>
#include <ctime>

using namespace std;
using Eigen::Matrix3f;
using Eigen::Vector3f;
using namespace cimg_library;
using namespace Eigen;

typedef Eigen::Matrix<unsigned char, 3, 4>  Color4; 		//Matrix  with 4 colomns RGB unsigned char.
typedef Eigen::Matrix<unsigned char, 3, 1>  Color1; 		//RGB unsigned char.
typedef Eigen::Matrix<int, 2, 4>  			Coordinates_Sq; //Matrix with 4 colums of X and Y coordinates.

#ifndef THREAD_NUMBER
#define THREAD_NUMBER	4
#endif

#define LENGTH_TRIANGLE_ARRAY 100

#ifndef MIN_SQUARE
#define MIN_SQUARE		4
#endif


struct MobiusSphere{            //Centre is with reference to the centre of the image (0,0,50) refers to a position on top of the image and in its centre
    float diameter, radius;
    float     centre[3];
    };
	
struct Parameters{
    float       	*ptrC[2];
    unsigned char   *ptrIm[2];
    MobiusSphere 	mobius;
    int         	rowThreadEntrances[THREAD_NUMBER][2];
    int         	imS[2];
    Matrix3f    	rotMatrix;
};
     
enum Coordinates {CARTH_3D, CARTH_2D_R};
enum Image {IMAGE, IMAGE_T};
enum Mode {BARE, TRIANGLES, QUADS};

MobiusSphere ReadMobiusSphere(void);
Vector3f ReadRotationAxisGraphical(CImg<unsigned char> image);
//ReadRotationalAxisGraphical provides a graphical interface to define the rotational axis. 
//The image is displayed and using the cursor and mouse a rotation axis is selected. 

string GetImage(void);

void ProjectImage2MobiusSphere(int const imSize[2], MobiusSphere &sphere, float *ptr3DX);
// Maps Carthesian 2D coordinates of the image to 3D coordinates on a sphere 
// *ptr3DX references to a predeclared CImg image type float, size imSize and 3 channels containing the (X,Y,Z)
// Euclidean coordinates of the projected pixels.
// input:	Width (imSize[0]) and Height (imSize[1]) of the image,
// input:	Struct defining the Mobius Sphere
// output:	Pixelcoordinates projected on the Mobius Sphere are written to the CImg image referenced by ptr3DX


void MobiusTransformation(unsigned char *ptrImR, float *ptrX, int const imSize[2], int const rows[2], Matrix3f rotationMatrix, MobiusSphere &sphere, float *ptrP2DX, unsigned char *ptrImTR);
//Rotates the Sphere and performs a reverse projection from the sphere to the 2D image plane.
//In addition, the pixels from the input image are copied to the appropriate positions in the output image.
// input:	ptrImR that references a CImg RGB image of Width imSize[0] and Height imSize[1],
// input:	ptrX that references to pixelcoordinates projected on the Mobius Sphere, 
//          stored in a CImg image of Width imSize[0] and Height imSize[1],
// input:	Width (imSize[0]) and Height (imSize[1]) of the image,
// input:	rows with rows[0] the start row and row[1] the end row refers to the rows to be transformed.
// input:	the rotation matrix that defines the rotation of the Mobius Sphere
// input:	Struct defining the Mobius Sphere
// output:	ptrP2DX that references a CImg image with 2 channels of type float and size imSize. The coordinates of 
//			pixels after the Mobius Transformation are written to this image
// output:	ptrImTR refers to new (output) image of format CImg RGB image of size (2 * imSize[0], 2 * imSize[1]).


void FillQuandrangulars_Recursive(unsigned char *ptrImTR, unsigned char *ptrImR, float *ptrTX, int imSize[2], int rows[2]);
// The output image comprises sparse areas with few pixels. Because of the continuity property of the transformation,
// the projected pixels form quadrangulars. To add the missing pixels, these quadrangulars are split in 2 or 4 subquadrangulars.
// This process is recursively applied. This functions works in conjunction withthe function "Segmentate"
// The RGB values of the new corners are calculated using the RGB values of the input corners.
// output:	ptrImTR references to new (output) image of format CImg RGB image of size (2 * imSize[0], 2 * imSize[1]).
// input:	ptrImR references to the input CImg RGB image of Width imSize[0] and Height imSize[1],
// input:	Width (imSize[0]) and Height (imSize[1]) of the image,
// input:	rows with rows[0] the start row and row[1] the end row refers to the rows to be transformed.


void SegmentateQuadrangular(unsigned char *ptrImTR, int imSize[2], Coordinates_Sq quadCoord, Color4 quadColors);
// This method segmentates a quadrangular and calculates the RGB color values of the new corners.
// This method works recursively.
// output:	ptrImTR references to new (output) image of format CImg RGB image of size (2 * imSize[0], 2 * imSize[1]).
//			Segmentate adds missing pixels in the output image of the Mobius-transformed image
// input:	Width (imSize[0]) and Height (imSize[1]) of the image,
// input:	quadCoord:	Eigen::Matrix<int, 2, 4> with the 4 corners of the quadrangular
// input:	quadColors:	Eigen::Matrix<unsigned char, 3, 4>  with the 4 RGB colors of the corners


void FillMissingPixels(unsigned char *ptrImTR, int const doubleSize[2], int rows[2]);
// FillMissingPixels adds missing pixels by assigning the color RGB values of the nearest neighbor with colors assigned.
// output:	ptrImTR references to new (output) image of format CImg RGB image of size (2 * imSize[0], 2 * imSize[1]).
//			FillMissingPixels adds missing pixels in the output image of the Mobius-transformed image
// input:	Width (imSize[0]) and Height (imSize[1]) of the image,
// input:	rows with rows[0] the start row and row[1] the end row refers to the rows to be transformed.



void FillQuandrangulars(unsigned char *ptrImTR, unsigned char *ptrImR, float *ptrImTX, int const imSize[2], int rows[2]);
// The output image comprises sparse areas with few pixels. Because of the continuity property of the transformation,
// the projected pixels form quadrangulars. To add the missing pixels, these quadrangulars are split in 2 triangles.
// Next, the method PaintTriangle is called to add the missing pixels in the triangles
// output:	ptrImTR references to new (output) image of format CImg RGB image of size (2 * imSize[0], 2 * imSize[1]).
// input:	ptrImR references to the input CImg RGB image of Width imSize[0] and Height imSize[1],
// input:	Width (imSize[0]) and Height (imSize[1]) of the image,
// input:	rows with rows[0] the start row and row[1] the end row refers to the rows to be transformed.


void FillTriangle(unsigned char *ptrImTR, int const imSize[2], int cTriangle[3][2], unsigned char cColours[3][3]);

// input:	ptrImR references to the input CImg RGB image of Width imSize[0] and Height imSize[1],



#endif