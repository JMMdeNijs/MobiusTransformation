/****************************************************************************************
Program performs a Mobius transformation on an image.

For information see http://www-users.math.umn.edu/~arnold/papers/moebius.pdf

Copyright (C), Jan de Nijs, December 2018

*****************************************************************************************/
#include "mobiuslib_v1.0.h"

//Compound functions for threads
void MobiusWrapper(Parameters *par, int i);
void MobiusTransformationQuadrangulars(Parameters *par, int i);
void MobiusTransformationTriangles(Parameters *par, int i);
	
/**************************************************************************************************************************

            Main Program

***************************************************************************************************************************/

int main(int argc,char **argv) {
    int rowThreadEntrances[THREAD_NUMBER][2];
	int rowRange[2], doubleRowRange[2];
	
	//Gather image
	string imDir(GetImage());
	
    //load image imName from directory imDir
    CImg<unsigned char> image(&imDir[0]);
	
    int const imageSize[2] ={image.width(), image.height()}, doubleSize[2] = {2*image.width(), 2*image.height()};
		
	rowThreadEntrances[0][0] = 0; rowThreadEntrances[0][1] = imageSize[1]/THREAD_NUMBER + 1 + imageSize[1]%THREAD_NUMBER;
	for(int i=1; i<THREAD_NUMBER; i++){
		rowThreadEntrances[i][0] = rowThreadEntrances[i-1][0] + imageSize[1]/THREAD_NUMBER;
		rowThreadEntrances[i][1] = rowThreadEntrances[i-1][1] + imageSize[1]/THREAD_NUMBER;
	}
	rowThreadEntrances[THREAD_NUMBER-1][1] -= 1; 
	
	rowRange[0] = 0; rowRange[1] = imageSize[1];
	doubleRowRange[0] = rowRange[0] *2; doubleRowRange[1] = rowRange[1] *2; 
	
	
    //declaration new variable for image after Mobius-Transformation
    CImg<unsigned char> imageTransformed(doubleSize[0],doubleSize[1],1,3,0), iFiltered(doubleSize[0],doubleSize[1],1,3,0);

    //declaration variables for the image coordinates after projection on the Mobius Sphere (carthesian3D)
	//and after rotation of the sphere and reverse projection to the 2D image plane (carthesian2DReverse)
    CImg<float> carthesian3D(imageSize[0],imageSize[1],1,3,0.0);
	CImg<float> carthesian2DReverse(imageSize[0],imageSize[1],1,2,0.0);
    
    //Input Mobius Sphere data
    MobiusSphere mobiusSphere = {100,50,0,0,50};//Definition Sphere for Mobius transformation Bol.Centre is with reference to the centre of the image (0,0,50) refers to a position on top of the image and in its centre
    
	//Declaration struct with thread-variables
	Parameters  par;
    par.ptrC[CARTH_3D] = carthesian3D.data();
    par.ptrC[CARTH_2D_R] = carthesian2DReverse.data();
    par.ptrIm[IMAGE]  = image.data();
    par.ptrIm[IMAGE_T] = imageTransformed.data();
    par.mobius = mobiusSphere;
    memcpy(par.rowThreadEntrances, rowThreadEntrances, sizeof(par.rowThreadEntrances));
    memcpy(par.imS, imageSize, sizeof(par.imS));
    
	//Timing for clock in image
	std::chrono::time_point<std::chrono::system_clock> clock_start, clock_end; 
	int clock_interval;
	char *textstring;
	string time_interval;
	
	ProjectImage2MobiusSphere(imageSize, mobiusSphere, carthesian3D.data());
	
    //input rotation axis of the Mobius sphere, see https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula#Matrix_notation
	Vector3f rotationAxis = ReadRotationAxisGraphical(image); 
    //Declaration matrices for rotation of the Mobius Sphere
	Matrix3f mat_K, mat_R, mat_I; mat_I << 1, 0, 0, 0, 1, 0, 0, 0, 1;

	//define a window to dispay transformed image
    CImgDisplay mobiusDisplay; 
	std::vector<std::thread> threadList;
	
	//Arrangements for rotation loop
	float angle(0), angleIncrement(1.0), angleIncrementPrevious;
	unsigned int fillNo(2);
	Mode  filling(TRIANGLES);
	unsigned char const yellow[3] = {0, 255, 255}, black[3] = {0, 0, 0};
	bool show_hide(true);
	float rotationAngle;
	
	while (angle <181){
        clock_start = std::chrono::system_clock::now(); 
    	rotationAngle=(angle*M_PI/90);
       
		//Calculation of the rotation matrix
        mat_K << 0,-rotationAxis(2), rotationAxis(1), rotationAxis(2), 0, -rotationAxis(0), -rotationAxis(1), rotationAxis(0), 0;
        mat_R = mat_I + mat_K*sin(rotationAngle) + mat_K*mat_K*(1-cos(rotationAngle));
        par.rotMatrix = mat_R;
		imageTransformed.fill(0);
		
		//Start processing threads 
		for (int i=0; i<THREAD_NUMBER; i++){
			if (filling == TRIANGLES) threadList.push_back(std::thread(MobiusTransformationTriangles, &par, i));
			if (filling == BARE) threadList.push_back(std::thread(MobiusWrapper, &par, i));
			if (filling == QUADS) threadList.push_back(std::thread(MobiusTransformationQuadrangulars, &par, i));
		}
        for (int i=0; i<THREAD_NUMBER; i++)
            (threadList.at(i)).join();
        threadList.erase (threadList.begin(),threadList.end());		
		
		//Post processing (add missing pixels)
		for (int ii=0; ii<fillNo; ii++) FillMissingPixels(imageTransformed.data(), doubleSize, doubleRowRange);
		
		//Display the image
		clock_interval = ceil((std::chrono::system_clock::now() - clock_start).count()/1e6); 
		time_interval = "Time per frame : " + to_string(clock_interval) + "  ms";
		imageTransformed.draw_text(10, 10, &time_interval[0], yellow, black, 1.0, 20);
		if (show_hide) {
			imageTransformed.draw_text(10, 2* imageSize[1] - 60, "Hide Show legend: (H), Proceed/Pause: P ", yellow, black, 1.0, 20);
			imageTransformed.draw_text(10, 2* imageSize[1] - 30, "Modes: BARE(B), QUADS (Q), TRIANGLES (T), Fill Missing Pixels: Arrow Up and Down", yellow, black, 1.0, 20);
		}
		mobiusDisplay.display(imageTransformed);
		mobiusDisplay.wait(100);
		
		//Process control: Pause / Continue of the stream	
		if (mobiusDisplay.is_keyP()) {
			if (angleIncrement > 0.01){
				angleIncrementPrevious = angleIncrement;
				angleIncrement = 0.;
			}
			else 
				angleIncrement = angleIncrementPrevious;
		}
		
		//Process control: Change Mode BARE, Traingles, QUADS,	
		if (mobiusDisplay.is_keyB()) filling = BARE;
		if (mobiusDisplay.is_keyT()) filling = TRIANGLES;
		if (mobiusDisplay.is_keyQ()) filling = QUADS;
		
		//Process control: Increase and decrease filling:
		if (mobiusDisplay.is_keyARROWUP()) fillNo++;
		if (mobiusDisplay.is_keyARROWDOWN() && fillNo >=1) fillNo--;
		
		//Process control: Hide (show) legend:
		if (mobiusDisplay.is_keyH()) show_hide = !show_hide;
		
		angle += angleIncrement;
    }
  while (!mobiusDisplay.is_closed() ) mobiusDisplay.wait();

  return 0;
}


//Wrappers for threads
void MobiusWrapper(Parameters *par, int i){
	MobiusTransformation(par->ptrIm[IMAGE], par->ptrC[CARTH_3D],  par->imS, par->rowThreadEntrances[i], par->rotMatrix, par->mobius, par->ptrC[CARTH_2D_R], par->ptrIm[IMAGE_T]);
	//void MobiusTransformationII(unsigned char *ptrImR, float *ptrX, int const imSize[2], int const rows[2], Matrix3f rotationMatrix, MobiusSphere &sphere, float *ptrPX2D, unsigned char *ptrImTR)
}


void MobiusTransformationQuadrangulars(Parameters *par, int i){
	MobiusTransformation(par->ptrIm[IMAGE], par->ptrC[CARTH_3D],  par->imS, par->rowThreadEntrances[i], par->rotMatrix, par->mobius, par->ptrC[CARTH_2D_R], par->ptrIm[IMAGE_T]);
	//Rough (incomplete) filling of all missing pixels
	FillQuandrangulars_Recursive(par->ptrIm[IMAGE_T], par->ptrIm[IMAGE], par->ptrC[CARTH_2D_R], par->imS, par->rowThreadEntrances[i]);
}

void MobiusTransformationTriangles(Parameters *par, int i){
	MobiusTransformation(par->ptrIm[IMAGE], par->ptrC[CARTH_3D],  par->imS, par->rowThreadEntrances[i], par->rotMatrix, par->mobius, par->ptrC[CARTH_2D_R], par->ptrIm[IMAGE_T]);
	//Rough (incomplete) filling of all missing pixels
	FillQuandrangulars(par->ptrIm[IMAGE_T], par->ptrIm[IMAGE], par->ptrC[CARTH_2D_R], par->imS, par->rowThreadEntrances[i]);
}

