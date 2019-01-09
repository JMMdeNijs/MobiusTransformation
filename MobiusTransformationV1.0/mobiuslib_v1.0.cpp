/**************************************************************************************************************************
      DEFINITIONS OF FUNCTIONS of MOBIUS TRANSFORMATION
        
		Copyright (C), Januari 2019, Jan de Nijs

***************************************************************************************************************************/


#include "mobiuslib_v1.0.h"


MobiusSphere ReadMobiusSphere(void){
        MobiusSphere sphere={100,50,0,0,50};
        char YN;
        cout << endl << endl << "The default Sphere parameters are:" << endl;
        cout << "\t Diameter:   100" << endl << "\t Centre position [0,0,50]" << endl;
        cout << "Change parameters (Y/N)?   : ";
        cin >> YN;
        if (YN=='Y' || YN == 'y'){
            cout << "Diameter : ";              cin >> sphere.diameter;
            cout << endl << "X-position :";     cin >> sphere.centre[0];
            cout << endl << "Y-position :";     cin >> sphere.centre[1];
            sphere.radius = sphere.centre[2] = sphere.diameter/2;
        }
        return sphere;
    }

	
Vector3f ReadRotationAxisGraphical(CImg<unsigned char> image){
        Vector3f res, res2;
		const unsigned char white[3] = {255, 255, 255}, yellow[3] = {0, 255, 255}, black[3] = {0, 0, 0};
		float opacity(0.5), rotX, rotY, rotZ, outOfCentre;
		int width(image.width()), height(image.height()), xPosition(-10), yPosition(-10);
		
		//Display image with lines to guide selection Mobius rotation axis
		for (int ii=-1; ii<2 ; ii++){
			image.draw_line(0, ii + height/2, width, ii+ height/2, white, opacity);
			image.draw_line(ii + width/2, 0, ii+ width/2, height, white, opacity);
			image.draw_circle(width/2, height/2, ii+ 0.9*(min(width/2, height/2)), white, opacity, 1);
			image.draw_circle(width/2, height/2, ii+ 0.1*(min(width/2, height/2)), white, opacity, 1);	
		}
		
		opacity = 1.0;
		image.draw_line(0, height*0.45, width, height*0.45, white, opacity);
		image.draw_line(0, height*0.55, width, height*0.55, white, opacity);
		image.draw_line(width*0.45, 0, width*0.45, height, white, opacity);
		image.draw_line(width*0.55, 0, width*0.55, height, white, opacity);
			
		image.draw_text(10, image.height() - 30, "Click on image to specify rotation axis", yellow, black, 1.0, 20);
		
		CImgDisplay image_display;
		image_display.display(image);	
		image_display.wait();
	
		//Make selections and calulate Mobius rotation axis	
		while (!image_display.is_closed() && xPosition < 0 && yPosition < 0)
			if (image_display.button()&1) {
				xPosition = image_display.mouse_x();
				yPosition = image_display.mouse_y();
			}	
		//calculate value rotX of rotation axis
		if (yPosition >= 0 && yPosition < 0.45 * height)
			rotX = yPosition / (0.45 * height);
		else if (yPosition > 0.55*height && yPosition <= height)
			rotX = (yPosition - height) / (0.45 * height);
		else rotX = 1;
		//calculate value rotY of rotation axis
		if (xPosition >= 0 && xPosition < 0.45 * width)
			rotY = xPosition / (0.45 * width);
		else if (xPosition > 0.55* width && xPosition <= width)
			rotY = (xPosition - width) / (0.45 * width);
		else rotY = 1;
		//calculate value rotZ of rotation axis
		outOfCentre = sqrt( pow(xPosition - width/2, 2) + pow(yPosition - height/2,2));
		if (outOfCentre >= 0 && outOfCentre < min(width, height)/20)
			if (xPosition < width/2)
				rotZ = 1;
			else rotZ = -1;
		else if (outOfCentre >= min(width, height)/20 && outOfCentre < 9*min(width, height)/20)
			if (xPosition < width/2)
				 rotZ = -(2*outOfCentre - min(width, height)) /min(width, height);
			else rotZ = 2*(outOfCentre - min(width, height)) /min(width, height);
		else
			rotZ = 0;
		//special cases
		if (xPosition > 0.45 * width && xPosition < 0.55 * width && yPosition > 0.45 * height && yPosition < 0.55 * height)
			if (xPosition < 0.5) res << 0, 0, 1;
			else res << 0, 0, -1;
		else if	(xPosition >= 0.45 * width  && xPosition < 0.5 * width)		res << 0,  1, 0;
		else if	(xPosition >= 0.5  * width  && xPosition < 0.55 * width)	res << 0, -1, 0;
		else if	(yPosition >= 0.45 * height && yPosition < 0.5 * height)	res <<  1, 0, 0;
		else if (yPosition >= 0.5  * height && yPosition < 0.55 * height)	res << -1, 0, 0;
		else res << rotX, rotY, rotZ;
		
        res.normalize();
        return res;
}	

	
string GetImage(void){
	std::string imDir, imName;

	bool found(false);
	while (!found) {
		std::cout << std::endl << "Give image name              : ";	std::cin >> imName;
		std::cout << std::endl << "Give path to image directory : "; 	std::cin >> imDir;
		for(unsigned int i=0; i < imDir.length(); i++)
			if (imDir[i] ==92) imDir[i]=47; //replaces '\'  by '/'
		imDir+= ("/" + imName);
		std::ifstream filecheck(imDir);
		found = (bool)filecheck;
		found? std::cout << "Image found \n" : std::cout << "Image not found \n\n";
	}
	
	return imDir;
}

void ProjectImage2MobiusSphere(int const imSize[2], MobiusSphere &sphere, float *ptr3DX){
    // Mapping Carthesian 2D coordinates of the image to 3D coordinates on a sphere 
    const int siz(imSize[0]*imSize[1]);
	float *ptr3DY(ptr3DX + siz), *ptr3DZ(ptr3DY + siz);
	float xx, yy;
    float halfWidth (float(imSize[0])/2), halfHeight (float(imSize[1])/2);
	const float corrX = 0.5 - float(imSize[0]%2)/2, corrY = 0.5 - float(imSize[1]%2)/2;
	float lengthSqXY, diamSQ(sphere.diameter*sphere.diameter), kPar;
	
    for (int i=0; i < imSize[1]; i++)
        for (int j=0; j < imSize[0]; j++){
			//calculate (x,y) cooridnate of the pixel
            xx = j - halfWidth - sphere.centre[0] + corrX;
            yy = halfHeight - i - sphere.centre[1] - corrY;
			//project (x,y) on the Mobius Sphere
			lengthSqXY = xx * xx + yy * yy;
            kPar = lengthSqXY/(lengthSqXY + diamSQ);
            *ptr3DZ++ = kPar * sphere.diameter - sphere.radius;
            *ptr3DY++ = (1.0 - kPar)* yy;
            *ptr3DX++ = (1.0 - kPar)* xx;
        }
}
	

void MobiusTransformation(unsigned char *ptrImR, float *ptrX, int const imSize[2], int const rows[2], Matrix3f rotationMatrix, MobiusSphere &sphere, float *ptrP2DX, unsigned char *ptrImTR){
	//Rotates the Sphere and performs a reverse projection from the sphere to the 2D image plane
	//In addition, the pixels from the input image are copied to the appropriate positions in the output image
	const int siz(imSize[0]*imSize[1]), addressRange((rows[1]-rows[0])*imSize[0]);
	const int imDoubleSize[2] = {imSize[0]*2, imSize[1]*2};
	ptrImR 	+= rows[0]*imSize[0];
	ptrX += rows[0]*imSize[0];
	ptrP2DX += rows[0]*imSize[0];

	unsigned char *ptrImTG(ptrImTR + 4*siz), *ptrImTB(ptrImTG + 4*siz), *ptrImG(ptrImR +siz), *ptrImB(ptrImG +siz);
	float *ptrY(ptrX + siz), *ptrZ(ptrY + siz);
	float *ptrP2DY(ptrP2DX + siz);
	
    Vector3f vec, vecRotated;
	float kPar;
	int x0(imSize[0]), y0(imSize[1]), xx, yy;
	
    for (int ii=rows[0]; ii < rows[1]; ii++)
        for (int jj=0; jj < imSize[0]; jj++){
			//Rotation of the sphere
			vec << *ptrX++, *ptrY++, *ptrZ++;
			vecRotated = rotationMatrix * vec;
			//Reverse projection to the plane of the new image. The result is written to an output variable of type CImg	
			kPar = sphere.diameter/(sphere.radius - vecRotated(2));
       		*ptrP2DX = kPar * vecRotated(0);
        	*ptrP2DY = kPar * vecRotated(1);
			//In the new image, the pixel from the input image are copied to their new pixel coordinates
			xx = x0 + int(*ptrP2DX++);
			yy = y0 - int(*ptrP2DY++); 		
			if (xx > 0 && xx < imDoubleSize[0]-1 && yy > 0 && yy < imDoubleSize[1]-1){
				*(ptrImTR + yy * imDoubleSize[0] + xx) = *ptrImR;
				*(ptrImTG + yy * imDoubleSize[0] + xx) = *ptrImG;
				*(ptrImTB + yy * imDoubleSize[0] + xx) = *ptrImB;
			}
		ptrImR++ ; ptrImG++ ; ptrImB++;
		
    }
}
	

void FillQuandrangulars_Recursive(unsigned char *ptrImTR, unsigned char *ptrImR, float *ptrTX, int imSize[2], int rows[2]){
	//Out:		ImTR		Transformed RGB Image 				2Nx2M
	//In:		ImR			Original RGB Image 					NxM
	//In:		TX			Transformed coordinate				NxM
	
	const int 	siz(imSize[0] * imSize[1]);
    int 		doubleSize[2] ={ imSize[0]*2, imSize[1]*2};
	unsigned char *ptrImTG(ptrImTR + 4*siz), *ptrImTB(ptrImTG + 4*siz), *ptrImG(ptrImR +siz), *ptrImB(ptrImG +siz);
	float 		  *ptrTY(ptrTX + siz);
	Coordinates_Sq	quadCoordinates;
	Color4			quadColors;
	int 		cornerIndex, sequence[4] = {0,1,3,2};
	int			xMax(imSize[0]), xMin(-imSize[0]), yMax(imSize[1]), yMin(-imSize[1]);
	
	for (int ii=rows[0]; ii < (rows[1]-1); ii++)
        for (int jj=0; jj < imSize[0]-1; jj++){
			//Define indices of corners
			for (int kk=0; kk<4; kk++){
				cornerIndex = (ii + kk/2) * imSize[0] + jj + kk%2;
				if (cornerIndex <= siz) {
					quadCoordinates(0,sequence[kk]) = (int) *(ptrTX + cornerIndex);
					quadCoordinates(1,sequence[kk]) = (int) *(ptrTY + cornerIndex);
					quadColors(0,sequence[kk]) = *(ptrImR + cornerIndex);
					quadColors(1,sequence[kk]) = *(ptrImG + cornerIndex);
					quadColors(2,sequence[kk]) = *(ptrImB + cornerIndex);
				}
			}
			
			if (quadCoordinates(0,0) > xMin && quadCoordinates(0,0) < xMax && quadCoordinates(1,0) > yMin && quadCoordinates(1,0) < yMax && 
				quadCoordinates(0,1) > xMin && quadCoordinates(0,1) < xMax && quadCoordinates(1,1) > yMin && quadCoordinates(1,1) < yMax &&
				quadCoordinates(0,2) > xMin && quadCoordinates(0,2) < xMax && quadCoordinates(1,2) > yMin && quadCoordinates(1,2) < yMax &&
				quadCoordinates(0,3) > xMin && quadCoordinates(0,3) < xMax && quadCoordinates(1,3) > yMin && quadCoordinates(1,3) < yMax 	) {
				SegmentateQuadrangular(ptrImTR, doubleSize, quadCoordinates, quadColors);
			}
		}
	
}


void SegmentateQuadrangular(unsigned char *ptrImTR, int imSize[2], Coordinates_Sq quadCoord, Color4 quadColors) {
	//In-Out:	ImR		Transformed RGB Image 				2Nx2M	
	//in:		imSize										2Nx2M
	Eigen::Vector2i centerPoint, leftTop;
	int dsiz = imSize[0]*imSize[1];
	int xRef(imSize[0]/2), yRef(imSize[1]/2);
	Coordinates_Sq subQuad, halfPoints;
	Color4		subQuadColors, halfPointsColor;
	Color1 		centerPointColor;
	bool		topSide(abs(quadCoord(0,0) - quadCoord(0,1)) > MIN_SQUARE || abs(quadCoord(1,0) - quadCoord(1,1))> MIN_SQUARE); //true if 'length side' > MIN_SQUARE;
	bool		rightSide(abs(quadCoord(0,1) - quadCoord(0,2)) > MIN_SQUARE || abs(quadCoord(1,1) - quadCoord(1,2))> MIN_SQUARE);
	
	bool		topSide_min(abs(quadCoord(0,0) - quadCoord(0,1)) > MIN_SQUARE -1|| abs(quadCoord(1,0) - quadCoord(1,1))> MIN_SQUARE-1);
	bool		rightSide_min(abs(quadCoord(0,1) - quadCoord(0,2)) > MIN_SQUARE -1 || abs(quadCoord(1,1) - quadCoord(1,2))> MIN_SQUARE -1);
	
	bool		topSide_minmin(abs(quadCoord(0,0) - quadCoord(0,1)) > MIN_SQUARE -2|| abs(quadCoord(1,0) - quadCoord(1,1))> MIN_SQUARE-2);
	bool		rightSide_minmin(abs(quadCoord(0,1) - quadCoord(0,2)) > MIN_SQUARE -2 || abs(quadCoord(1,1) - quadCoord(1,2))> MIN_SQUARE -2);
	
	//Option 1: square with both sides > MIN_SQUARE => segment is split in 4 subsquares
	if (topSide && rightSide) {	
		for (int ii=0; ii< 4; ii++){
			int jj = (ii+1)%4;
			halfPoints.col(ii) = (quadCoord.col(ii) + quadCoord.col(jj) ) / 2;		
			for (int kk=0; kk<3 ; kk++) {
				halfPointsColor(kk,ii) = (quadColors(kk, ii) >> 1) + (quadColors(kk, jj)>>1);
				*(ptrImTR + halfPoints(0,ii) + xRef + (yRef - halfPoints(1,ii))*imSize[0] + kk*dsiz) = halfPointsColor(kk,ii);
			}
		}
		
		centerPoint = (halfPoints.col(0) + halfPoints.col(2)) / 2;
		for (int kk=0; kk<3 ; kk++) {
			centerPointColor(kk) = (halfPointsColor(kk,0) >> 1) + (halfPointsColor(kk, 2) >>1);
			*(ptrImTR + centerPoint(0) + xRef + (yRef - centerPoint(1))*imSize[0] + kk*dsiz) = centerPointColor(kk);
		}
		
		//recurrent call segmentation
		for (int ii=0; ii< 4; ii++){
			subQuad.col(ii) 		= quadCoord.col(ii);				subQuadColors.col(ii) 		= quadColors.col(ii);
			subQuad.col((ii+1)%4) = halfPoints.col(ii);			subQuadColors.col((ii+1)%4) 	= halfPointsColor.col(ii);
			subQuad.col((ii+2)%4) = centerPoint;					subQuadColors.col((ii+2)%4) 	= centerPointColor;
			subQuad.col((ii+3)%4) = halfPoints.col((ii+3)%4);	subQuadColors.col((ii+3)%4)	= halfPointsColor.col((ii+3)%4);
			SegmentateQuadrangular(ptrImTR, imSize, subQuad, subQuadColors);
		}
	}
	
	//Option 2: square with side Sq[0]-Sq[1] > MIN_SQUARE, Sq[1]-Sq[2] < MIN_SQUARE  => segment is split in 2 subsquares
	else if (topSide) {
		for (int ii=0; ii< 2; ii++){
			int jj = 2*ii+1;
			halfPoints.col(ii) = (quadCoord.col(2*ii) + quadCoord.col(jj)) / 2;
			for (int kk=0; kk<3 ; kk++) {
				halfPointsColor(kk, ii) = (quadColors(kk, 2*ii) >> 1) + (quadColors(kk, jj) >>1);
				*(ptrImTR + halfPoints(0,ii) + xRef + (yRef - halfPoints(1,ii))*imSize[0] + kk*dsiz) = halfPointsColor(kk,ii);
			}
		}
		//recurrent call segmentation on subsquares
		for (int ii=0; ii< 2; ii++){
			subQuad.col(0) = quadCoord.col(2*ii);			subQuadColors.col(0) 	= quadColors.col(2*ii);
			subQuad.col(1) = halfPoints.col(ii);			subQuadColors.col(1)	= halfPointsColor.col(ii);
			subQuad.col(2) = halfPoints.col((ii+1)%2);	subQuadColors.col(2) 	= halfPointsColor.col((ii+1)%2);
			subQuad.col(3) = quadCoord.col((2*ii+3)%4);		subQuadColors.col(3) 	= quadColors.col((2*ii+3)%4);
			SegmentateQuadrangular(ptrImTR, imSize, subQuad, subQuadColors);
		}
	}
	
	//Option 3: square with side Sq[0]-Sq[1] <= MIN_SQUARE, Sq[1]-Sq[2] > MIN_SQUARE  => segment is split in 2 subsquares
	else if (rightSide) {
		for (int ii=0; ii< 2; ii++){
			int jj = (2*ii+2)%4;
			halfPoints.col(ii) = (quadCoord.col(2*ii+1) + quadCoord.col(jj)) / 2;
			for (int kk=0; kk<3 ; kk++){
				halfPointsColor(kk, ii) = (quadColors(kk, 2*ii+1) >> 1) + (quadColors(kk, jj) >> 1); 
				*(ptrImTR + halfPoints(0,ii) + xRef + (yRef - halfPoints(1,ii))*imSize[0] + kk*dsiz) = halfPointsColor(kk,ii);
			}
		}
		
		//recurrent call segmentation on subsquares
		for (int ii=0; ii< 2; ii++){
			subQuad.col(0) = quadCoord.col(2*ii);			subQuadColors.col(0) 	= quadColors.col(2*ii);
			subQuad.col(1) = quadCoord.col(2*ii+1);			subQuadColors.col(1) 	= quadColors.col(2*ii+1);
			subQuad.col(2) = halfPoints.col(ii);			subQuadColors.col(2) 	= halfPointsColor.col(ii);
			subQuad.col(3) = halfPoints.col((ii+1)%2);	subQuadColors.col(3) 	= halfPointsColor.col((ii+1)%2);
			SegmentateQuadrangular(ptrImTR, imSize, subQuad, subQuadColors);
		}
	}
	

	//Option 4: square with side Sq[0]-Sq[1] > MIN_SQUARE-1, Sq[1]-Sq[2] >= MIN_SQUARE-1  => center of the square is painted with the average color
	else if (topSide_minmin && rightSide_minmin) {
		for (int kk=0; kk<3 ; kk++) centerPointColor(kk) = (quadColors(kk,0) >> 2) + (quadColors(kk, 1) >>2) + (quadColors(kk,2) >> 2) + (quadColors(kk, 3) >>2);
		leftTop(0) = (4*xRef + quadCoord(0,0) + quadCoord(0,1) + quadCoord(0,2) + quadCoord(0,3))/4;
		leftTop(1) = (4*yRef - quadCoord(1,0) - quadCoord(1,1) - quadCoord(1,2) - quadCoord(1,3))/4;
		
		for(int ii=0; ii<4; ii++) 
			for(int kk=0; kk<3; kk++) *(ptrImTR + leftTop(0) + ii%2 + (leftTop(1) + ii/2)*imSize[0] + kk*dsiz) = centerPointColor(kk);
	}
	
	//Option 5: square with side Sq[0]-Sq[1] > MIN_SQUARE-2, Sq[1]-Sq[2] >= MIN_SQUARE-2  => center of the square is painted with the average color
	else if (topSide_minmin && rightSide_minmin) {
		for (int kk=0; kk<3 ; kk++) centerPointColor(kk) = (quadColors(kk,0) >> 2) + (quadColors(kk, 1) >>2) + (quadColors(kk,2) >> 2) + (quadColors(kk, 3) >>2);
		leftTop(0) = (4*xRef + quadCoord(0,0) + quadCoord(0,1) + quadCoord(0,2) + quadCoord(0,3))/4;
		leftTop(1) = (4*yRef - quadCoord(1,0) - quadCoord(1,1) - quadCoord(1,2) - quadCoord(1,3))/4;
		
		for(int kk=0; kk<3; kk++) *(ptrImTR + leftTop(0)+ leftTop(1)*imSize[0] + kk*dsiz) = centerPointColor(kk);
	}
}
	
void FillQuandrangulars(unsigned char *ptrImTR, unsigned char *ptrImR, float *ptrImTX, int const imSize[2], int rows[2]){
	
	int const imDoubleSize[2] = {imSize[0]*2, imSize[1]*2}, siz(imSize[0] * imSize[1]);
	
	ptrImTX += rows[0]*imSize[0];
	ptrImR += rows[0]*imSize[0];
	float *ptrImTY(ptrImTX + siz);
	unsigned char *ptrImG(ptrImR+siz), *ptrImB(ptrImG+siz);
	int x0(imSize[0]), y0(imSize[1]); 
	int topTriangleCorners[3][2], bottomTriangleCorners[3][2];
	unsigned char topTriangleColors[3][3], bottomTriangleColors[3][3];
	bool paint, inside;
	
	for (int ii=rows[0]; ii < rows[1]-1; ii++){
        for (int jj=0; jj < imSize[0]-1; jj++){
			topTriangleCorners[0][0] = x0 + int(*ptrImTX); 				topTriangleCorners[0][1] = y0 - int(*ptrImTY);
			topTriangleCorners[1][0] = x0 + int(*(ptrImTX+1)); 			topTriangleCorners[1][1] = y0 - int(*(ptrImTY+1));
			topTriangleCorners[2][0] = x0 + int(*(ptrImTX+x0)); 		topTriangleCorners[2][1] = y0 - int(*(ptrImTY+x0));
			bottomTriangleCorners[0][0] = topTriangleCorners[2][0]; 	bottomTriangleCorners[0][1] = topTriangleCorners[2][1];
			bottomTriangleCorners[1][0] = topTriangleCorners[1][0]; 	bottomTriangleCorners[1][1] = topTriangleCorners[1][1];
			bottomTriangleCorners[2][0] = x0 + int(*(ptrImTX++ +x0 +1)); 	bottomTriangleCorners[2][1] = y0 - int(*(ptrImTY++ +x0+1));
			
			paint = (abs(topTriangleCorners[2][0] - topTriangleCorners[1][0]) > 3 || abs(topTriangleCorners[2][1] - topTriangleCorners[1][1]) > 3) && (abs(topTriangleCorners[1][0] - bottomTriangleCorners[2][0]) > 3 || abs(topTriangleCorners[2][0] - topTriangleCorners[1][0]) > 3);
			inside = true;
			
			for (int kk=0; kk<3; kk++){
				inside = inside && topTriangleCorners[kk][0] > 1 && topTriangleCorners[kk][0] < imDoubleSize[0] -1 && topTriangleCorners[kk][1] > 1 && topTriangleCorners[kk][1] < imDoubleSize[1] -1;
			}
			inside = inside && bottomTriangleCorners[2][0] > 1 && bottomTriangleCorners[2][0] < imDoubleSize[0] -1 && bottomTriangleCorners[2][1] > 1 && bottomTriangleCorners[2][1] < imDoubleSize[1] -1;
					
			if (paint && inside){
							
				topTriangleColors[0][0] = *(ptrImR);			topTriangleColors[0][1] = *(ptrImG);			topTriangleColors[0][2] = *(ptrImB);
				topTriangleColors[1][0] = *(ptrImR + 1);		topTriangleColors[1][1] = *(ptrImG + 1);		topTriangleColors[1][2] = *(ptrImB + 1);
				topTriangleColors[2][0] = *(ptrImR + x0);		topTriangleColors[2][1] = *(ptrImG + x0);		topTriangleColors[2][2] = *(ptrImB + x0);
				
				bottomTriangleColors[0][0] = *(ptrImR + x0);	bottomTriangleColors[0][1] = *(ptrImG + x0);	bottomTriangleColors[0][2] = *(ptrImB + x0 );
				bottomTriangleColors[1][0] = *(ptrImR + 1);		bottomTriangleColors[1][1] = *(ptrImG + 1);		bottomTriangleColors[1][2] = *(ptrImB + 1);
				bottomTriangleColors[2][0] = *(ptrImR + x0+1);	bottomTriangleColors[2][1] = *(ptrImG + x0 +1);	bottomTriangleColors[2][2] = *(ptrImB + x0 + 1);
				
				FillTriangle(ptrImTR, imDoubleSize, topTriangleCorners, topTriangleColors);
				FillTriangle(ptrImTR, imDoubleSize, bottomTriangleCorners, bottomTriangleColors);

			}
			ptrImR++; ptrImG++; ptrImB++;
		}
		ptrImTX++; ptrImTY++; ptrImR++; ptrImG++; ptrImB++;
	}
	
}

void FillTriangle(unsigned char *ptrImTR, int const imSize[2], int triangleCorners[3][2], unsigned char triangleColors[3][3]){  
	
	enum sides {AB, BC, CA};
	
	const int sidesAB_BC_CA[3] = {0, 1, 2};
	int sidesHight[3]; //hight of the sides AB, BC and CA
    float displacementVectorsABC[3][2]; //displacement y for 
    int borderPixels[LENGTH_TRIANGLE_ARRAY][3];//borderPixels[y][AB] Contains x-value of pixel of row y along side AB that lies inside the triangle ABC
    int insidePixels[LENGTH_TRIANGLE_ARRAY][2];//insidePixels[y][2] Contains first and last pixels of row y inside triangle ABC
	int longestSide, topRow, triangleHeight, oppositeCorner;//triangle parameters - topRow is the row associated with the top corner (smallest y pixel index)
		
	//Initialize borderPixels[][] and insidePixels[][]
	for (int ii=0; ii<LENGTH_TRIANGLE_ARRAY; ii++){
        borderPixels[ii][0] =-1; borderPixels[ii][1]=-1; borderPixels[ii][2]=-1;
        insidePixels[ii][0] =-1; insidePixels[ii][1]=-1;
    }
     
    //Calculation of "vertical" length (hight) of the sides (measured along y-axis)
    for (auto ii : sidesAB_BC_CA) sidesHight[ii] = triangleCorners[(ii+1)%3][1]-triangleCorners[ii][1];
    
	//Calculate triangle parameters
	longestSide = (abs(sidesHight[AB]) > abs(sidesHight[BC]) )? 0:1;
    longestSide = (abs(sidesHight[CA]) > abs(sidesHight[longestSide]) )? 2: longestSide;
	oppositeCorner =(longestSide+2)%3;
	topRow = (triangleCorners[longestSide][1] < triangleCorners[(longestSide+1)%3][1])?  triangleCorners[longestSide][1] : triangleCorners[(longestSide+1)%3][1];    
    triangleHeight = abs(sidesHight[longestSide]);
	if (LENGTH_TRIANGLE_ARRAY < triangleHeight - 5) std::cout << "LENGTH_TRIANGLE_ARRAY to short !!!!" << std::endl;
	
    //For each triangle side, the displacement of y is calculated upon increment x by 1;
    for (auto ii : sidesAB_BC_CA) {     
        if (sidesHight[ii]!= 0) {
            displacementVectorsABC[ii][0] = float(triangleCorners[(ii+1)%3][0]-triangleCorners[ii][0])/float(sidesHight[ii]);
            displacementVectorsABC[ii][1] = float(sidesHight[ii])/abs(sidesHight[ii]);
        }
        else {
            displacementVectorsABC[ii][0] = 10000.;
            displacementVectorsABC[ii][1] = -1;
        }
    }
 
    //Search borderPixels at the border sides AB, BC and CA inside the triangle. 
    for (auto kk : sidesAB_BC_CA) {
        float xx(triangleCorners[kk][0]);
        if (displacementVectorsABC[kk][1] >0)
            for (int ii = triangleCorners[kk][1]+1; ii < triangleCorners[(kk+1)%3][1]; ii++){
                xx += displacementVectorsABC[kk][0];
                borderPixels[ii-topRow][kk]= int(floor(xx));
            }
        else
            for (int ii = triangleCorners[kk][1]-1; ii > triangleCorners[(kk+1)%3][1]; ii--){
                xx -= displacementVectorsABC[kk][0];
                borderPixels[ii-topRow][kk] = int(ceil(xx));
            }
 
        borderPixels[triangleCorners[kk][1]-topRow][kk] = triangleCorners[kk][0];
        borderPixels[triangleCorners[(kk+1)%3][1]-topRow][kk] = triangleCorners[(kk+1)%3][0];  
    }       
     
    //Find first and last pixel in each row inside the triangle - write to insidePixels[y][]
    for (int ii = 0; ii< triangleHeight+1; ii++) {
        if (triangleCorners[oppositeCorner][0] >  borderPixels[triangleCorners[oppositeCorner][1] - topRow][longestSide]) {
            insidePixels[ii][0] = borderPixels[ii][longestSide];
            insidePixels[ii][1] = (borderPixels[ii][(longestSide+1)%3]> borderPixels[ii][(longestSide+2)%3])? borderPixels[ii][(longestSide+1)%3] : borderPixels[ii][(longestSide+2)%3];
        }
        else {
            insidePixels[ii][1] = borderPixels[ii][longestSide];
            insidePixels[ii][0] = (borderPixels[ii][(longestSide+1)%3] > borderPixels[ii][(longestSide+2)%3])? borderPixels[ii][(longestSide+1)%3] : borderPixels[ii][(longestSide+2)%3];
        }
        if (triangleCorners[oppositeCorner][0] >  borderPixels[triangleCorners[oppositeCorner][1] - topRow][longestSide]) 
            insidePixels[triangleCorners[oppositeCorner][1] - topRow][1] =  triangleCorners[oppositeCorner][0]-1;
        else
            insidePixels[triangleCorners[oppositeCorner][1] - topRow][0] =  triangleCorners[oppositeCorner][0]+1;
         
    }
     
    //Accidentally, the x values of insidePixels can be reversed
    int dummyVar;
    for (int ii = 1; ii< triangleHeight; ii++) 
        if (insidePixels[ii][0] > insidePixels[ii][1]){
            dummyVar 			= insidePixels[ii][1];
            insidePixels[ii][1] = insidePixels[ii][0];
            insidePixels[ii][0] = dummyVar;
        }
       
    //Declare variable colour gradients along the sides of the triangle
    double colorGradient[3][2];
    for (int ii=0; ii<3; ii++){
        colorGradient[ii][0] = 0.0;
        colorGradient[ii][1] = 0.0;
    }
     
    //Reassign corner coordinates; corner B wil be the top (smallest y-value)
    int help;
	int cornersReAssigned[3][2];			//cornersReAssigned: Top of the triangle is corner 'B' / with ref '0', 'A' = 2, 'C' = 1  !!!
	unsigned char colorsReAssigned[3][3]; 	//colorsReAssigned : Top of the triangle is corner 'B' / with ref '0'!!!
    if (triangleCorners[0][1] < triangleCorners[1][1])
        if (triangleCorners[0][1] < triangleCorners[2][1])	help = 0;
		else 	help = 2;
    else 		help = 1;
	
	for (auto ii : sidesAB_BC_CA){
		cornersReAssigned[ii][0] = triangleCorners[(help+ii)%3][0]; 
		cornersReAssigned[ii][1] = triangleCorners[(help+ii)%3][1];
		for (auto kk : sidesAB_BC_CA) colorsReAssigned[ii][kk] = triangleColors[(help+ii)%3][kk];
	}
 
    //Calculate colorGradientients
    int xAB(cornersReAssigned[2][0]-cornersReAssigned[0][0]), yAB(cornersReAssigned[2][1]-cornersReAssigned[0][1]), xBC(cornersReAssigned[0][0]-cornersReAssigned[1][0]), yBC(cornersReAssigned[0][1]-cornersReAssigned[1][1]);
    double Denominator(double(xAB*yBC - xBC*yAB));
	
	for (int ii=0; ii<3; ii++){
        colorGradient[ii][0] = ((int(colorsReAssigned[2][ii]) - int(colorsReAssigned[0][ii])) * yBC - (int(colorsReAssigned[0][ii]) - int(colorsReAssigned[1][ii])) * yAB)/Denominator;
        colorGradient[ii][1] = ((int(colorsReAssigned[0][ii]) - int(colorsReAssigned[1][ii])) * xAB - (int(colorsReAssigned[2][ii]) - int(colorsReAssigned[0][ii])) * xBC)/Denominator;
    }
     
    //Calculate pixel RGB values
	int Size(imSize[0]*imSize[1]), RelativePosition, x2B, y2B;
	unsigned char *ptrImTG(ptrImTR+Size), *ptrImTB(ptrImTG+Size); 
    for (int ii = 0; ii< triangleHeight; ii++)
        for (int jj = insidePixels[ii][0]; jj<=insidePixels[ii][1]; jj++){
			RelativePosition = (ii + topRow)*imSize[0] + jj;
			x2B = cornersReAssigned[0][0] - jj;
			y2B = cornersReAssigned[0][1] - ii - topRow;
            *(ptrImTR + RelativePosition) = int(int(colorsReAssigned[0][0]) - x2B*colorGradient[0][0] - y2B*colorGradient[0][1]);
			*(ptrImTG + RelativePosition) = int(int(colorsReAssigned[0][1]) - x2B*colorGradient[1][0] - y2B*colorGradient[1][1]);
			*(ptrImTB + RelativePosition) = int(int(colorsReAssigned[0][2]) - x2B*colorGradient[2][0] - y2B*colorGradient[2][1]);
        }
}


void FillMissingPixels(unsigned char *ptrImTR, int const doubleSize[2], int rows[2]){
	int dsiz(doubleSize[0] * doubleSize[1]);
	unsigned char *ptrImR, *ptrImG, *ptrImB;
	
	//color black pixel [0, 0, 0] if pixel above has a color
	ptrImR = ptrImTR + (rows[1]+1)*doubleSize[0]-1; ptrImG = ptrImR+dsiz, ptrImB=ptrImG+dsiz; //pointer to last pixel left bottom of image /rowRange
	for (int ii=rows[0]+1; ii<rows[1]; ii++)
		for (int jj = 0; jj <doubleSize[0]; jj++) {
			if ((*ptrImR ==0 && *ptrImG ==0 && *ptrImB ==0) && ( *(ptrImR - doubleSize[0]) > 0 || *(ptrImG  - doubleSize[0]) > 0 || *(ptrImB  - doubleSize[0]) > 0)) {
				 *ptrImR = *(ptrImR - doubleSize[0]);
				 *ptrImG = *(ptrImG - doubleSize[0]);
				 *ptrImB = *(ptrImB - doubleSize[0]);
			}
			ptrImR--; ptrImG--; ptrImB--;
		}
	
	//color black pixel [0, 0, 0] if left pixel has a color
	ptrImR = ptrImTR + (rows[1]+1)*doubleSize[0]-1; ptrImG = ptrImR+dsiz, ptrImB=ptrImG+dsiz; //pointer to last pixel left bottom of image / rowRange
	for (int ii=rows[0]+1; ii<rows[1]; ii++)
		for (int jj = 1; jj <doubleSize[0]; jj++) {
			if ((*ptrImR ==0 && *ptrImG ==0 && *ptrImB ==0) && ( *(ptrImR - 1) > 0 || *(ptrImG  - 1) > 0 || *(ptrImB  - 1) > 0)) {
				 *ptrImR = *(ptrImR - 1);
				 *ptrImG = *(ptrImG - 1);
				 *ptrImB = *(ptrImB - 1);
			}
			ptrImR--; ptrImG--; ptrImB--;
		}
	
	//color black pixel [0, 0, 0] if pixel below has a color
	ptrImR = ptrImTR + rows[0]*doubleSize[0]; ptrImG = ptrImR+dsiz, ptrImB=ptrImG+dsiz;
	for (int ii=rows[0]; ii<rows[1]-1; ii++)
		for (int jj = 0; jj <doubleSize[0]; jj++) {
			if ((*ptrImR ==0 && *ptrImG ==0 && *ptrImB ==0) && ( *(ptrImR + doubleSize[0]) > 0 || *(ptrImG  + doubleSize[0]) > 0 || *(ptrImB  + doubleSize[0]) > 0)) {
				 *ptrImR = *(ptrImR + doubleSize[0]);
				 *ptrImG = *(ptrImG + doubleSize[0]);
				 *ptrImB = *(ptrImB + doubleSize[0]);
			}
			ptrImR++; ptrImG++; ptrImB++;
		}	 
	
	//color black pixel [0, 0, 0] if right pixel has a color
	ptrImR = ptrImTR + rows[0]*doubleSize[0]; ptrImG = ptrImR+dsiz, ptrImB=ptrImG+dsiz;
	for (int ii=rows[0]; ii<rows[1]; ii++){
		for (int jj = 0; jj <doubleSize[0]-1; jj++) {
			if ((*ptrImR ==0 && *ptrImG ==0 && *ptrImB ==0) && ( *(ptrImR + 1) > 0 || *(ptrImG  + 1) > 0 || *(ptrImB  + 1) > 0)) {
				 *ptrImR = *(ptrImR + 1);
				 *ptrImG = *(ptrImG + 1);
				 *ptrImB = *(ptrImB + 1);
			}
			ptrImR++; ptrImG++; ptrImB++;
		}
	 	ptrImR++; ptrImG++; ptrImB++;
	}		 
}

