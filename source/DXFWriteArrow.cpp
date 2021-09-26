/*****************************************************************//**
 * \file   DXFWriteArrow.cpp
 *
 * \author Rainer_Jordan@<very, very warm>mail.com
 * \date   5. November 2017
 *********************************************************************/
 //#include "stdafx.h"

#include <cstdlib>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

extern int DXFWriteLine(ostream& outData, unsigned long long& handle,
	const string& LayerName, unsigned short color, double xStart, double yStart,
	double zStart, double xEnd, double yEnd, double zEnd);

/**
* \fn int TransformWrite(ostream& outData, double x0, double y0, double z0,
	double xMidPoint, double yMidPoint, double zMidPoint, double sin_phiX, double cos_phiX,
	double sin_phiZ, double cos_phiZ,
	unsigned long long& handle, const string& LayerName, unsigned short color) 
* \brief transforms arrow line to direction of tube
*
* \param outData file handle
* \param x0 x-coordinate of arrow line second point
* \param y0 y-coordinate of arrow line second point
* \param z0 z-coordinate of arrow line second point
* \param xMidPoint x-coordinate of line (tube) center point, it's the tip of arrow
* \param yMidPoint y-coordinate of line (tube) center point, it's the tip of arrow
* \param zMidPoint z-coordinate of line (tube) center point, it's the tip of arrow
* \param sin_phiX sinus of angle line (tube) to x-plane
* \param cos_phiX cosinus of angle line (tube) to x-plane
* \param sin_phiZ sinus of angle line (tube) to z-plane
* \param cos_phiZ cosinus of angle line (tube) to z-plane
* \param handle graphic object handle (has to be a unique number)
* \param LayerName layer name
* \param color according CAD numbering scheme
* \return int error code
*/
int TransformWrite(ostream& outData, double x0, double y0, double z0,
	double xMidPoint, double yMidPoint, double zMidPoint, double sin_phiX, double cos_phiX,
	double sin_phiZ, double cos_phiZ,
	unsigned long long& handle, const string& LayerName, unsigned short color) {
	double y1, xE, yE, zE;
	//rotation around z-Axis and x-Axis
	y1 = x0 * sin_phiZ + y0 * cos_phiZ;
	xE = xMidPoint + x0 * cos_phiZ - y0 * sin_phiZ;
	yE = yMidPoint + y1 * cos_phiX - z0 * sin_phiX;
	zE = zMidPoint + y1 * sin_phiX + z0 * cos_phiX;

	return DXFWriteLine(outData, handle, LayerName, color,
		xMidPoint, yMidPoint, zMidPoint, xE, yE, zE);
}

/**
 *int DXFWriteArrow(ostream& outData, unsigned long long& handle, unsigned short color,
        double xMidPoint, double yMidPoint, double zMidPoint, double dx, double dy, double dz)
* \brief Shows an arrow at middle of line/arc to indicate direction of line (from line start point to line end point)
*
* Layer is always "Arrow"\n
* The baseline itself has to be handled elsewhere\n
* An arrow consists of 8 lines to be recognized in 3D view\n
* in this function the orientation is in direction of x-axis. it will be transformed in TransformWrite
* \param outData File handle for .dxf file to be written to
* \param handle Graphic object handle
* \param color of arrow
* \param xMidPoint x-coordinate of line (tube) mid point [mm]
* \param yMidPoint y-coordinate of line (tube) mid point [mm]
* \param zMidPoint z-coordinate of line (tube) mid point [mm]
* \param dx difference of x-coordinate of line (tube) start and end point [mm]
* \param dy difference of y-coordinate of line (tube) start and end point [mm]
* \param dz difference of z-coordinate of line (tube) start and end point [mm]
* \return int error code
*/
int DXFWriteArrow(ostream& outData, unsigned long long& handle, unsigned short color,
        double xMidPoint, double yMidPoint, double zMidPoint, double dx, double dy, double dz) {

	double length, lsin15;
	double x0, y0, z0;
	double sin_phiX, cos_phiX, sin_phiZ, cos_phiZ;

	string LayerName = "Arrow";
	int error = 0;
	const double sin15 = 0.2588, cos15 = 0.9659, sin30 = 0.5, cos30 = 0.866,
		sin60 = 0.866, cos60 = 0.5;

	length = sqrt(dx * dx + dy * dy + dz * dz);
	if (length < 1e-6) return 0;

	if (fabs(dy) < 1e-6 && fabs(dz) < 1e-6) { // x-axis
		sin_phiX = 0.;
		cos_phiX = 1.;

		sin_phiZ = 0.;
		cos_phiZ = dx / fabs(dx);
	}
	else {
		// rotation to xy plane
		double y1 = sqrt(dz * dz + dy * dy);
		sin_phiX = dz / y1;
		cos_phiX = dy / y1;

		sin_phiZ = y1 / length;
		cos_phiZ = dx / length;
	}

	length /= 20.;
	lsin15 = length * sin15;

	x0 = -length * cos15;
	y0 = lsin15;
	z0 = 0.;
	error = TransformWrite(outData, x0, y0, z0, 
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;
	y0 = -lsin15;
	error = TransformWrite(outData, x0, y0, z0, 
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;
	//-----------------
	y0 = 0.;
	z0 = lsin15;
	error = TransformWrite(outData, x0, y0, z0, 
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;

	z0 = -z0; //-length * sin15;
	error = TransformWrite(outData, x0, y0, z0, 
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;
	//--------------------------------
	y0 = lsin15 * cos30;
	z0 = lsin15 * sin30;
	error = TransformWrite(outData, x0, y0, z0, 
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;

	y0 = -y0; //-length * sin15*cos30;
	error = TransformWrite(outData, x0, y0, z0, 
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;
	//---------------------------
	y0 = lsin15 * cos60;
	z0 = lsin15 * sin60;
	error = TransformWrite(outData, x0, y0, z0, 
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;

	y0 = -y0; //-length * sin15*cos60;
	error = TransformWrite(outData, x0, y0, z0,
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;
	//---------------------------

	y0 = lsin15 * cos30;
	z0 = -lsin15 * sin30;
	error = TransformWrite(outData, x0, y0, z0, 
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
//	handle++;
	if (error) return error;

	y0 = -y0; //-length * sin15*cos30;
	error = TransformWrite(outData, x0, y0, z0,
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;
	//---------------------------
	y0 = lsin15 * cos60;
	z0 = -lsin15 * sin60;
	error = TransformWrite(outData, x0, y0, z0, 
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;

	y0 = -y0; //-length * sin15*cos60;
	error = TransformWrite(outData, x0, y0, z0, 
		xMidPoint, yMidPoint, zMidPoint,
		sin_phiX, cos_phiX, sin_phiZ, cos_phiZ,
		handle, LayerName, color);
	if (error) return error;

	return 0;
}
