/*****************************************************************//**
 * \file   DXFWriteTube.cpp
 * \author Rainer_Jordan@<very, very warm>mail.com
 * ## Licence
 * Licensed under the EUPL, Version 1.2 or later
* \date   November 2017
 *********************************************************************/
 /** 
 * \fn int DXFWriteTube(ostream& outData, _tube& iTube,  unsigned long long& handle,
	const string& LayerName, unsigned short color, double& PointMidX, double& PointMidY,
	double& PointMidZ, double& dx, double& dy, double& dz)
 * \brief writing tube (straight or bend) geometry to .dxf file 
 *
 * \param [in] outData outstream of .dxf file
 * \param [in] iTube reference of tube
 * \param [in, out] handle Graphic object handle (has to be unique number)
 * \param [in] LayerName Layer name
 * \param [in] color line color
 * \param [out] PointMidX x-coordinate of mid point
 * \param [out] PointMidY x-coordinate of mid point
 * \param [out] PointMidZ x-coordinate of mid point
 * \param [out] dx difference of x-coordinate between start and end point 
 * \param [out] dy difference of y-coordinate between start and end point
 * \param [out] dz difference of z-coordinate between start and end point
 * \return int error
 */
//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

using namespace std;

extern int DXFWriteLine(ostream& outData, unsigned long long& handle,
	const string& LayerName, unsigned short color, double xStart, double yStart,
	double zStart, double xEnd, double yEnd, double zEnd);

extern int DXFWriteArc(ostream& outData, unsigned long long& handle, const string& layer, unsigned short color,
	double OCSCenterX, double OCSCenterY, double OCSCenterZ, double Radius, double Nx, double Ny, double Nz,
	double OCSStartAngle, double OCSEndAngle);

extern _point OCS2WCS(double OCSCenterX, double OCSCenterY, double OCSCenterZ,
	double Nx, double Ny, double Nz,
	double OCSRadius, double OCSAngle);
 
extern int DXFWriteArrow(ostream& outData, unsigned long long& handle,
	unsigned short color, double PointMidX, double PointMidY,
	double PointMidZ, double dx, double dy, double dz);

int DXFWriteTube(ostream& outData, const _tube& iTube,  unsigned long long& handle,
	const string& LayerName, unsigned short color, bool writeArrow) {
	double PointMidX,
	PointMidY,
	PointMidZ,
		OCSMidAngle;

   double dx, dy, dz;

	int error = 0;
	_point* PtIn = &Points[iTube.PointIn];
	_point* PtOut= &Points[iTube.PointOut];

	if (iTube.RadiusBend < 1e-3) { //straight tube
		error = DXFWriteLine(outData, handle, LayerName, color,
			PtIn->xCoord, PtIn->yCoord, PtIn->zCoord,
			PtOut->xCoord, PtOut->yCoord, PtOut->zCoord);
		if (error) return error;
		PointMidX = (PtIn->xCoord + PtOut->xCoord) / 2.;
		PointMidY = (PtIn->yCoord + PtOut->yCoord) / 2.;
		PointMidZ = (PtIn->zCoord + PtOut->zCoord) / 2.;
	}
	else { //bend
		error = DXFWriteArc(outData, handle, LayerName, color,
			iTube.OCSCenterX, iTube.OCSCenterY, iTube.OCSCenterZ,
			iTube.RadiusBend*1e3, iTube.Nx, iTube.Ny, iTube.Nz,
			iTube.OCSStartAngle, iTube.OCSEndAngle);
		if (error) return error;
		if (iTube.OCSEndAngle - iTube.OCSStartAngle < 0.) {
			OCSMidAngle = (iTube.OCSStartAngle + iTube.OCSEndAngle + 360.) / 2.;
			if (OCSMidAngle > 360.) OCSMidAngle -= 360.;
		}
		else {
			OCSMidAngle = (iTube.OCSStartAngle + iTube.OCSEndAngle) / 2.;
		}
		_point WCSMidPoint = OCS2WCS(iTube.OCSCenterX, iTube.OCSCenterY, iTube.OCSCenterZ,
			iTube.Nx, iTube.Ny, iTube.Nz, iTube.RadiusBend*1e3, OCSMidAngle);
		PointMidX = WCSMidPoint.xCoord;
		PointMidY = WCSMidPoint.yCoord;
		PointMidZ = WCSMidPoint.zCoord;
	}
	dx = PtOut->xCoord - PtIn->xCoord;
	dy = PtOut->yCoord - PtIn->yCoord;
	dz = PtOut->zCoord - PtIn->zCoord;

	if (writeArrow) {
		error = DXFWriteArrow(outData, handle, color,
			PointMidX, PointMidY, PointMidZ, dx, dy, dz);
		if (error) return error;
	}
	return error;
}

