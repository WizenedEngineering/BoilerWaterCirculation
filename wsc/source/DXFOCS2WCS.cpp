/*****************************************************************//**
 * \file   OCS2WCS.cpp
 *
 * \author Rainer_Jordan@<very, very warm>mail.com
 * ## Licence
 * Licensed under the EUPL, Version 1.2 or later
 * \date   August 2021
 *********************************************************************/
 //#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

using namespace std;

 /**
 * \brief converts OCS data of arcs to WCS data.
 *
 * \param OCSCenterX x-coordinate of center point in OCS
 * \param OCSCenterY y-coordinate of center point in OCS
 * \param OCSCenterZ z-coordinate of center point in OCS
 * \param Nx x-component of extrusion vector
 * \param Ny y-component of extrusion vector
 * \param Nz z-component of extrusion vector
 * \param OCSRadius radius of arc
 * \param OCSAngle angle of point in OCS
 * \return
 */
_point OCS2WCS(double OCSCenterX, double OCSCenterY, double OCSCenterZ,
	double Nx, double Ny, double Nz,
	double OCSRadius, double OCSAngle) {

	// Autocad help gives a limit of 1/64 for Nx and Ny to switch to different calculation
	// otherwise division by 0
	_point WCS;

	if (fabs(Nx) < 1. / 64. && fabs(Ny) < 1. / 64.) {
		double root = sqrt(1. - Ny * Ny);
		WCS.xCoord = OCSCenterX * Nz / root - OCSCenterY * Ny * Nx / root + OCSCenterZ * Nx;
		WCS.yCoord = OCSCenterY * root + OCSCenterZ * Ny;
		WCS.zCoord = -OCSCenterX * Nx / root - OCSCenterY * Ny * Nz / root + OCSCenterZ * Nz;
		if (OCSRadius < 1e-6) return WCS;
		double OCSx = cos(OCSAngle / 180. * M_PI) * OCSRadius;
		double OCSy = sin(OCSAngle / 180. * M_PI) * OCSRadius;
		WCS.xCoord += OCSx * Nz / root - OCSy * Ny * Nx / root;
		WCS.yCoord += OCSy * root;
		WCS.zCoord += -OCSx * Nx / root - OCSy * Ny * Nz / root;
	}
	else {
		double root = sqrt(1. - Nz * Nz);
		WCS.xCoord = -OCSCenterX * Ny / root - OCSCenterY * Nz * Nx / root + OCSCenterZ * Nx;
		WCS.yCoord = OCSCenterX * Nx / root - OCSCenterY * Nz * Ny / root + OCSCenterZ * Ny;
		WCS.zCoord = OCSCenterY * root + OCSCenterZ * Nz;
		if (OCSRadius < 1e-6) return WCS;
		double OCSx = cos(OCSAngle / 180. * M_PI) * OCSRadius;
		double OCSy = sin(OCSAngle / 180. * M_PI) * OCSRadius;
		WCS.xCoord += -OCSx * Ny / root - OCSy * Nz * Nx / root;
		WCS.yCoord += OCSx * Nx / root - OCSy * Nz * Ny / root;
		WCS.zCoord += OCSy * root;
	}
	return WCS;
}
