/*****************************************************************//**
 * \file   DXFWriteArc.cpp
 * \author Rainer_Jordan@<very, very warm>mail.com
 * ## Licence
 * Licensed under the EUPL, Version 1.2 or later
 * \date   August 2021
 *********************************************************************/
//#include "stdafx.h"
#undef MAINFUNCTION
//#include "CommonHeader.h"
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

/*!
*\fn int DXFWriteArc(ostream& outData, unsigned long long& handle, const string& layer, unsigned short color,
   double OCSCenterX, double OCSCenterY, double OCSCenterZ, double Radius, double Nx, double Ny, double Nz,
   double OCSStartAngle, double OSCEndAngle)
*  \brief writes data for an arc to .dxf file
*
* \param outData File handle for .dxf file to be written to
* \param handle Graphic object handle (has to be unique number)
* \param layer Layer name
* \param color color of line
* \param OCSCenterX x-coordinate of center point (OCS)
* \param OCSCenterY y-coordinate of center point (OCS)
* \param OCSCenterZ z-coordinate of center point (OCS)
* \param Radius x-coordinate of end point
* \param Nx x-component of extrusion vector
* \param Ny y-component of extrusion vector
* \param Nz z-component of extrusion vector
* \param OCSStartAngle start angle in OCS
* \param OSCEndAngle end angle in OCS
* \return int error code
*/
int DXFWriteArc(ostream& outData, unsigned long long& handle, const string& layer, unsigned short color,
   double OCSCenterX, double OCSCenterY, double OCSCenterZ, double Radius, double Nx, double Ny, double Nz,
   double OCSStartAngle, double OSCEndAngle) {

   outData << "  0\nARC\n  5\n";
   outData << std::hex << ++handle << "\n100\nAcDbEntity\n  8\n" << layer << "\n 62\n";
   outData << std::dec << color;
   outData << "\n100\nAcDbCircle\n 10\n" 
      << setw(18) << setprecision(15) << OCSCenterX << "\n 20\n"
      << setw(18) << setprecision(15) << OCSCenterY << "\n 30\n"
      << setw(18) << setprecision(15) << OCSCenterZ << "\n";
   outData << " 40\n" << Radius << "\n210\n" << Nx << "\n220\n" << Ny << "\n230\n" << Nz << "\n100\nAcDbArc\n 50\n" 
      << OCSStartAngle << "\n 51\n" << OSCEndAngle << "\n";
   return 0;
}
