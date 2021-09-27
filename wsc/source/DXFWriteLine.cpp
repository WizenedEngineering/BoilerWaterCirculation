/*****************************************************************//**
 * \file   DXFWriteLine.cpp
 * \author Rainer_Jordan@<very, very warm>mail.com
 * ## Licence
 * Licensed under the EUPL, Version 1.2 or later
 * \date   November 2017
 *********************************************************************/
 //#include "stdafx.h"
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
using namespace std;

/*
* \brief writes data for a line to .dxf file
*
* \param outData File handle for .dxf file to be written to
* \param handle Graphic object handle (has to be unique number)
* \param LayerName Layer name
* \param color color of line
* \param xStart x-coordinate of start point
* \param yStart y-coordinate of start point
* \param zStart z-coordinate of start point
* \param xEnd x-coordinate of end point
* \param yEnd y-coordinate of end point
* \param zEnd z-coordinate of end point
* \return int error code
*/
int DXFWriteLine(ostream& outData, unsigned long long& handle, const string& LayerName, unsigned short color,
	double xStart, double yStart, double zStart, double xEnd, double yEnd, double zEnd) {

	outData << "  0\nLINE\n  5\n";
	outData << std::hex << ++handle << "\n100\nAcDbEntity\n  8\n" << LayerName << "\n 62\n" << std::dec << color;
	outData << "\n100\nAcDbLine\n 10\n" 
		<< setw(18) << setprecision(15) << xStart << "\n 20\n"
		<< setw(18) << setprecision(15) << yStart << "\n 30\n"
		<< setw(18) << setprecision(15) << zStart << "\n";
	outData << " 11\n" 
		<< setw(18) << setprecision(15) << xEnd << "\n 21\n"
		<< setw(18) << setprecision(15) << yEnd << "\n 31\n"
		<< setw(18) << setprecision(15) << zEnd << "\n";
	return 0;
}
