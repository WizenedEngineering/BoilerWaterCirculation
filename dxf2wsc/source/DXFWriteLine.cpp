/*****************************************************************//**
 * \file   DXFWriteLine.cpp
 * \author rainer_jordan-at-<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
 * \date   November 2017
 *********************************************************************/
#include "common.h"

using namespace std;

/*!
* \brief writes data for a line to .dxf file
*
* \param outData File handle for .dxf file to be written to
* \param handle Graphic object handle (has to be unique number)
* \param layer Layer name
* \param color color of line
* \param xStart x-coordinate of start point
* \param yStart y-coordinate of start point
* \param zStart z-coordinate of start point
* \param xEnd x-coordinate of end point
* \param yEnd y-coordinate of end point
* \param zEnd z-coordinate of end point
* \return int error code
*/
int DXFWriteLine(ostream& outData, unsigned long long& handle, const string& layer, unsigned short color,
        double xStart, double yStart, double zStart, double xEnd, double yEnd, double zEnd) {

   outData << "  0\nLINE\n  5\n";
   outData << std::hex << ++handle << "\n100\nAcDbEntity\n  8\n" << layer << "\n 62\n" ;
   outData<< std::dec<<color;
   outData << "\n100\nAcDbLine\n 10\n" << xStart << "\n 20\n" << yStart << "\n 30\n" << zStart << "\n";
   outData << " 11\n" << xEnd << "\n 21\n" << yEnd << "\n 31\n" << zEnd << "\n";
   return 0;
}
