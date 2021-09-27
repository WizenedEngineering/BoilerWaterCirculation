/*****************************************************************//**
 * \file   DXFWritePoint.cpp
 * \author rainer_jordan-at-<very, very warm>mail.com
 * #### Licence
 * Licensed under the European Union Public Licence (EUPL), Version 1.2 or later* \author rainer
 * \date   16. November 2017
 *********************************************************************/

#include "common.h"

using namespace std;

extern int DXFWriteLine(ostream& outData, unsigned long long& handle,
        const string& LayerName, unsigned short color, double xStart, double yStart,
        double zStart, double xEnd, double yEnd, double zEnd);

/*!
* \brief Indicates a point by a small 3D cross
*
* In Autocad 2015 the points are indicated by crosshairs\n
* Those crosshairs can be bigger than the lines, looks ugly\n
*
* Points are drawn as small 3D crosses (each direction 10 mm)\n
* \param outData File handle for .dxf file to be written to
* \param handle Graphic object handle (has to be unique number)
* \param LayerName Layer name
* \param color color of point
* \param xPos x-coordinate of point
* \param yPos y-coordinate of point
* \param zPos z-coordinate of point
* \return int error code
*/
int DXFWritePoint(ostream& outData, unsigned long long& handle, const string& LayerName, unsigned short color,
        double xPos, double yPos, double zPos) {
   //
   DXFWriteLine(outData, handle, LayerName, color,
           xPos, yPos, zPos, xPos + 10., yPos, zPos);
   DXFWriteLine(outData, handle, LayerName, color,
           xPos, yPos, zPos, xPos - 10., yPos, zPos);
   DXFWriteLine(outData, handle, LayerName, color,
           xPos, yPos, zPos, xPos, yPos + 10., zPos);
   DXFWriteLine(outData, handle, LayerName, color,
           xPos, yPos, zPos, xPos, yPos - 10., zPos);
   DXFWriteLine(outData, handle, LayerName, color,
           xPos, yPos, zPos, xPos, yPos, zPos + 10.);
   DXFWriteLine(outData, handle, LayerName, color,
           xPos, yPos, zPos, xPos, yPos, zPos - 10.);

   return 0;
}
