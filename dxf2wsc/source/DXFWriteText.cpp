/*****************************************************************//**
 * \file   DXFWriteText.cpp
 * \author rainer_jordan-at-<very, very warm>mail.com
 * #### Licence
 * Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
 * \date   November 2017
 * \brief writes text at Position on x-z plane to .dxf
*********************************************************************/
#include "common.h"

using namespace std;

/*!
* \brief writes data for a text to .dxf file
*
* \param outData File handle for .dxf file to be written to
* \param handle Graphic object handle (has to be unique number)
* \param layer Layer name
* \param color color of text
* \param xPos x-coordinate of position
* \param yPos y-coordinate of position
* \param zPos z-coordinate of position
* \param height text height 
* \param content text to be written
* \return int error code
*/
int DXFWriteText(ostream& outData, unsigned long long& handle, const string& layer, unsigned short color,
        double xPos, double yPos, double zPos, double height, const string& content) {

   outData << "  0\nTEXT\n  5\n";
   outData << std::hex << ++handle << "\n100\nAcDbEntity\n  8\n" << layer << "\n 62\n";
   outData << std::dec << color;
   // text is written in x-z plane -> yPos and zPos have to be changed 
   outData << "\n100\nAcDbText\n 10\n" << xPos << "\n 20\n" << zPos << "\n 30\n" << -yPos << "\n";
   outData << " 40\n" << height << "\n  1\n" << content << "\n";
   outData << "210\n 0.\n220\n -1.\n230\n 0.\n";
   outData << "100\nAcDbText\n";
   return 0;
}
