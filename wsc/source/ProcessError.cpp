/*****************************************************************//**
 * \file ProcessError.cpp
 * \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*********************************************************************/
//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"
#include <list>

extern int DXFWriteLine(ostream& outData, unsigned long long& handle,
	const string& LayerName, unsigned short color, double xStart, double yStart,
	double zStart, double xEnd, double yEnd, double zEnd);


void processError(size_t iPt, const string& text) {
	ofstream outData;
	string  writefile = PathFile + "PathError.dxf";
	cout << "\nerror on " << text << ". Starting point " << iPt << ". Please check " << writefile << " for position ";
	prot << "\nerror on " << text << ". Starting point " << iPt << ". Please check " << writefile << " for position ";

	outData.open(writefile.c_str());

	if (!outData.good()) {
		cout << "\n cannot open file " << writefile << endl;
	}
	else {
		outData << "999\nWSC_Error\n  0\nSECTION\n  2\nCLASSES\n  0\nENDSEC\n  0\nSECTION\n  2\nENTITIES\n";
		unsigned long long handle = 100ULL;
		unsigned short color = 6;
		DXFWriteLine(outData, handle, "PathError", color,
			100000, 100000, 100000, Points[iPt].xCoord, Points[iPt].yCoord, Points[iPt].zCoord);

		outData << "  0\nENDSEC\n  0\nSECTION\n  2\nOBJECTS\n  0\nENDSEC\n  0\nEOF" << std::endl;
		outData.close();
	}//file opened correctly
	prot.close();
	exit(1);
}
