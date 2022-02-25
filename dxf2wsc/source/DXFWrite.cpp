/*!
 * \file DXFWrite.cpp
 * \fn int DXFWrite ( ostream& outData, vector <_tube> &Tubes, size_t mTb, vector <_point> &Points, size_t mPt)
 * \brief writing mesh geometry (tubes and points) to .dxf file
 * \param outData outstream of .dxf file
 * \param &Tubes vector of tubes
 * \param mTb max. number of tubes
 * \param &Points vector of points
 * \param mPt max number of points
 * \return int error code
 * \author rainer_jordan-at-<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
 * \date   February 2021
 *********************************************************************/
#include "common.h"

using namespace std;

extern int DXFWriteLine(ostream& outData, unsigned long long& handle, const string& layer,
	unsigned short color, double xStart, double yStart, double zStart,
	double xEnd, double yEnd, double zEnd);

extern int DXFWriteArc(ostream& outData, unsigned long long& handle, const string& layer, unsigned short color,
	double OCSCenterX, double OCSCenterY, double OCSCenterZ, double Radius, double Nx, double Ny, double Nz,
	double OCSStartAngle, double OCSEndAngle);

extern int DXFWritePoint(ostream& outData, unsigned long long& handle, const string& LayerName,
	unsigned short color, double xPos, double yPos, double zPos);

extern int DXFWriteText(ostream& outData, unsigned long long& handle, const string& layer,
	unsigned short color, double PointMidX, double PointMidY, double PointMidZ,
	double height, const string& content);

extern int DXFWriteArrow(ostream& outData, unsigned long long& handle,
	unsigned short color, double PointMidX, double PointMidY,
	double PointMidZ, double dx, double dy, double dz);

int DXFWrite(ostream& outData,
	vector <_tube>& Tubes, size_t mTb,
	vector <_point>& Points, size_t mPt) {
	string LayerName, content;
	int error = 0;
	size_t iTb, iPt;
	_tube* iTube;
	//center point of line
	double PointMidX, PointMidY, PointMidZ, dx, dy, dz, OCSMidAngle;
	_point* PtIn;
	_point* PtOut;
	// font height
	double height;
	//using new handles 
	unsigned long long handle = 100uLL;
	// set color
	unsigned short color = 0;
	//determine max and min extension 
	double zExtMax = -200000.,
		zExtMin = 200000.;
	for (auto& element : Points) {
		if (element.zCoord > zExtMax) zExtMax = element.zCoord;
		if (element.zCoord < zExtMin) zExtMin = element.zCoord;
	}
	height = (zExtMax - zExtMin) / 100.;

	outData << "999\nDXF2WC\n  0\nSECTION\n  2\nCLASSES\n  0\nENDSEC\n  0\nSECTION\n  2\nENTITIES\n";

	for (iTb = 0; iTb <= mTb; iTb++) {
		LayerName = "Tb " + to_string(iTb);
		iTube = &Tubes[iTb];
		PtIn = &Points[iTube->PointIn];
		PtOut = &Points[iTube->PointOut];
		if (iTube->RadiusBend < 1.) { //straight tube
			error = DXFWriteLine(outData, handle, LayerName, color,
				PtIn->xCoord, PtIn->yCoord, PtIn->zCoord,
				PtOut->xCoord, PtOut->yCoord, PtOut->zCoord);
			PointMidX = (PtIn->xCoord + PtOut->xCoord) / 2.;
			PointMidY = (PtIn->yCoord + PtOut->yCoord) / 2.;
			PointMidZ = (PtIn->zCoord + PtOut->zCoord) / 2.;
		}
		else { //bend
			error = DXFWriteArc(outData, handle, LayerName, color,
				iTube->OCSCenterX, iTube->OCSCenterY, iTube->OCSCenterZ,
				iTube->RadiusBend, iTube->Nx, iTube->Ny, iTube->Nz,
				iTube->OCSStartAngle, iTube->OCSEndAngle);
			if (iTube->OCSEndAngle - iTube->OCSStartAngle < 0.) {
				OCSMidAngle = (iTube->OCSStartAngle + iTube->OCSEndAngle + 360.) / 2.;
				if (OCSMidAngle > 360.) OCSMidAngle -= 360.;
			}
			else {
				OCSMidAngle = (iTube->OCSStartAngle + iTube->OCSEndAngle) / 2.;
			}
			_point WCSMidPoint = OCS2WCS(iTube->OCSCenterX, iTube->OCSCenterY, iTube->OCSCenterZ,
				iTube->Nx, iTube->Ny, iTube->Nz, iTube->RadiusBend, OCSMidAngle);
			PointMidX = WCSMidPoint.xCoord;
			PointMidY = WCSMidPoint.yCoord;
			PointMidZ = WCSMidPoint.zCoord;
		}
		dx = PtOut->xCoord - PtIn->xCoord;
		dy = PtOut->yCoord - PtIn->yCoord;
		dz = PtOut->zCoord - PtIn->zCoord;

		error = DXFWriteArrow(outData, handle, color,
			PointMidX, PointMidY, PointMidZ, dx, dy, dz);
		error = DXFWriteText(outData, handle, "Line number", color,
			PointMidX, PointMidY, PointMidZ, height, content);
	}
	for (iPt = 0; iPt <= mPt; iPt++) {
		LayerName = "Pt " + to_string(iPt);
		error = DXFWritePoint(outData, handle, LayerName, color,
			Points[iPt].xCoord,
			Points[iPt].yCoord,
			Points[iPt].zCoord);
	}

	outData << "  0\nENDSEC\n  0\nSECTION\n  2\nOBJECTS\n  0\nENDSEC\n  0\nEOF\n";

	return error;
}
