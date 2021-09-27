/**
 *\file main.cpp
 *\brief This program is part of the Boiler Water Circulation suite
 *
 * it reads a .dxf file, checks for correct input and write data to .dat file to be read by circulation program
 * \author rainer_jordan-at-<very, very warm>mail.com
 * #### Licence
 * Licensed under the European Union Public Licence (EUPL), Version 1.2 or later* \author rainer
 * \date   February 2021
 */
 /* Created on November 23, 2016, 7:09 PM  */

#include "common.h"

size_t mTb; ///< max index of tubes (max number of tubes - 1)
size_t mPt; ///< max index of Points (number of points - 1)  

using namespace std;

extern int DXFWrite(ostream& outData,
	vector <_tube>& Tubes, size_t mTb,
	vector <_point>& Points, size_t mPt);
extern int DXFWriteError(ostream& outData,
	vector <size_t>& OrphanedPoints,
	vector <_point>& Points,
	vector<_point>& TubeErrors,
	vector<_point>& Intersect);
extern int check_intersection(vector <_tube>& Tubes, vector <_point>& Points, vector<_point>& Intersect);

/**
 *\fn OCS2WCS(double OCSCenterX, double OCSCenterY, double OCSCenterZ,
	double Nx, double Ny, double Nz,
	double OCSRadius, double OCSAngle)
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
 * \return point, in x,y,z coordinates
 */
_point OCS2WCS(double OCSCenterX, double OCSCenterY, double OCSCenterZ,
	double Nx, double Ny, double Nz,
	double OCSRadius, double OCSAngle) {

	/// Autocad (TM) help gives a limit of 1/64 for Nx and Ny to switch to different calculation
	/// otherwise division by 0 possible
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


/**
 * \fn ProcessError (int NoLineDxf,  double p1x, double p2x, double p1y, double p2y, double p1z, double p2z,
 *	vector<_point>& ErrorsVector, bool& dxfError)
 * \brief writes an error message and stores center point in Errors vector
 *
 * \param NoLineDxf line number in .dxf file where the error occurred
 * \param p1x x-coordinate of line start point
 * \param p2x x-coordinate of line end point
 * \param p1y y-coordinate of line start point
 * \param p2y y-coordinate of line end point
 * \param p1z z-coordinate of line start point
 * \param p2z z-coordinate of line end point
 * \param ErrorsVector vector of points pointing to a faulty line (tube)
 * \param dxfError an error is detected, setting this variable to true
 */
void ProcessError(int NoLineDxf, double p1x, double p2x, double p1y, double p2y, double p1z, double p2z,
	std::vector<_point>& ErrorsVector, bool& dxfError)
{
	_point center;
	cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file";
	cout << "\n Please correct dxf-file" << endl;
	center.xCoord = (p1x + p2x) / 2.;
	center.yCoord = (p1y + p2y) / 2.;
	center.zCoord = (p1z + p2z) / 2.;
	ErrorsVector.push_back(center);
	dxfError = true;
}

/**
 * \fn AddPointIfNotExist (size_t& mPt, double px, double py, double pz, vector<_point>& Points)
 * \brief checks if a point with given coordinates exist, if not add it to Points vector
 *
 * \param px x-coordinate of point
 * \param py y-coordinate of point
 * \param pz z-coordinate of point
 * \param Points vector of points
 * \return index of point in Points vector
 */
size_t AddPointIfNotExist(double px, double py, double pz, vector<_point>& Points) {
	for (size_t iPt = 0; iPt <= mPt; iPt++) {
		if (fabs(Points[iPt].xCoord - px) < 1. &&
			fabs(Points[iPt].yCoord - py) < 1. &&
			fabs(Points[iPt].zCoord - pz) < 1.) {
			return iPt;
		}
	}
	// not found in Point vector -> add
	Points.push_back(_point());
	mPt++;
	Points[mPt].xCoord = px;
	Points[mPt].yCoord = py;
	Points[mPt].zCoord = pz;
	return mPt;
}

/**
 * \fn isPointOnDrumShell(size_t CenterPoint, size_t iPoint, vector<_point> Points, double Diameter, double Thickness)
 * \brief Checks whether given point is on drum shell
 *
 * \param CenterPoint center point of drum (steam drum or mud drum)
 * \param iPoint index of point in Points vector
 * \param Points vector of points
 * \param Diameter outside diameter of drum
 * \param Thickness thickness of drum shell
 * \return true if point is on drum shell, false if not
 */
bool isPointOnDrumShell(size_t CenterPoint, size_t iPoint, vector<_point> Points, double Diameter, double Thickness) {
	double dist, diffx, diffy, diffz;
	// checking in x-z plane 
	diffx = Points[CenterPoint].xCoord - Points[iPoint].xCoord;
	diffz = Points[CenterPoint].zCoord - Points[iPoint].zCoord;
	dist = sqrt(diffx * diffx + diffz * diffz);

	if (dist <= Diameter / 2. + 1. && dist >= Diameter / 2. - Thickness - 1.) return true;
	// checking in y-z plane 
	diffy = Points[CenterPoint].yCoord - Points[iPoint].yCoord;
	diffz = Points[CenterPoint].zCoord - Points[iPoint].zCoord;
	dist = sqrt(diffy * diffy + diffz * diffz);

	if (dist <= Diameter / 2. + 1. && dist >= Diameter / 2. - Thickness - 1.) return true;
	return false;
}

/**
 * \fn makeDrumTube(size_t CenterPoint, size_t iPoint, size_t& mTb,
 *	vector<_tube>& Tubes, double Diameter, double Thickness)
 * \brief generates connection tubes between drum center and points on shell
 *
 * \param CenterPoint point number of drum (per default steam drum = DRUM=0)
 * \param iPoint point on drum shell (start or end point of a tube)
 * \param Tubes vector of tubes
 * \param Diameter diameter of drum
 * \param Thickness shell thickness of drum
 */
void makeDrumTube(size_t CenterPoint, size_t iPoint, 
	vector<_tube>& Tubes, double Diameter, double Thickness) {
	mTb++;
	Tubes.push_back(_tube());
	Tubes[mTb].PointIn = iPoint;
	Tubes[mTb].PointOut = CenterPoint;
	Tubes[mTb].dOut = Diameter;
	Tubes[mTb].Thickness = Thickness;
	Tubes[mTb].NoParallel = 1.;
	if (CenterPoint == DRUM) {
		Tubes[mTb].Name = "Steam Drum";
	}
	else {
		Tubes[mTb].Name = "Mud Drum";
	}

	return;
}
/**
 * \fn ReadCheckCoord(ifstream& inData, int NoLineDxf, size_t iTb, string error)
 * \brief reads a text line from .dxf file and converts it to double value, if something wrong stop of program
 *
 * \param inData input stream (.dxf file)
 * \param NoLineDxf line number in .dxf file
 * \param iTb number of tube (index in Tubes vector)
 * \param error string indicating what should be read
 * \return value of read .dxf text line
 */
double ReadCheckCoord(ifstream& inData, int NoLineDxf, size_t iTb, string error) {
	string row;
	string::size_type tailptr;
	double data;
	if (getline(inData, row).good()) {
		try {
			data = stod(row, &tailptr);
		}
		catch (const std::invalid_argument& ia) {
			cout << "Invalid " << error << "point of line (tube) " << iTb;
			cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file";
			cout << "\nPlease correct dxf-file" << endl;
			exit(1);
		}
	}
	else {
		cout << "Invalid " << error << "point of line (tube) " << iTb;
		cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
		cout << "\nPlease correct dxf-file" << endl;
		exit(1);
	}
	return data;
}

/**
 * \fn ReadCheckLayerInfo(istringstream& f, const string& LayerInf, int NoLineDxf, size_t iTb,
 *	double p1x, double p1y, double p1z,	double p2x, double p2y, double p2z,
 *	const string& ErrorString, double& data, bool isCritical, bool& dxfError, vector <_point>& TubeErrors)
 * \brief  read and checks one part of the layer info at a time
 *
 * \param f subset of LayerInf that is not processed yet
 * \param LayerInf string containing the layer information, only needed for error message
 * \param NoLineDxf number of line in .dxf file, only needed for error message
 * \param iTb tube number, only needed for error message
 * \param p1x x-coordinate of tube inlet point
 * \param p1y y-coordinate of tube inlet point
 * \param p1z z-coordinate of tube inlet point
 * \param p2x x-coordinate of tube outlet point
 * \param p2y y-coordinate of tube outlet point
 * \param p2z z-coordinate of tube outlet point
 * \param ErrorString needed for error message
 * \param data the number to be extracted from string
 * \param isCritical this number is needed further on, if not available causes an error
 * \param dxfError indicator if an error occurred
 * \param TubeErrors vector of points where an error occurred
 * \return
 */
bool ReadCheckLayerInfo(istringstream& f, const string& LayerInf, int NoLineDxf, size_t iTb,
	double p1x, double p1y, double p1z, double p2x, double p2y, double p2z,
	const string& ErrorString, double& data, bool isCritical, bool& dxfError, vector <_point>& TubeErrors) {
	//
	string token = "";
	string::size_type tailptr;
	_point center;
	// reading the next token from istringstream, separator is _ 
	getline(f, token, '_');
	if (token.empty()) {
		// if nothing read but a value needed -> error
		if (isCritical) {
			cout << "\n No " << ErrorString << " in Layer of tube # " << iTb << " " << LayerInf;
			ProcessError(NoLineDxf - 23, p1x, p2x, p1y, p2y, p1z, p2z, TubeErrors, dxfError);
		}
		return false;
	}
	else {
		// if something read but can not be converted to a number -> error
		try {
			data = stod(token, &tailptr);
		}
		catch (const std::invalid_argument& ia) {
			cout << "Invalid " << ErrorString << ": " << LayerInf;
			ProcessError(NoLineDxf - 23, p1x, p2x, p1y, p2y, p1z, p2z, TubeErrors, dxfError);
			return false;
		}
	}
	return true;
}
/**
 * \fn isTubeCorrect(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z,
 *	vector<_point>& Points, vector<_tube>& Tubes, 
 *	int NoLineDxf, bool& dxfError, vector<_point>& TubeErrors, string LayerInf)
 * \brief checks if tube is correct 
 *
 * \param p1x x-coordinate of tube inlet point
 * \param p1y y-coordinate of tube inlet point
 * \param p1z z-coordinate of tube inlet point
 * \param p2x x-coordinate of tube outlet point
 * \param p2y y-coordinate of tube outlet point
 * \param p2z z-coordinate of tube outlet point
 * \param &Points vector of point(s)
 * \param &Tubes vector of tube(s)
 * \param NoLineDxf number of line in .dxf file, only needed for error message
 * \param &dxfError indicator if an error occurred
 * \param &TubeErrors vector of points where an error occurred
 * \param LayerInf string containing the layer information, only needed for error message
 * \return
 */
bool isTubeCorrect(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z,
	vector<_point>& Points,  vector<_tube>& Tubes, 
	int NoLineDxf, bool& dxfError, vector<_point>& TubeErrors, string LayerInf)
{
	double d__1, d__2, d__3, length,
		qperm; // heat absorption of all parallel tubes per m in kW/m
	char confirm;
	_point center;
	string token;
	// --------------------------------------------
	// check for 0 length tube (line) and ignore it
	// --------------------------------------------
	if (fabs(p2x - p1x) < 1. &&
		fabs(p2y - p1y) < 1. &&
		fabs(p2z - p1z) < 1.) {
		cout << "zero length line (tube) found and ignored";
		cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
		do {
			cout << "please confirm [y] :";
			cin >> confirm;
		} while (confirm != 'y');
		mTb--;
		Tubes.pop_back();
		return false; //read next line
	}
	// -------------------------------------------            
	//check if this line already exists (2 lines on top of each other)
	// -------------------------------------------
	for (size_t iTb = 0; iTb < mTb; iTb++) {
		size_t iPtIn = Tubes[iTb].PointIn;
		size_t iPtOut = Tubes[iTb].PointOut;
		if ((fabs(Points[iPtIn].xCoord - p1x) < 1. &&
			fabs(Points[iPtIn].yCoord - p1y) < 1. &&
			fabs(Points[iPtIn].zCoord - p1z) < 1. &&
			fabs(Points[iPtOut].xCoord - p2x) < 1. &&
			fabs(Points[iPtOut].yCoord - p2y) < 1. &&
			fabs(Points[iPtOut].zCoord - p2z) < 1.) ||
			(fabs(Points[iPtIn].xCoord - p2x) < 1. &&
				fabs(Points[iPtIn].yCoord - p2y) < 1. &&
				fabs(Points[iPtIn].zCoord - p2z) < 1. &&
				fabs(Points[iPtOut].xCoord - p1x) < 1. &&
				fabs(Points[iPtOut].yCoord - p1y) < 1. &&
				fabs(Points[iPtOut].zCoord - p1z) < 1.)) {
			if (LayerInf == Tubes[iTb].LayerInfo) {
				//if lines (Tubes) are identical. ignore it but ask for confirmation
				// just to remind that the design is not 100% clean
				cout << "\nidentical lines (tubes) found and ignored" << endl;
				do {
					cout << "please confirm (Yes/Show in error file)[y/s] :";
					cin >> confirm;
					confirm = confirm & '\xdf';
				} while (confirm != 'Y' && confirm != 'S');
				mTb--;
				Tubes.pop_back();
				if (confirm == 'S') {
					dxfError = true;
					center.xCoord = (p1x + p2x) / 2.;
					center.yCoord = (p1y + p2y) / 2.;
					center.zCoord = (p1z + p2z) / 2.;
					TubeErrors.push_back(center);
				}
				break;
			}
			else {
				/* lines have same coordinates but different layer information -> 2 different tubes on top of each other
				 dxf file need to be corrected*/
				cout << "\nLines with coordinates\n " <<
					p1x << " " << p1y << " " << p1z << " " <<
					p2x << " " << p2y << " " << p2z <<
					"\n and layername \n " <<
					Tubes[iTb].LayerInfo << "\n "
					<< LayerInf <<
					"\nlay on top of each other" << endl;
				ProcessError(NoLineDxf, p1x, p2x, p1y, p2y, p1z, p2z, TubeErrors, dxfError);
			}
			return false;
		}
	}
	// ----------------------------------     
	// end of checking for faulty lines      
	// ----------------------------------

	//set number of inlet and outlet point, check if they already exist
	Tubes[mTb].PointIn = AddPointIfNotExist(p1x, p1y, p1z, Points);
	Tubes[mTb].PointOut = AddPointIfNotExist(p2x, p2y, p2z, Points);

	//parse of LayerInfo
	//the layername contains all information 
	//separator is underscore
	// name, diameter and thickness are needed to define the tube
	// if something wrong save center point in TubeError vector
	// if the other values are not given, standard values are used 
	istringstream f(LayerInf);
	if (!getline(f, token, '_').good()) {
		cout << "\nNo tube name in Layer of tube # " << mTb << " " << LayerInf;
		ProcessError(NoLineDxf - 23, p1x, p2x, p1y, p2y, p1z, p2z, TubeErrors, dxfError);
		return false;
	}
	else {
		Tubes[mTb].Name = token;
	}

	ReadCheckLayerInfo(f, LayerInf, NoLineDxf, mTb, p1x, p1y, p1z, p2x, p2y, p2z,
		"tube diameter", Tubes[mTb].dOut, true, dxfError, TubeErrors);
	ReadCheckLayerInfo(f, LayerInf, NoLineDxf, mTb, p1x, p1y, p1z, p2x, p2y, p2z,
		"tube thickness", Tubes[mTb].Thickness, true, dxfError, TubeErrors);
	qperm = 0.;
	if (ReadCheckLayerInfo(f, LayerInf, NoLineDxf, mTb, p1x, p1y, p1z, p2x, p2y, p2z,
		"heat per m tube length", qperm, false, dxfError, TubeErrors)) {
		if (ReadCheckLayerInfo(f, LayerInf, NoLineDxf, mTb, p1x, p1y, p1z, p2x, p2y, p2z,
			"number of parallel tubes", Tubes[mTb].NoParallel, false, dxfError, TubeErrors)) {
			if (ReadCheckLayerInfo(f, LayerInf, NoLineDxf, mTb, p1x, p1y, p1z, p2x, p2y, p2z,
				"additional resistance factor", Tubes[mTb].ksiAdd, false, dxfError, TubeErrors)) {
				if (ReadCheckLayerInfo(f, LayerInf, NoLineDxf, mTb, p1x, p1y, p1z, p2x, p2y, p2z,
					"orifice diameter at tube inlet", Tubes[mTb].DiaOrificeIn, false, dxfError, TubeErrors)) {
					if (ReadCheckLayerInfo(f, LayerInf, NoLineDxf, mTb, p1x, p1y, p1z, p2x, p2y, p2z,
						"orifice diameter tube at outlet", Tubes[mTb].DiaOrificeOut, false, dxfError, TubeErrors)) {
						if (ReadCheckLayerInfo(f, LayerInf, NoLineDxf, mTb, p1x, p1y, p1z, p2x, p2y, p2z,
							"heating factor", Tubes[mTb].FactHeat, false, dxfError, TubeErrors)) {
							ReadCheckLayerInfo(f, LayerInf, NoLineDxf, mTb, p1x, p1y, p1z, p2x, p2y, p2z,
								"given inlet enthalpy", Tubes[mTb].EnthGiven, false, dxfError, TubeErrors);
						}
					}
				}
			}
		}
	}
	if (Tubes[mTb].RadiusBend < 1.) {
		d__1 = p2x - p1x;
		d__2 = p2y - p1y;
		d__3 = p2z - p1z;
		length = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	}
	else {
		length = fabs(Tubes[mTb].OCSEndAngle - Tubes[mTb].OCSStartAngle) / 180. * M_PI * Tubes[mTb].RadiusBend ;
	}
	Tubes[mTb].q = qperm * length / 1e3;
	// further checking inside diameter, short length     
	if (length < 10.) { //should be longer 10 mm
		cout << "Invalid tube length: ";
		ProcessError(NoLineDxf - 23, p1x, p2x, p1y, p2y, p1z, p2z, TubeErrors, dxfError);
		return false;
	}
	if ((Tubes[mTb].dOut - 2. * Tubes[mTb].Thickness) < 5.) {
		cout << "Invalid tube inside diameter: " << LayerInf;
		ProcessError(NoLineDxf - 23, p1x, p2x, p1y, p2y, p1z, p2z, TubeErrors, dxfError);
		return false;
	}
	return true;
}

/**
 *
 * @param argc
 * @param argv can contain the project name, if empty user will be asked for file name
 * @return
 */
int main(int argc, char** argv) {
	/* Local variables */
	bool dxfError = false, entities = false;
	char confirm, Method;
	int maxit, NoLineDxf = 0, control, error;
	size_t iTb, iPt;
	string filename, readfile, writefile, row, LayerInf, token, Project;
	double p1x =0., p1y=0., p1z=0., p2x=0., p2y=0., p2z=0.;
	double tolerance;
	double pressure;
	double roughness;
	double LevelFact;
	double LevelWDrum;
	double dpDynDrum;
	double CRStart;
	vector <_tube> Tubes;
	vector <_point> Points;
	vector <size_t> OrphanedPoints;
	vector <_point> TubeErrors;
	vector <_point> Intersect;
	_point center;
	ifstream inData;
	ofstream outData;
	ofstream prot;

	/**
	 * -# defining point for drum first\n
	 * 	to clearly define the drum the center has to sit at 100 000, 100 000, 100 000 mm
	 * -# reading project name if not already passed as argument and determine the name of input and output files
	 * -# open the input file (.dxf file)
	 * -# read all lines until content is "EOF" (End Of File)
	 *    - there can be lines in some definition like blocks\n
	 *    those should not be processed\n
	 *    only the lines in ENTITIES section should be taken
	 * -# read data of new tube. If error: indicate at which line (approximately)
	 * -# data for this line is read, now it must be handled
	 *    - check for 0 length tube (line) and ignore it
	 *    - check if this line already exists (2 lines on top of each other)\n
	 *		if lines (Tubes) are identical. ignore it but ask for confirmation\n
	 *		just to remind that the design is not 100% clean
	 *		- lines have same coordinates but different layer information -> 2 different tubes on top of each other\n
	 *		dxf file need to be corrected
	 * -# check if p1 or p2 already exist\n
	 *		if not found in Point vector -> add
	 * -# parse of LayerInfo\n
	 *		the layername contains all information\n
	 *		separator is underscore\n
	 *		name, diameter and thickness are needed to define the tube. If something wrong save center point in TubeError vector\n
	 *		if the other values are not given, standard values are used
	 * -# check: each point need to have at least 2 tubes connected to this point\n
	 *		 if one point only has 1 line (tube) connected, save this point in vector OrphanedPoint\n
	 *     check for intersection of lines that are not start or end point (2 lines crossing)
	 * -# to allow different ways to produce or edit circulation calculation input file,the check for branches, arriving/leaving branches at node and flow direction should be done in water circulation program
	 * -# if there is something wrong like missing data in layer info, orphaned point or 2 lines(tubes) crossing write the coordinates to error dxf-file\n
	 *		 this dxf-file can be overlaid to the original faulty drawing. The lines in this file point to the erroneous positions
	 * -# if geometry is ok, input base data like drum pressure, calculation method etc\n
	 *		and save to circulation data file
	 * -# write a new dxf-file with tube numbers and point numbers as layer names\n
	 * 		the numbering of points and tubes was done in this program. The dxf file helps to identfy the position of the tubes and points
	 */

	 /*defining point for drum first
	  * to clearly define the drum the center has to sit at 100 000, 100 000, 100 000 mm
	  */
	mPt = 0;
	Points.push_back(_point());
	Points[DRUM].xCoord = 100000.;
	Points[DRUM].yCoord = 100000.;
	Points[DRUM].zCoord = 100000.;

	/*
	 * reading project name if not already passed as argument and
	 *  determine the name of input and output file
	 */
	if (argc > 1) {
		Project = argv[1];
	}
	else {
		cout << "\nProject name : ";
		cin >> Project;
	}
//		filename = "..\\WSCDATA\\" + Project + "\\" + Project;
 //    filename = "..\\..\\WSCDATA\\" + Project + "\\" + Project;
	filename = "WSCDATA\\" + Project + "\\" + Project;
	readfile = filename + ".dxf";
	writefile = filename + ".dat";
	/*
	 * open the input file (.dxf file)
	 */
	inData.open(readfile.c_str());
	if (!inData.good()) {
		cout << "\n cannot open file " << readfile << endl;
		exit(1);
	}
	else {
		mTb = MINUS1;
		while (true) {
			if (getline(inData, row).good()) {
				NoLineDxf++;
				/*
				 * read all lines until content is "EOF" (End Of File)
				 */
				if ((entities && row == "ENDSEC") || row == "EOF") {
					break;
				}
				else if (row == "ENTITIES") {
					entities = true;
					/*
					 * there can be lines in some definition like blocks
					 * those should not be processed
					 * only the lines in ENTITIES section should be taken
					 */
				}
				else if (entities && row == "LINE") {
					//new tube 
					mTb++;
					Tubes.push_back(_tube());
					while (true) {
						NoLineDxf++;
						if (getline(inData, row).good()) {
							if (row == "  0") {
								break;
							}
							else if (row == "  8") {
								NoLineDxf++;
								if (getline(inData, LayerInf).good()) {
									Tubes[mTb].LayerInfo = LayerInf;
								}
								else {
									cout << "Error reading layer name of line (tube) " << mTb;
									cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
								}
							}
							else if (row == " 10") {
								NoLineDxf++;
								p1x = ReadCheckCoord(inData, NoLineDxf, mTb, "x-coordinate of starting point");
							}
							else if (row == " 20") {
								NoLineDxf++;
								p1y = ReadCheckCoord(inData, NoLineDxf, mTb, "y-coordinate of starting point");
							}
							else if (row == " 30") {
								NoLineDxf++;
								p1z = ReadCheckCoord(inData, NoLineDxf, mTb, "z-coordinate of starting point");
							}
							else if (row == " 11") {
								NoLineDxf++;
								p2x = ReadCheckCoord(inData, NoLineDxf, mTb, "x-coordinate of end point");
							}
							else if (row == " 21") {
								NoLineDxf++;
								p2y = ReadCheckCoord(inData, NoLineDxf, mTb, "y-coordinate of end point");
							}
							else if (row == " 31") {
								NoLineDxf++;
								p2z = ReadCheckCoord(inData, NoLineDxf, mTb, "z-coordinate of end point");
							}
						}
						else {
							cout << "Error reading line (tube) " << mTb;
							cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
							cout << "\nPlease correct dxf-file" << endl;
							exit(1);
						}
					}

					//data for this line (tube) is read, now it must be handled
					isTubeCorrect(p1x, p1y, p1z, p2x, p2y, p2z, Points,  Tubes, 
						NoLineDxf, dxfError, TubeErrors, LayerInf);
				}
				else if (entities && row == "ARC") {
					//new arc
					mTb++;
					Tubes.push_back(_tube());
					while (true) {
						NoLineDxf++;
						if (getline(inData, row).good()) {
							if (row == "  0") {
								break;
							}
							else if (row == "  8") {
								NoLineDxf++;
								if (getline(inData, LayerInf).good()) {
									Tubes[mTb].LayerInfo = LayerInf;
								}
								else {
									cout << "Error reading layer name of line (tube) " << mTb;
									cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
								}
							}
							else if (row == " 10") {
								NoLineDxf++;
								Tubes[mTb].OCSCenterX = ReadCheckCoord(inData, NoLineDxf, mTb, "x-coordinate of OCS center point");
							}
							else if (row == " 20") {
								NoLineDxf++;
								Tubes[mTb].OCSCenterY = ReadCheckCoord(inData, NoLineDxf, mTb, "y-coordinate of OCS center point");
							}
							else if (row == " 30") {
								NoLineDxf++;
								Tubes[mTb].OCSCenterZ = ReadCheckCoord(inData, NoLineDxf, mTb, "z-coordinate of OCS center point");
							}
							else if (row == " 40") {
								NoLineDxf++;
								Tubes[mTb].RadiusBend = ReadCheckCoord(inData, NoLineDxf, mTb, "Radius of arc");
							}
							else if (row == " 50") {
								NoLineDxf++;
								Tubes[mTb].OCSStartAngle = ReadCheckCoord(inData, NoLineDxf, mTb, "Start angle of arc");
							}
							else if (row == " 51") {
								NoLineDxf++;
								Tubes[mTb].OCSEndAngle = ReadCheckCoord(inData, NoLineDxf, mTb, "End angle of arc");
							}
							else if (row == "210") {
								NoLineDxf++;
								Tubes[mTb].Nx = ReadCheckCoord(inData, NoLineDxf, mTb, "x-component of extrusion vector");
							}
							else if (row == "220") {
								NoLineDxf++;
								Tubes[mTb].Ny = ReadCheckCoord(inData, NoLineDxf, mTb, "y-component of extrusion vector");
							}
							else if (row == "230") {
								NoLineDxf++;
								Tubes[mTb].Nz = ReadCheckCoord(inData, NoLineDxf, mTb, "z-component of extrusion vector");
							}
						}
						else {
							cout << "\nError reading arc " << mTb;
							cout << "\nThe error occurred at about line " << NoLineDxf << " of dxf-file" << endl;
							cout << "\nPlease correct dxf-file" << endl;
							exit(1);
						}
					}

					_point WCSPointIn = OCS2WCS(Tubes[mTb].OCSCenterX, Tubes[mTb].OCSCenterY, Tubes[mTb].OCSCenterZ,
						Tubes[mTb].Nx, Tubes[mTb].Ny, Tubes[mTb].Nz, Tubes[mTb].RadiusBend, Tubes[mTb].OCSStartAngle);
					_point WCSPointOut = OCS2WCS(Tubes[mTb].OCSCenterX, Tubes[mTb].OCSCenterY, Tubes[mTb].OCSCenterZ,
						Tubes[mTb].Nx, Tubes[mTb].Ny, Tubes[mTb].Nz, Tubes[mTb].RadiusBend, Tubes[mTb].OCSEndAngle);
					//data for this line (tube) is read, now it must be handled
					isTubeCorrect(WCSPointIn.xCoord, WCSPointIn.yCoord, WCSPointIn.zCoord,
						WCSPointOut.xCoord, WCSPointOut.yCoord, WCSPointOut.zCoord, Points, Tubes,
						NoLineDxf, dxfError, TubeErrors, LayerInf);

				} // end of this tube (line)
			} // valid line read
		} // while loop for line reading 
		inData.close();
	}//file opened correctly
	cout << "\n All " << mTb << " lines read" << endl;

	// check: each point need to have at least 2 tubes connected to this point 
	// exception: drum, it can have no tubes if it is a bidrum boiler and the pseudo tubes in drum have to be determined
	for (iPt = 1; iPt <= mPt; iPt++) {
		Points[iPt].NoTb = 0;
		for (auto& iTube : Tubes) {
			if (iPt == iTube.PointIn ||
				iPt == iTube.PointOut) {
				Points[iPt].NoTb++;
			}
		}
		// if one point only has 1 line (tube) connected, save this point in vector OrphanedPoint 
		if (Points[iPt].NoTb < 2) {
			OrphanedPoints.push_back(iPt);
		}
	}

	if (OrphanedPoints.size() > 10) {
		cout << "\n big number of points with only one tube connected";
		cout << "\n in bi-drum boilers there is a big number of pseudo tubes in the drums ";
		cout << "\n those can be generated automatically ";
		do {
			cout << "\n Is it a model for bi-drum boiler [y/n]?";
			cin >> confirm;
			confirm = confirm & '\xdf';
		} while (!(confirm == 'Y' || confirm == 'N'));
		if (confirm == 'Y') {
			double SteamDrumDiameter, SteamDrumThickness, MudDrumDiameter, MudDrumThickness,
				mdx, mdy, mdz;
			size_t mdp;
			// center point of steam drum already given; orientation of drum unknown 
			cout << "\n please give steam drum outside diameter [mm] :";
			cin >> SteamDrumDiameter;
			cout << "\n please give steam drum shell thickness [mm] :";
			cin >> SteamDrumThickness;

			cout << "\n please give x-coordinate of mud drum center point [mm] :";
			cin >> mdx;
			cout << "\n please give y-coordinate of mud drum center point [mm] :";
			cin >> mdy;
			cout << "\n please give z-coordinate of mud drum center point [mm] :";
			cin >> mdz;
			cout << "\n please give mud drum outside diameter [mm] :";
			cin >> MudDrumDiameter;
			cout << "\n please give mud drum shell thickness [mm] :";
			cin >> MudDrumThickness;
			mdp = AddPointIfNotExist(mdx, mdy, mdz, Points);

			// checking steam drum first  
			for (size_t i = 0; i < OrphanedPoints.size(); i++) {
				iPt = OrphanedPoints[i];
				if (isPointOnDrumShell(DRUM, iPt, Points, SteamDrumDiameter, SteamDrumThickness)) {
					makeDrumTube(DRUM, iPt, Tubes, SteamDrumDiameter, SteamDrumThickness);
					OrphanedPoints.erase(OrphanedPoints.begin() + i);
					i--;
					continue;
				}
				if (isPointOnDrumShell(mdp, iPt, Points, MudDrumDiameter, MudDrumThickness)) {
					makeDrumTube(mdp, iPt, Tubes, MudDrumDiameter, MudDrumThickness);
					OrphanedPoints.erase(OrphanedPoints.begin() + i);
					i--;
				}
			}
		}
	}
	if (!OrphanedPoints.empty()) {
		dxfError = true;
		for (size_t iPto : OrphanedPoints) {
			cout << "\n Point " << iPto << " with coordinates \n x " << Points[iPto].xCoord <<
				" y " << Points[iPto].yCoord << " z " << Points[iPto].zCoord << "\n has "
				<< Points[iPto].NoTb << " tube connected" <<
				"\n Please correct dxf-file" << endl;
		}
	}

	check_intersection(Tubes, Points, Intersect);

	if (!Intersect.empty()) {
		dxfError = true;
		cout << "\n Intersection of tubes beside start or end points " <<
			"\n Please see error file for position and correct dxf-file" << endl;
	}
	//to allow different sources check for branches, arriving/leaving branches at node 
	// and flow direction should be done in water circulation program            
	if (dxfError) {
		/* if there is something wrong like missing data in layer info or orphaned point
		 * write the coordinates to error dxf-file
		 * this dxf-file can be overlayed to the original faulty drawing
		 * the lines in this file point to the erroneous positions */
		writefile = filename + "_Error.dxf";

		outData.open(writefile.c_str());

		if (!outData.good()) {
			cout << "\n cannot open file " << writefile << endl;
			exit(1);
		}
		else {
			error = DXFWriteError(outData, OrphanedPoints, Points, TubeErrors, Intersect);
			outData.close();
			cout << "\n dxf file contains errors \nPlease check " << writefile << " for position of the errors ";
		}//file opened correctly

	}
	else {
		// if geometry is ok, input base data like drum pressure, calculation method etc
		cout << "\n maximum number of iterations ";
		cin >> maxit;
		cout << "\n tolerance ";
		cin >> tolerance;
		cout << "\n control parameter for protocol file :";
		cout << "\n show nothing: 0 ";
		cout << "\n show info about mesh/network : 1 ";
		cout << "\n show info about calculation of enthalpy at branches inlet : 2 ";
		cout << "\n show info about pressure drop in tubes : 3 ";
		cout << "\n show info about pressure drop in branches :5 ";
		cout << "\n show info about nodal pressure calculation : 6 ";
		cout << "\n show info and .dxf about flow for next iteration step : 7 ";
		cout << "\n show info reversal of flow direction in branches and tubes : 8 ";
		cout << "\n show all info (careful, protocol file can get HUGE) : 10  : ";
		cout << "\n show more details about mesh/network : 11 ";
		cout << "\n show more details about pressure drop in tubes : 13 ";
		cout << "\n show more info reversal of flow direction in branches and tubes : 18 ";
		cout << "\n show everything (careful, protocol file can get really HUGE) : 20 : ";
		cin >> control;
		do {
			cout << "\n calculation method for 2-phase density and pressure drop"
				<< "\n J for Jirous/Jirous"
				<< "\n R for Rouhani/Becker"
				<< "\n E for Chexal/VDIHeatAtlas "
				<< "\n W for Welzer/Welzer "
				<< "\n C for Woldesemayat/VDIHeatAtlas ";
			cin >> Method;
			Method = Method & '\xdf';
		} while (!(Method == 'J' || Method == 'R' || Method == 'E' || Method == 'W' || Method == 'C'));
		cout << "\n drum pressure [MPag] ";
		cin >> pressure;
		cout << "\n tube roughness [mm] ";
		cin >> roughness;
		cout << "\n LevelFactor between highest and lowest heated point ";
		cout << "\n downcomers have to pass through this level, 0 ... 1 ";
		cin >> LevelFact;
		cout << "\n Water level in drum, difference to drum center [mm] ";
		cin >> LevelWDrum;
		cout << "\n Resistance of drum internals [kPa] ";
		cin >> dpDynDrum;
		cout << "\n Circulation ratio for start values \n     a number between 10 and 30 would be a good initial guess ";
		cin >> CRStart;
		//  save to circulation data file 
		outData.open(writefile.c_str());
		outData << filename << endl; // project name == filename
		outData << maxit << " " << tolerance << " " << control << " " << Method << endl;
		outData << pressure + 0.1 << " " << roughness << " " << LevelFact
			<< " " << LevelWDrum << " " << dpDynDrum << " " << CRStart << endl;

		for (iPt = 0; iPt <= mPt; iPt++) {
			outData.setf(ios::fixed, ios::floatfield);
			outData << setw(7) << iPt << " "
			   << setw(12) << setprecision(10) << Points[iPt].xCoord << " "
				<< setw(12) << setprecision(10) << Points[iPt].yCoord << " "
				<< setw(12) << setprecision(10) << Points[iPt].zCoord << "\n";
      }
		outData << "TubeData\n";
		for (iTb = 0; iTb <= mTb; iTb++) {
			outData << Tubes[iTb].Name << endl;
			outData.setf(ios::fixed, ios::floatfield);
			outData << setw(7) << iTb << setw(7) << Tubes[iTb].PointIn << " "
				<< setw(7) << Tubes[iTb].PointOut << " "
				<< setw(8) << setprecision(2) << Tubes[iTb].dOut << " "
				<< setw(6) << setprecision(2) << Tubes[iTb].Thickness << " "
				<< setw(8) << setprecision(2) << Tubes[iTb].RadiusBend << " "
				<< setw(7) << setprecision(1) << Tubes[iTb].NoParallel << " "
				<< setw(8) << setprecision(2) << Tubes[iTb].ksiAdd << " "
				<< setw(8) << setprecision(2) << Tubes[iTb].q << " "
				<< setw(5) << setprecision(2) << Tubes[iTb].FactHeat << " "
				<< setw(6) << setprecision(2) << Tubes[iTb].DiaOrificeIn << " "
				<< setw(6) << setprecision(2) << Tubes[iTb].DiaOrificeOut << " "
				<< setw(6) << setprecision(1) << Tubes[iTb].EnthGiven << " "
			<< setw(8) << setprecision(6) << Tubes[iTb].OCSStartAngle << " "
			<< setw(8) << setprecision(6) << Tubes[iTb].OCSEndAngle << " "
			<< setw(18) << setprecision(15) << Tubes[iTb].OCSCenterX << " "
			<< setw(18) << setprecision(15) << Tubes[iTb].OCSCenterY << " "
			<< setw(18) << setprecision(15) << Tubes[iTb].OCSCenterZ << " "
			<< setw(9) << setprecision(6) << Tubes[iTb].Nx << " "
			<< setw(9) << setprecision(6) << Tubes[iTb].Ny << " "
			<< setw(9) << setprecision(6) << Tubes[iTb].Nz << "\n";
		}
		outData.close();

		//writing tube numbers, point numbers as layer names in dxf-file
		NoLineDxf = 0;
		writefile = filename + "_TubeNo.dxf";

		outData.open(writefile.c_str());

		if (!outData.good()) {
			cout << "\n cannot open file " << writefile << endl;
			exit(1);
		}
		else {
			error = DXFWrite(outData, Tubes, mTb, Points, mPt);
			outData.close();

		}//file opened correctly
	}
	return (EXIT_SUCCESS);
}

