/*****************************************************************//**
 * \file   ReadData.cpp
 * \brief reads data from input file
*
* \param PathFile file name including path
* \return int
 *
 * \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*********************************************************************/
//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

using namespace std;

void ErrorRead(ifstream& inData, size_t iTb, int NbLine, const string& text ){
	cout << "\n tube # " << iTb << text;
	cout << "\n error found in data file line " <<NbLine;
	prot << "\n tube # " << iTb << text;
	prot << "\n error found in data file line " << NbLine;
	inData.close();
	prot.close();
	exit(1);

   return;
}

void ErrorParsing(ifstream& inData, int NbLine) {
	cout << "\n error found in data file line " <<NbLine  << " (line parsing)";
	prot << "\n error found in data file line " <<NbLine  << " (line parsing)";
	inData.close();
	prot.close();
	exit(1);
	return;
}


int readData() {

	/* Local variables */
	string Project, iname, line, TubeName;
	int InDataLine = 3;
	size_t p1, p2, iTb, NbPt;
	int control = 0;
	double dx, dy, dz;
	double q1, x, y, z, npar1;
	double dOut, BendRadius, fAdd, fHeat, dOrificeIn, dOrificeOut, EnthGiven;
	double StartAngle, EndAngle, CenterX, CenterY, CenterZ, ArcNx, ArcNy, ArcNz;
	double thk, tSatDrum;
	double LevelHighestPoint = 0.;
	double LevelLowestPoint = 100.;
	double LevelFact = 0.;

	ifstream inData;
	// Set exceptions to be thrown on failure
	//inData.exceptions(std::ifstream::failbit | std::ifstream::badbit);

	/*     ---------------------------- */
	/*     initialize data */
	/*     ---------------------------- */
	Drum.LevelW = 0.;
	Drum.dpDyn = 0.;

	iname = PathFile + ".dat";

	try {
		inData.open(iname.c_str());
	}
	catch (system_error& e) {
		cout << e.code().message() << endl;
		exit(1);
	}
	if (getline(inData, line).good()) {
		istringstream iss(line);
		if (!(iss >> Project)) ErrorParsing(inData, 1); 
	}
	else ErrorParsing(inData, 1);
	/**the first 2 lines in input data file contain some general data
	*
	* first line
	*
	*
	* Base.maxit: maximum number of flow iterations\n
	* Base.tol: tolerance for final flow results\n
	* control: number, which data should be written to protocol file (for debugging)\n
	* Base.Method: character to indicate calculation method for 2-phase flow
	 */
	if (getline(inData, line).good()) {
		istringstream iss(line);
		if (!(iss >> Base.maxit >> Base.tol >> control >> Base.Method)) {
			ErrorParsing(inData, 2);
		} // error
	}
	else {
		ErrorParsing(inData, 2);
	}
	prot << " maxit " << Base.maxit << " control " << control;
	if (control != 0) {
		Base.showMesh = true;
		Base.showEnth = (control == 2 || control == 10 || control == 20);
		Base.showDPTube = (control == 3 || control == 10 || control == 13 || control == 20);
		Base.showDPBranch = (control == 5 || control == 10 || control == 20);
		Base.showNodePressure = (control == 6 || control == 10 || control == 20);
		Base.showFlow = (control == 7 || control == 10 || control == 20);
		Base.showReverse = (control == 8 || control == 10 || control == 18 || control == 20);
		Base.showMeshDetail = (control == 11 || control == 20);
		Base.showDPTubeDetail = (control == 13 || control == 20);
		Base.showReverseDetail = (control == 18 || control == 20);
	}
	/**
	 * second line
	 *
	 * Drum.pMPa: drum pressure [MPa]\n
	 * Base.Rough: tube roughness [mm]\n
	 * LevelFact: a factor to give a level between highest and lowest heated point, used to determine downcomers [0...1]\n
	 * Drum.LevelW: water level in drum relative to drum center, above positive, below negative [m]\n
	 * Drum.dpDyn: dynamic pressure drop in drum due to drum internals (baffles, cyclones etc.) [Pa]\n
	 * Base.CRStart: initial guess of circulation ratio for the first branch flow estimation initFlow()
	 */

	if (getline(inData, line).good()) {
		istringstream iss(line);
		if (!(iss >> Drum.pMPa >> Base.Rough >> LevelFact >> Drum.LevelW >> Drum.dpDyn >> Base.CRStart)) {
			ErrorParsing(inData, 3);
		} // error
	}
	else {
		ErrorParsing(inData, 3);
	}
	Base.Rough /= 1e3; // conversion to m
	tSatDrum = H2O::satTemp(Drum.pMPa);
	Drum.rhoW = 1. / H2O::specVol(tSatDrum, Drum.pMPa, WATER);
	Drum.enthW = H2O::enth(tSatDrum, Drum.pMPa, WATER);
	Drum.qSum = 0.;
	Drum.enthEvap = H2O::enth(tSatDrum, Drum.pMPa, STEAM) - Drum.enthW;
	Drum.dpDyn *= 1e3; // conversion to Pa
	mTb = 0;
	Tubes.push_back(_tube());
	//   cout << Tube[0].PointIn << " " << Tube[0].Points[1];
	mPt = 0;
	Points.push_back(_point());
	/**
	 * loop through data file until no more data left
	 */
	while (getline(inData, line).good()) {
		InDataLine++;
		istringstream iss(line);
/**
 * 1.) reading point data
 * 
 * NbPt -> point number\n
 * x -> x-coordinate of point[mm]\n
 * y -> y-coordinate of point[mm]\n
 * z -> z-coordinate of point[mm]\n
 * The point coordinates are kept in mm because they are needed later in .dxf file in mm\n
 */
		if (!(iss >> NbPt >> x >> y >> z)) {
			istringstream iss1(line);
			if (!(iss1 >> TubeName)) {
				ErrorParsing(inData, InDataLine);
			}
			else {
				if (TubeName == "TubeData") {
					break;
				}
			}
		}
// To avoid Points index to be different to Number, if we have gaps (1,2,3,5,...) 
		if (NbPt > mPt) {
			Points.resize(NbPt + 1);
			mPt = NbPt;
		}
		Points[NbPt].Number = NbPt;
		Points[NbPt].xCoord = x; // ->mm
		Points[NbPt].yCoord = y; // ->mm
		Points[NbPt].zCoord = z; // ->mm
		LevelHighestPoint = fmax(LevelHighestPoint, z / 1e3);// ->m
		LevelLowestPoint = fmin(LevelLowestPoint, z / 1e3);// ->m
	}
	while (getline(inData, line).good()) {
		InDataLine++;
		istringstream iss(line);
		/**
		 * 1.) reading of data
		 *
		 * first line of tube data
		 *
		 * TubeName -> Tubes[iTb].Name: TubeName;
		 *
		 */
		if (!(iss >> TubeName)) {
			ErrorParsing(inData, InDataLine);
		} // error
		if (getline(inData, line).good()) {
			InDataLine++;
			istringstream iss1(line);
			/**
			 * second line of tube data
			 *
			 * iTb Number of tube\n
			 * p1 -> Tubes[iTb].PointIn: Number of inlet point\n
			 * x1 -> Points[p1].xCoord: x-coordinate of inlet point \n
			 * y1 -> Points[p1].yCoord: y-coordinate of inlet point \n
			 * z1 -> Points[p1].zCoord: z-coordinate of inlet point  \n
			 * p2 -> Tubes[iTb].PointOut: Number of outlet point \n
			 * x2 -> Points[p2].xCoord: x-coordinate of outlet point  \n
			 * y2 -> Points[p2].yCoord: y-coordinate of outlet point  \n
			 * z2 -> Points[p2].zCoord: z-coordinate of outlet point  \n
			 * dOut: outside diameter of tube [mm]->[m]\n
			 * thk: tube thickness [mm] \n
			 * BendRadius -> Tubes[iTb].RadiusBend: bending radius of bend at tube inlet [mm]->[m] \n
			 * npar1 -> Tubes[iTb].NoParallel: number of tubes flown parallel [-]\n
			 * fAdd ->Tubes[iTb].ksiAdd: additional resistance factor [-]\n
			 * q1 ->Tubes[iTb].q: heat transferred to tube. In case of more than 1 tube flown in parallel, it is the total heat of all tubes [kW]\n
			 * fHeat factor for position of heating [-]\n
			 * dOrificeIn -> Tubes[iTb].DiaOrificeIn: diameter of orifice at tube inlet [mm]->[m]\n
			 * dOrificeOut -> Tubes[iTb].DiaOrificeOut: diameter of orifice at tube outlet [mm]->[m]\n
			 */
			if (!(iss1 >> iTb >> p1 >> p2 >> dOut >>thk >> BendRadius >> npar1 >> fAdd >> q1 >> fHeat >>
				 dOrificeIn >> dOrificeOut >> EnthGiven >>StartAngle >> EndAngle >>CenterX >>CenterY >>CenterZ >>
				ArcNx >>				ArcNy >>				ArcNz 				)) {
				ErrorParsing(inData, InDataLine);
			} // error
			while (iTb > mTb) {
				Tubes.push_back(_tube());
				mTb = iTb;
			}
			_tube* iTube = &Tubes[iTb];
				/**
			 * 2.) setting max. number of points and tubes
			 *
			 * 3.) error checking
			 *
			 * 4.) converting from mm to m
			 *
			 * 5.) writing temporary input variable into Points and Tubes vectors
			 *
			 * 6.) calculating some data like tube inside cross-section area, tube length etc. and writing to Tubes vector
			 */
			if (iTube->PointIn != MINUS1 && iTube->PointOut != MINUS1) {
				processError(iTube->PointIn, "tube already defined");
				ErrorRead(inData, iTb,InDataLine, " : this tube is already defined");
			}
			if (p1 > mPt || p2 > mPt) {
				ErrorRead(inData, iTb, InDataLine, " : tube point number higher than maximum number");
			}
			iTube->Number = iTb;
			iTube->PointIn = p1;
			iTube->PointOut = p2;

			if (dOut - thk * 2. < 1.) {
				ErrorRead(inData, iTb, InDataLine, " inside diameter too small");
			}
			iTube->Name = TubeName;
			iTube->Dia = (dOut - thk * 2.) / 1e3; // ->m
			iTube->area = M_PI_4 * (iTube->Dia * iTube->Dia);// ->m2

 			iTube->RadiusBend = BendRadius / 1e3;   /* -> m */
			_point* PtIn = &Points[iTube->PointIn];
			_point* PtOut = &Points[iTube->PointOut];
			dz = PtOut->zCoord - PtIn->zCoord;
			iTube->Height = dz / 1e3; // ->m
			if (Tubes[iTb].RadiusBend < 1e-6) {
				dx = PtOut->xCoord - PtIn->xCoord;
				dy = PtOut->yCoord - PtIn->yCoord;
				iTube->Length = sqrt(dx * dx + dy * dy + dz * dz) / 1e3;// ->m
			}
			else {
				iTube->OCSStartAngle = StartAngle; //in deg
				iTube->OCSEndAngle = EndAngle;     //in deg
				iTube->OCSCenterX = CenterX;       // ->mm
				iTube->OCSCenterY = CenterY;       // ->mm
				iTube->OCSCenterZ = CenterZ;       // ->mm
				iTube->Nx = ArcNx;
				iTube->Ny = ArcNy;
				iTube->Nz = ArcNz;
				if (fabs(EndAngle) < 1e-3) EndAngle = 360.;
				iTube->beta = EndAngle - StartAngle; //in deg
				iTube->Length = iTube->beta * M_PI / 180. * iTube->RadiusBend;// in m
			}
			if (iTube->Length < 1e-3) {
				ErrorRead(inData, iTb, InDataLine, " length too short");				
			}
			PtIn->NoTb++;
			PtIn->NbTb.push_back(iTb);
			PtOut->NoTb++;
			PtOut->NbTb.push_back(iTb);

			iTube->q = q1;
			Drum.qSum += q1;
			if (npar1 < 1.) {
				ErrorRead(inData, iTb, InDataLine, " too few number of parallel tubes");	
			}
			iTube->NoParallel = npar1;
			iTube->ksiAdd = fAdd;
			if (fHeat < .001 || fHeat > 2.) {
				ErrorRead(inData, iTb, InDataLine, " wrong heating position factor");	
			}
			iTube->FactHeat = fHeat;
			//         Tube[iTb].NoConnect = (double) nk1;
			if (dOrificeIn / 1e3 > iTube->Dia || dOrificeIn < 0.) {
				ErrorRead(inData, iTb, InDataLine," wrong inlet orifice diameter");	
			}
			iTube->DiaOrificeIn = dOrificeIn / 1e3;// ->m
			if (dOrificeOut / 1e3 > iTube->Dia || dOrificeOut < 0.) {
				ErrorRead(inData, iTb, InDataLine, " wrong outlet orifice diameter");
			}
			iTube->DiaOrificeOut = dOrificeOut / 1e3;// ->m
			iTube->EnthInGiven = EnthGiven;
			iTube->HeatFlux = q1 /
				(M_PI * iTube->Dia * iTube->Length * iTube->NoParallel);
		}
		else {
			break;
		}
	}
	inData.close();
	Base.DowncomerLevel = (LevelHighestPoint - LevelLowestPoint) * LevelFact
		+ LevelLowestPoint;
	Drum.FlowSteam = Drum.qSum / Drum.enthEvap;
	prot << "\n highest " << LevelHighestPoint << " LevelLowestPoint " << LevelLowestPoint;
	prot << "\n downcomerlevel " << Base.DowncomerLevel;
	Points.shrink_to_fit();
	Tubes.shrink_to_fit();

	return 0;
}
