//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"
/*!
 * \file Print2dxf.cpp
 * \brief opens .dxf file according to mode
 *
 * writes branch numbers, flow directions and different results to .dxf file
 *
 * \param PathFile project name
 * \param mode determines what should be printed as .dxf
 * \return int error code
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/

int Print2dxf(ShowMode mode) {
	/* Local variables */
	int error = 0;
	string writefile;
	ofstream outData;

	switch (mode) {
	case ShowMode::TubesPoints:
		writefile = PathFile + "TubesPoints.dxf";
		break;
	case ShowMode::BranchesNodes:
		writefile = PathFile + "Branches.dxf";
		break;
	case ShowMode::Flow:
		writefile = PathFile + "Flow" + to_string(Base.iterg) + ".dxf";
		break;
	case ShowMode::Arrows:
		writefile = PathFile + "step" + to_string(Base.iterg) + ".dxf";
		break;
	case ShowMode::SafetyFactor:
		writefile = PathFile + "_result_safety.dxf";
		break;
	case  ShowMode::SteamQuality:
		writefile = PathFile + "_result_X.dxf";
		break;
	case ShowMode::VoidFraction:
		writefile = PathFile + "_result_void.dxf";
		break;
	case ShowMode::MassVel:
		writefile = PathFile + "_result_massvel.dxf";
		break;
	case ShowMode::Velocity:
		writefile = PathFile + "_result_velocity.dxf";
		break;
	case ShowMode::FlowPattern:
		writefile = PathFile + "_result_pattern.dxf";
		break;
	case ShowMode::DPdynPerLength:
		writefile = PathFile + "_result_dpdynamic.dxf";
		break;
	case ShowMode::ResistanceFactor:
		writefile = PathFile + "_result_resistance.dxf";
		break;
	}
	outData.open(writefile.c_str());
	if (!outData.good()) {
		cout << "\n cannot open file " << writefile << endl;
		exit(1);
	}
	else {
		error = DXFWrite(outData, mode);
		outData.close();
	}//file opened correctly
	return error;
}

