/*! \file   SaveResults.cpp */


//#include "stdafx.h"
#include <cstdlib>
#undef MAINFUNCTION
#include "CommonHeader.h"

using namespace std;

int SaveResults() {
	ofstream result, dchange;

	string rname = PathFile + ".res";
	// Set exceptions to be thrown on failure
	result.exceptions(ofstream::failbit | ofstream::badbit);
	/// 1.) opening and writing data to result file (....res)
	try {
		result.open(rname.c_str());
	}
	catch (system_error& e) {
		cerr << rname << " : " << e.code().message() << endl;
		exit(1);
	}
	result << "  iTb PtIn PtOut    g       xIn    xOut    void    beta     UorS HeatFlux" <<
		"    dpIn     dpOut     dpdyn    dpstat   rhoIn  rhoOut" << //rhoMean
		"     w    safety-  flow\n";
	result << "                 [kg/s]     [%]     [%]     [%]                 [kW/m2]" <<
		"      [Pa]      [Pa]      [Pa]      [Pa] [kg/m3] [kg/m3]" //[kg/m3]
		<< "   [m/s]  factor  chart ";
	result.setf(ios::fixed, ios::floatfield);
	/// for each tube following data is saved:\n
	for (auto& iTube : Tubes) {
		//         prot << "iTb " << iTube.Number << " q " << iTube.q << " xOut "<<iTube.xOut<< endl;
		//         if (iTube.xIn < 0.) {
		//            iTube.xIn = 0.;
		//         }
		//         if (iTube.xOut < 0.) {
		//            iTube.xOut = 0.;
		//         }
		iTube.Criteria();
		iTube.Steiner(iTube.pPaOut, iTube.xOut);

		/// tube number (it's the same as in input file),\n
		/// point number of inlet (it's the same as in input file),\n
		/// point number of outlet (it's the same as in input file),\n
		/// flow in kg/s,\n
		/// steam quality at inlet [%],\n
		/// steam quality at outlet [%],\n 
		/// void fraction at outlet [%],\n
		/// angle to preceding tube [deg],\n
		/// U or S arrangement [-],\n
		/// HeatFlux [kW/m2],\n 
		/// pressure difference at inlet (dpIn) [Pa],\n 
		/// pressure difference at outlet (dpOut) [Pa],\n 
		/// dynamic pressure difference incl. acceleration pressure difference (dpdyn) [Pa],\n
		/// static pressure difference (- on downward flow) (dpstat) [Pa],\n 
		/// density at tube inlet (rhoIn) [kg/m3],\n
		/// density at tube outlet (rhoOut) [kg/m3],\n
		/// actual velocity at tube outlet (velOut) [m/s],\n
		/// lowest safety factor against different criteria (SafetyFactor)[-],\n
		/// flow pattern according to Steiner.

		result << "\n" << setw(5) << iTube.Number
			<< setw(5) << iTube.PointIn
			<< setw(5) << iTube.PointOut
			<< setw(8) << setprecision(2) << Branches[iTube.NbBr].g
			<< setw(8) << setprecision(2) << iTube.xIn * 100.
			<< setw(8) << setprecision(2) << iTube.xOut * 100.
			<< setw(8) << setprecision(2) << iTube.VoidFractionOut * 100.
			<< setw(8) << setprecision(2) << iTube.beta;
		if (iTube.UorS == USArrangement::S) {
			result << "       S";
		}
		else if (iTube.UorS == USArrangement::U) {
			result << "       U";
		}
		else {
			result << "       -";
		}
		result << setw(8) << setprecision(2) << iTube.HeatFlux
			<< setw(10) << setprecision(2) << iTube.dpIn
			<< setw(10) << setprecision(2) << iTube.dpOut
			<< setw(10) << setprecision(2) << iTube.dpdyn
			<< setw(10) << setprecision(2) << iTube.dpstat
			<< setw(8) << setprecision(2) << iTube.rhoIn
			<< setw(8) << setprecision(2) << iTube.rhoOut
			//                << setw(8) << setprecision(2) << iTube.rhoMean
			//<< setw(8) << setprecision(2) << iTube.velWater
			<< setw(8) << setprecision(2) << iTube.velOut;
		if (iTube.q > 0.) {
			result << setw(8) << setprecision(2) << iTube.SafetyFactor;
		}
		else {
			result << "       -";
		}
		if (iTube.xOut > 1e-3 && fabs(iTube.Height / iTube.Length) < 0.174) {
			switch (iTube.steiner) {
			case FlowPattern::stratified:
				result << " stratified";
				break;
			case FlowPattern::wavy:
				result << " wavy";
				break;
			case FlowPattern::bubble:
				result << " bubble";
				break;
			case FlowPattern::plugSlug:
				result << " plug or slug";
				break;
			case FlowPattern::mist:
				result << " mist";
				break;
			case FlowPattern::annular:
				result << " annular";
				break; //annular flow
			default:
				result << " -------";
			}
		}
	}
	result.close();

	/// 2.)saving geometry to dxf-file coloring according SafetyFactor\n
	/// 3.)saving geometry to dxf-file coloring according steam quality\n
	/// 4.)saving geometry to dxf-file coloring according void fraction\n
	/// 5.)saving geometry to dxf-file coloring according mass velocity\n
	/// 6.)saving geometry to dxf-file coloring according mass velocity\n
	/// 7.)saving geometry to dxf-file coloring according flow pattern(Steiner)\n
	/// 8.)saving geometry to dxf-file coloring according dpdyn / length(it helps to identify bottle necks)\n
	/// 9.)saving geometry to dxf-file coloring according dpdyn / length / mass velocity^2(it helps to identify bottle necks)\n
	/// in dxf-files the flow direction is indicated by arrows in the centers of the tubes

	Print2dxf(ShowMode::Flow);              // showing tubes, coloring according Flow 
	Print2dxf(ShowMode::SafetyFactor);     // showing tubes, coloring according SafetyFactor
	Print2dxf(ShowMode::SteamQuality);     // showing tubes, coloring according steam quality
	Print2dxf(ShowMode::VoidFraction);     // showing tubes, coloring according void fraction
	Print2dxf(ShowMode::MassVel);          // showing tubes, coloring according mass velocity
	Print2dxf(ShowMode::Velocity);         // showing tubes, coloring according velocity
	Print2dxf(ShowMode::FlowPattern);      // showing tubes, coloring according flow pattern(Steiner)
	Print2dxf(ShowMode::DPdynPerLength);       // showing tubes, coloring according dpdyn / length(to find bottle neck)
	Print2dxf(ShowMode::ResistanceFactor); // showing tubes, coloring according dpdyn / length / mass velocity^2(to find bottle neck)

	return 0;
}
