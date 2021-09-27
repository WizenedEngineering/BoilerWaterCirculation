//#include "stdafx.h"
#include "CommonHeader.h"

int _tube::dPSection(// returns error code: 0 all ok, 1 phase change, 2 out of turbulent range
		  // no phase change in tube or single place
	TbRegion region, // single place (inlet, outlet, orifice=1, bend, = valve) or tube = 0
	double FrictCoeff, //friction coefficient [-]
	// if single place: friction coeff of this particular place
	// if tube: additional friction coeff
	double enthInSect, // inlet enthalpy [kJ/kg]
	double enthOutSect, // outlet enthalpy [kJ/kg]
	double LengthOrifice, // length (thickness) of orifice [m]
	bool& PhaseChange,
	double& tSatOutSect, //saturation temperature at outlet [K]
	double& volWOutSect, //spec. volume of saturated water at outlet  [m3/kg]
	double& volSOutSect, //spec. volume of saturated steam at outlet [m3/kg]
	double& SurfTensOutSect, //surface tension at outlet [N/m]
	double& dynVisSOutSect, // dynamic viscosity of saturated steam at outlet [Pa s]
	double& dynVisWOutSect, //dynamic viscosity of saturated water at outlet [Pa s]
	double& dynVisOutSect, //dynamic viscosity at outlet [Pa s], single phase
	double& enthWSatOutSect, //spec. enthalpy of saturated water at outlet [kJ/kg]
	double& enthSSatOutSect, //spec. enthalpy of saturated steam at outlet [kJ/kg]
	double& pPaOutSect, // absolute pressure at outlet [Pa]
	double& dpDynSect, // dynamic pressure drop [Pa]
	double& dpStatSect, // static pressure difference [Pa]
	double& xInSect, //steam quality at inlet
	double& xOutSect, //steam quality at outlet
	double& rhoInSect, // density at inlet [kg/m3]
	double& rhoOutSect, // density at outlet [kg/m3]
	// if no previous section, it has to be 0.
	double& VoidInSect, // void fraction at inlet
	double& VoidOutSect, // void fraction at outlet 
	double& rhoMeanSect, // mean velocity [m/s]
	double& VelSect) // mean velocity [m/s]
{
	int error = 0;
	double tSatIn, pMPaIn;
	double enthWSatInSect, //enthalpy of saturated water at inlet (from previous section) [kJ/kg]
		enthSSatInSect; //enthalpy of saturated steam at inlet (from previous section) [kJ/kg]
	double pPaInSect; // absolute pressure at inlet [Pa]
/**
 * Most of the function parameters are properties at outlet of tube section.
 *
 * to avoid double calculation of the values at inlet the outlet values of preceding section are used.
 *
 * only for the first section in a tube the properties have be calculated indicated by rhoOutSect = 0
 */

 /**
  * if a phase change (f.i from single phase water to 2 phase) occurred in the previous section the inlet conditions have to be calculated anew.\n
  */
	if (PhaseChange) {
		rhoOutSect = 0.;
		PhaseChange = false;
	}
	pPaInSect = pPaOutSect;
	if (rhoOutSect < 1. ||
		enthWSatOutSect < 1. ||
		enthSSatOutSect < 1.) {
		pMPaIn = pPaInSect * 1e-6;
		tSatIn = H2O::satTemp(pMPaIn);
		enthWSatInSect = H2O::enth(tSatIn, pMPaIn, WATER);
		enthSSatInSect = H2O::enth(tSatIn, pMPaIn, STEAM);
	}
	else {
		enthWSatInSect = enthWSatOutSect;
		enthSSatInSect = enthSSatOutSect;
	}
	//  	prot << "\n enthInSect " << enthInSect << " enthWSatInSect " << enthWSatInSect;
	if (enthInSect <= enthWSatInSect + 1e-3) {
		// assuming short tube sections -> the section with phase change will
		// be calculated as single phase at inlet (water ->boiling)
		//		prot << "single phase water" << endl;
		xInSect = 0.;
		xOutSect = 0.;
		VoidInSect = 0.;
		VoidOutSect = 0.;
		error = dPSinglePhase(// returns error code: 0 all ok,
				  // no phase change in tube or single place
			region, // single place (inlet, outlet, valve :1, orifice, valve) or tube = 0
			FrictCoeff, //friction coefficient [-]
			// if single place: friction coeff of this particular place
			// if tube: additional friction coeff like deflection
			enthInSect, // inlet enthalpy [kJ/kg]
			enthOutSect, // outlet enthalpy [kJ/kg]
			LengthOrifice, // Length of orifice [m]
			WATER, // index for water or steam
			PhaseChange, // if outlet above saturation = true
			dynVisOutSect, // dyn. viscosity at inlet [Pa s]
			pPaOutSect, // absolute pressure at outlet [Pa]
			dpDynSect, // dynamic pressure drop [Pa]
			dpStatSect, // static pressure difference [Pa]
			xOutSect, //steam quality at outlet
			rhoInSect, // density at inlet [kg/m3]
			rhoOutSect, // density at outlet [kg/m3]
			rhoMeanSect, // mean density [kg/m3]
			VelSect); // mean velocity [m/s]
		if (Base.showDPTubeDetail) {
			prot << "\nsingle phase water" << endl;
		}
	}
	else if (enthInSect >= enthSSatInSect - 1e-3) { // single phase steam
		xInSect = 1.;
		xOutSect = 1.;
		VoidInSect = 1.;
		VoidOutSect = 1.;
		error = dPSinglePhase(// returns error code: 0 all ok,
				  // no phase change in tube or single place
			region, // single place (inlet, outlet, valve :1, orifice, valve) or tube = 0
			FrictCoeff, //friction coefficient [-]
			// if single place: friction coeff of this particular place
			// if tube: additional friction coeff like deflection
			enthInSect, // inlet enthalpy [kJ/kg]
			enthOutSect, // outlet enthalpy [kJ/kg]
			LengthOrifice, // Length of orifice [m]
			STEAM, // index for water or steam
			PhaseChange, // if outlet below saturation = true
			dynVisOutSect, // dyn. viscosity at inlet [Pa s]
			pPaOutSect, // absolute pressure at outlet [Pa]
			dpDynSect, // dynamic pressure drop [Pa]
			dpStatSect, // static pressure difference [Pa]
			xOutSect, //steam quality at outlet
			rhoInSect, // density at inlet [kg/m3]
			rhoOutSect, // density at outlet [kg/m3]
			rhoMeanSect, // mean density [kg/m3]
			VelSect); // mean velocity [m/s]
		if (Base.showDPTubeDetail) {
			prot << "\nsingle phase steam" << endl;
		}
	}
	else { //two-phase flow
		error = dPTwoPhase(region, FrictCoeff,
			enthInSect, enthOutSect, LengthOrifice, PhaseChange,
			tSatOutSect, volWOutSect, volSOutSect, SurfTensOutSect, dynVisSOutSect,
			dynVisWOutSect, enthWSatOutSect, enthSSatOutSect, pPaOutSect, dpDynSect, dpStatSect,
			xInSect, xOutSect, rhoInSect, rhoOutSect, VoidInSect, VoidOutSect, rhoMeanSect, VelSect);
	} // end two-phase
	return error;
}
