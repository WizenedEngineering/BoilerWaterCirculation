//#include "stdafx.h"

#include "CommonHeader.h"

/*
*  calculates the pressure difference in a tube
*
*  ResFactIn resistance factor at inlet [-]
*  ResFactOut resistance factor at outlet [-]
*  g Flow in all parallel tubes [kg/s]
* return int error code
*/
int _tube::dpTube(double ResFactIn, double ResFactOut, double g) {
	int error = 0;
	double rhoInSect = 0.; // density at inlet [kg/m3]
	double rhoOutSect = 0.; // density at outlet [kg/m3]
	double rhoMeanSect = 0.; // mean density [kg/m3]
	//   double pPaInSect; // absolute pressure at inlet [Pa]
	double FrictCoeffAdd = 0.; //friction coefficient [-]
	// if single place: friction coefficient of this particular place
	// if tube: additional friction coefficient
	double enthInSect = 0.; // inlet enthalpy [kJ/kg]
	double enthOutSect = 0.; // outlet enthalpy [kJ/kg]
	double LengthOrifice = 16.e-3; // length (thickness) of orifice [m]
	double tSatOutSect = 0.; //saturation temperature at outlet [K]
	double volWOutSect = 0.; //spec. volume of saturated water at outlet  [m3/kg]
	double volSOutSect = 0.; //spec. volume of saturated steam at outlet [m3/kg]
	double SurfTensOutSect = 0.; //surface tension at outlet [N/m]
	double dynVisSOutSect = 0.; // dynamic viscosity of saturated steam at outlet [Pa s]
	double dynVisWOutSect = 0.; //dynamic viscosity of saturated water at outlet [Pa s]
	double dynVisOutSect = 0.; //dynamic viscosity at outlet [Pa s], single phase
	double enthWSatOutSect = 0.; //spec. enthalpy of saturated water at outlet [kJ/kg]
	double enthSSatOutSect = 0.; //spec. enthalpy of saturated steam at outlet [kJ/kg]
	double pPaOutSect = 0.; // absolute pressure at outlet [Pa]
	double dpDynSect = 0.; // dynamic pressure drop [Pa]
	double dpStatSect = 0.; // static pressure difference [Pa]
	double xInSect = 0.; //steam quality at inlet
	double xOutSect = 0.; //steam quality at outlet
	double VoidInSect = 0.; // void fraction at inlet
	double VoidOutSect = 0.; // void fraction at outlet 
	double VelSect = 0.; // mean velocity [m/s]
	bool PhaseChange = false;

	Flow = g / NoParallel;
	MassVel = Flow / area;
	//splitting tube in shorter sections to fulfill PMean = (PIn+POut)/2 and
	dpdyn = 0.;
	dpstat = 0.;
	dpIn = 0.;
	dpOut = 0.;
	// calculation of inlet
	enthInSect = EnthIn;
	pPaOutSect = pPaIn;
	if (Base.showDPTube) {
		prot << "\n tube number " << Number << " Flow " << Flow << " massvel " << MassVel;
		prot << "\nInlet : PPaIn " << pPaOutSect << endl;
	}
	enthOutSect = EnthIn;
	rhoIn = 0.;
	rhoOutSect = 0.;
	ksiIn = 0.;
	ksiTube = 0.;
	ksiOut = 0.;

	/**
	 * The inlet, outlet, inlet/outlet orifice and bend pressure losses are assumed at a single place.\n
	 * This places have "zero" length and height, therefore dynamic pressure difference only.
	 *
	 *  1) first handling of inlet resistance factor, if it is not zero
	 */
	if (fabs(ResFactIn) > 1e-6) {
		if (Base.showDPTube) {
			prot << "\n ResFactIn" << ResFactIn;
		}
		/**
		 * the inlet resistance factor adds to ksiIn
		 */
		ksiIn = ResFactIn;
		xOutSect = 0.;
		error = dPSection(TbRegion::SinglePlace, // single place (inlet, outlet,orifice, valve) or tube = 0
			ResFactIn, //friction coefficient [-]
			// if single place: friction coefficient of this particular place
			// if tube: additional friction coefficient
			enthInSect, // inlet enthalpy [kJ/kg]
			enthOutSect, // outlet enthalpy [kJ/kg]
			0., // length (thickness) of orifice [m]
			PhaseChange,
			tSatOutSect, //saturation temperature at outlet [K]
			volWOutSect, //spec. volume of saturated water at outlet  [m3/kg]
			volSOutSect, //spec. volume of saturated steam at outlet [m3/kg]
			SurfTensOutSect, //surface tension at outlet [N/m]
			dynVisSOutSect, // dynamic viscosity of saturated steam at outlet [Pa s]
			dynVisWOutSect, //dynamic viscosity of saturated water at outlet [Pa s]
			dynVisOutSect, //dynamic viscosity at outlet [Pa s], single phase
			enthWSatOutSect, //spec. enthalpy of saturated water at outlet [kJ/kg]
			enthSSatOutSect, //spec. enthalpy of saturated steam at outlet [kJ/kg]
			pPaOutSect, // absolute pressure at outlet [Pa]
			dpDynSect, // dynamic pressure drop [Pa]
			dpStatSect, // static pressure difference [Pa]
			xInSect, //steam quality at inlet
			xOutSect, //steam quality at outlet
			rhoInSect, // density at inlet [kg/m3]
			rhoOutSect, // density at outlet [kg/m3]
			// if no previous section, it has to be 0.
			VoidInSect, // void fraction at inlet
			VoidOutSect, // void fraction at outlet 
			rhoMeanSect, // mean density [kg/m3]
			VelSect); // mean velocity [m/s]
		if (Base.showDPTube) {
			prot << "\nksiIn dpDynSect " << dpDynSect << " dpStatSect " << dpStatSect;
			prot << "\nPPaOutSect " << pPaOutSect << " voidIn " << VoidInSect;
		}
		dpIn = dpDynSect;
		rhoIn = rhoInSect;
		xIn = xInSect;
		VoidFractionIn = VoidInSect;
	}
	/**
	 * 2) if there is an orifice at inlet adding the pressure loss
	 */
	if (DiaOrificeIn > 0.) {
		if (Base.showDPTube) {
			prot << "\n Orifice In";
		}
		error = dPSection(TbRegion::Orifice, // single place (inlet, outlet,orifice , valve) or tube = 0
			0., //friction coefficient [-]
			// if single place: friction coefficient of this particular place
			// if tube: additional friction coefficient
			enthInSect, // inlet enthalpy [kJ/kg]
			enthOutSect, // outlet enthalpy [kJ/kg]
			LengthOrifice, // length (thickness) of orifice [m]
			PhaseChange,
			tSatOutSect, //saturation temperature at outlet [K]
			volWOutSect, //spec. volume of saturated water at outlet  [m3/kg]
			volSOutSect, //spec. volume of saturated steam at outlet [m3/kg]
			SurfTensOutSect, //surface tension at outlet [N/m]
			dynVisSOutSect, // dynamic viscosity of saturated steam at outlet [Pa s]
			dynVisWOutSect, //dynamic viscosity of saturated water at outlet [Pa s]
			dynVisOutSect, //dynamic viscosity at outlet [Pa s], single phase
			enthWSatOutSect, //spec. enthalpy of saturated water at outlet [kJ/kg]
			enthSSatOutSect, //spec. enthalpy of saturated steam at outlet [kJ/kg]
			pPaOutSect, // absolute pressure at outlet [Pa]
			dpDynSect, // dynamic pressure drop [Pa]
			dpStatSect, // static pressure difference [Pa]
			xInSect, //steam quality at inlet
			xOutSect, //steam quality at outlet
			rhoInSect, // density at inlet [kg/m3]
			rhoOutSect, // density at outlet [kg/m3]
			// if no previous section, it has to be 0.
			VoidInSect, // void fraction at inlet
			VoidOutSect, // void fraction at outlet 
			rhoMeanSect, // mean density [kg/m3]
			VelSect); // mean velocity [m/s]
		if (Base.showDPTube) {
			prot << "\norificeIn: dpDynSect " << dpDynSect << " dpStatSect " << dpStatSect;
			prot << "\nPPaOutSect " << pPaOutSect << " voidIn " << VoidInSect;
		}
		dpdyn += dpDynSect;
		dpstat += dpStatSect;
		if (rhoIn < 1.) {
			rhoIn = rhoInSect;
			xIn = xInSect;
			VoidFractionIn = VoidInSect;
		}
		/**
		 * the inlet orifice resistance adds to ksiIn
		 */
		ksiIn += dpDynSect * 2. * rhoMeanSect / MassVel / MassVel;
	}
	/**
	 * 3) if there is an angle between this tube and the preceding one or if the tube itself is a bend/elbow,\n
	 *  the pressure difference of the bend is taken into account\n
	 *  the length of the bend is handled in the tube section
	 * 
	 */

	if (beta > 1.) {
		if (Base.showDPTube) {
			prot << "\n bend";
		}
		error = dPSection(TbRegion::Bend, // single place (inlet, outlet,orifice , valve) or tube = 0
			0., //friction coefficient [-]
			// if single place: friction coefficient of this particular place
			// if tube: additional friction coefficient
			enthInSect, // inlet enthalpy [kJ/kg]
			enthOutSect, // outlet enthalpy [kJ/kg]
			0, // length (thickness) of orifice [m]
			PhaseChange,
			tSatOutSect, //saturation temperature at outlet [K]
			volWOutSect, //spec. volume of saturated water at outlet  [m3/kg]
			volSOutSect, //spec. volume of saturated steam at outlet [m3/kg]
			SurfTensOutSect, //surface tension at outlet [N/m]
			dynVisSOutSect, // dynamic viscosity of saturated steam at outlet [Pa s]
			dynVisWOutSect, //dynamic viscosity of saturated water at outlet [Pa s]
			dynVisOutSect, //dynamic viscosity at outlet [Pa s], single phase
			enthWSatOutSect, //spec. enthalpy of saturated water at outlet [kJ/kg]
			enthSSatOutSect, //spec. enthalpy of saturated steam at outlet [kJ/kg]
			pPaOutSect, // absolute pressure at outlet [Pa]
			dpDynSect, // dynamic pressure drop [Pa]
			dpStatSect, // static pressure difference [Pa]
			xInSect, //steam quality at inlet
			xOutSect, //steam quality at outlet
			rhoInSect, // density at inlet [kg/m3]
			rhoOutSect, // density at outlet [kg/m3]
			// if no previous section, it has to be 0.
			VoidInSect, // void fraction at inlet
			VoidOutSect, // void fraction at outlet 
			rhoMeanSect, // mean density [kg/m3]
			VelSect); // mean velocity [m/s]
		if (Base.showDPTube) {
			prot << "\nbend: dpDynSect " << dpDynSect << " dpStatSect " << dpStatSect;
			prot << "\nPPaOutSect " << pPaOutSect << " voidIn " << VoidInSect;
		}
		dpdyn += dpDynSect;
		dpstat += dpStatSect;
		if (rhoIn < 1.) {
			rhoIn = rhoInSect;
			xIn = xInSect;
			VoidFractionIn = VoidInSect;
		}
		/**
		* however calculated separately the bend resistance factor adds to ksiTube
		*/
		ksiTube += dpDynSect * 2. * rhoMeanSect / MassVel / MassVel;
	}
	/**
	* 4) loop through sections and calculate pressure difference of tube without inlet and outlet
	*/
	FrictCoeffAdd = ksiAdd / NoSections;
	for (int i = 1; i <= NoSections; ++i) {
		enthInSect = enthOutSect;
		enthOutSect += HeatSection / Flow;
		if (Base.showDPTube) {
			prot << "\n Tube Section " << i << " NoSections " << NoSections;
		}
		//cout<< " Tube Section " << i<<endl;
		error = dPSection(TbRegion::Tube, // single place (inlet, outlet,orifice, valve) or tube = 0
			FrictCoeffAdd, //friction coefficient [-]
			// if single place: friction coefficient of this particular place
			// if tube: additional friction coefficient
			enthInSect, // inlet enthalpy [kJ/kg]
			enthOutSect, // outlet enthalpy [kJ/kg]
			0., // length (thickness) of orifice [m]
			PhaseChange,
			tSatOutSect, //saturation temperature at outlet [K]
			volWOutSect, //spec. volume of saturated water at outlet  [m3/kg]
			volSOutSect, //spec. volume of saturated steam at outlet [m3/kg]
			SurfTensOutSect, //surface tension at outlet [N/m]
			dynVisSOutSect, // dynamic viscosity of saturated steam at outlet [Pa s]
			dynVisWOutSect, //dynamic viscosity of saturated water at outlet [Pa s]
			dynVisOutSect, //dynamic viscosity at outlet [Pa s], single phase
			enthWSatOutSect, //spec. enthalpy of saturated water at outlet [kJ/kg]
			enthSSatOutSect, //spec. enthalpy of saturated steam at outlet [kJ/kg]
			pPaOutSect, // absolute pressure at outlet [Pa]
			dpDynSect, // dynamic pressure drop [Pa]
			dpStatSect, // static pressure difference [Pa]
			xInSect, //steam quality at inlet
			xOutSect, //steam quality at outlet
			rhoInSect, // density at inlet [kg/m3]
			rhoOutSect, // density at outlet [kg/m3]
			// if no previous section, it has to be 0.
			VoidInSect, // void fraction at inlet
			VoidOutSect, // void fraction at outlet 
			rhoMeanSect, // mean density [kg/m3]
			VelSect); // mean velocity [m/s]
		if (Base.showDPTube) {
			prot << "\n section: dpDynSect " << dpDynSect << " dpStatSect " << dpStatSect;
			prot << "\n rhooutsect" << rhoOutSect << " vel " << MassVel / rhoOutSect;
			prot << "\n PPaOutSect " << pPaOutSect << " voidOut " << VoidOutSect;
		}
		dpdyn += dpDynSect;
		dpstat += dpStatSect;
		if (i == 1 && rhoIn < 1.) {
			rhoIn = rhoInSect;
			xIn = xInSect;
			VoidFractionIn = VoidInSect;
		}
		/**
		 * ksiTube is added up over all sections
		 */
		ksiTube += dpDynSect * 2. * rhoMeanSect / MassVel / MassVel;
		if (Base.showDPTube) {
			prot << "\n tube: dpDyn " << dpdyn << " dpStat " << dpstat;
		}

	}
	/**
	 * 5) if there is an orifice at outlet adding the pressure loss
	 */
	if (DiaOrificeOut > 0.) {
		if (Base.showDPTube) {
			prot << "\n Orifice out";
		}
		// cout<< " Orifice out"<<endl;
		error = dPSection(TbRegion::Orifice, // single place (inlet, outlet,orifice, valve) or tube = 0
			0., //friction coefficient [-]
			// if single place: friction coefficient of this particular place
			// if tube: additional friction coefficient
			enthInSect, // inlet enthalpy [kJ/kg]
			enthOutSect, // outlet enthalpy [kJ/kg]
			-0.016, // length (thickness) of orifice [m]
			PhaseChange,
			tSatOutSect, //saturation temperature at outlet [K]
			volWOutSect, //spec. volume of saturated water at outlet  [m3/kg]
			volSOutSect, //spec. volume of saturated steam at outlet [m3/kg]
			SurfTensOutSect, //surface tension at outlet [N/m]
			dynVisSOutSect, // dynamic viscosity of saturated steam at outlet [Pa s]
			dynVisWOutSect, //dynamic viscosity of saturated water at outlet [Pa s]
			dynVisOutSect, //dynamic viscosity at outlet [Pa s], single phase
			enthWSatOutSect, //spec. enthalpy of saturated water at outlet [kJ/kg]
			enthSSatOutSect, //spec. enthalpy of saturated steam at outlet [kJ/kg]
			pPaOutSect, // absolute pressure at outlet [Pa]
			dpDynSect, // dynamic pressure drop [Pa]
			dpStatSect, // static pressure difference [Pa]
			xInSect, //steam quality at inlet
			xOutSect, //steam quality at outlet
			rhoInSect, // density at inlet [kg/m3]
			rhoOutSect, // density at outlet [kg/m3]
			// if no previous section, it has to be 0.
			VoidInSect, // void fraction at inlet
			VoidOutSect, // void fraction at outlet 
			rhoMeanSect, // mean density [kg/m3]
			VelSect); // mean velocity [m/s]
		if (Base.showDPTube) {
			prot << "\n orificeout dpDynSect " << dpDynSect << " dpStatSect " << dpStatSect;
			prot << "\n PPaOutSect " << pPaOutSect << " voidOut " << VoidOutSect;
		}
		dpdyn += dpDynSect;
		dpstat += dpStatSect;
		/**
		* the outlet orifice resistance adds to ksiOut
		*/
		ksiOut += dpDynSect * 2. * rhoMeanSect / MassVel / MassVel;
	}
	/**
	* 6) finally handling of outlet resistance factor .
	*/
	if (fabs(ResFactOut) > 1e-6) {
		if (Base.showDPTube) {
			prot << "\n ResFactOut" << ResFactOut;
		}
		/**
		* the outlet resistance factor adds to ksiOut
		*/
		ksiOut += ResFactOut;
		// cout<< " ResFactOut"<<endl;
		error = dPSection(TbRegion::SinglePlace, // single place (inlet, outlet,orifice, valve) or tube = 0
			ResFactOut, //friction coefficient [-]
			// if single place: friction coefficient of this particular place
			// if tube: additional friction coefficient
			enthInSect, // inlet enthalpy [kJ/kg]
			enthOutSect, // outlet enthalpy [kJ/kg]
			0, // length (thickness) of orifice [m]
			PhaseChange,
			tSatOutSect, //saturation temperature at outlet [degC]
			volWOutSect, //spec. volume of saturated water at outlet  [m3/kg]
			volSOutSect, //spec. volume of saturated steam at outlet [m3/kg]
			SurfTensOutSect, //surface tension at outlet [N/m]
			dynVisSOutSect, // dynamic viscosity of saturated steam at outlet [Pa s]
			dynVisWOutSect, //dynamic viscosity of saturated water at outlet [Pa s]
			dynVisOutSect, //dynamic viscosity at outlet [Pa s], single phase
			enthWSatOutSect, //spec. enthalpy of saturated water at outlet [kJ/kg]
			enthSSatOutSect, //spec. enthalpy of saturated steam at outlet [kJ/kg]
			pPaOutSect, // absolute pressure at outlet [Pa]
			dpDynSect, // dynamic pressure drop [Pa]
			dpStatSect, // static pressure difference [Pa]
			xInSect, //steam quality at inlet
			xOutSect, //steam quality at outlet
			rhoInSect, // density at inlet [kg/m3]
			rhoOutSect, // density at outlet [kg/m3]
			// if no previous section, it has to be 0.
			VoidInSect, // void fraction at inlet
			VoidOutSect, // void fraction at outlet 
			rhoMeanSect, // mean density [kg/m3]
			VelSect); // mean velocity [m/s]
		if (Base.showDPTube) {
			prot << "\nKsiout: dpDynSect " << dpDynSect << " dpStatSect " << dpStatSect;
			prot << "\n PPaOutSect " << pPaOutSect << " voidOut " << VoidOutSect;
		}
		dpOut = dpDynSect;
	}
	velIn = MassVel / rhoIn;
	//   prot << "velIn " << velIn << endl;
	rhoOut = rhoOutSect;
	EnthOut = enthOutSect;
	pPaOut = pPaOutSect;
	velOut = MassVel / rhoOut;
	xOut = xOutSect;
	VoidFractionOut = VoidOutSect;
	if (xOut <= 1e-6) {
		velWater = velOut;
	}
	else {
		velWater = MassVel * volWOutSect;
	}
	if (fabs(Height) < 1e-3) {
		rhoMean = (rhoIn + rhoOut) / 2.;
	}
	else {
		rhoMean = dpstat / 9.80665 / Height;
	}

	return error;
}

