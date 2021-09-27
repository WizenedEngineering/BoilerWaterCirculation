//#include "stdafx.h"
#include "CommonHeader.h"


// It is assumed that the tube sections are small enough that the properties
// are mean value between inlet and outlet

int _tube::dPSinglePhase(// returns error code: 0 all ok,
		  // no phase change in tube or single place
	TbRegion region, // single place (inlet, outlet, valve :1, orifice, valve) or tube = 0
	double FrictCoeff, //friction coefficient [-]
	// if single place: friction coefficient of this particular place
	// if tube: additional friction coefficient like deflection
	double enthInSect, // inlet enthalpy [kJ/kg]
	double enthOutSect, // outlet enthalpy [kJ/kg]
	double LengthOrifice, // Length of orifice [m]
	int index, // index for water or steam
	bool& PhaseChange, // =0, if outlet above saturation = 1
	double& dynViscOutSect, // dyn. viscosity at inlet [Pa s]
	double& pPaOutSect, // absolute pressure at outlet [Pa]
	double& dpDynSect, // dynamic pressure drop [Pa]
	double& dpStatSect, // static pressure difference [Pa]
	double& xOutSect, //steam quality at outlet
	double& rhoInSect, // density at inlet [kg/m3]
	double& rhoOutSect, // density at outlet [kg/m3]
	// if no previous section, it has to be 0.
	double& rhoMeanSect, // mean density [kg/m3]
	double& VelSect) // mean velocity [m/s]
{
	double pPaOutPrev,
		volMean,
		viscMean,
		Reynolds,
		tSatOut,
		volIn;
	double friction, tempIn, tempOut, volOut, pMPaIn, pMPaOut;
	double dynViscInSect, // dyn. viscosity at inlet [Pa s]
		pPaInSect; // absolute pressure at inlet [Pa]
	double enthWSatOut,
		enthSSatOut;
	int i = 0;
	//   cout<<"\n single phase 181 pPaOut "<<pPaOutSect<<" rhoOut "<<rhoOutSect<<" dynVisout "<<dynViscOutSect<<" region "<<region<<" frictcoeff "<<FrictCoeff;
/**
 * Most of the function parameters are properties at outlet of tube section.\n
 * to avoid double calculation of the values at inlet the outlet values of preceding section are used.\n
 * only for the first section in a tube the properties have be calculated indicated by rhoOutSect = 0
 */
	pPaInSect = pPaOutSect;
	if (rhoOutSect < 1. ||
		dynViscOutSect > 1. ||
		dynViscOutSect < 0.) {
		pMPaIn = pPaInSect * 1e-6;
		tempIn = H2O::temp(enthInSect, pMPaIn);
		volIn = H2O::specVol(tempIn, pMPaIn, index);
		rhoInSect = 1. / volIn;
		dynViscInSect = H2O::dynVisc(tempIn, volIn);
		dynViscOutSect = dynViscInSect;
	}
	else {
		rhoInSect = rhoOutSect;
		volIn = 1. / rhoOutSect;
		dynViscInSect = dynViscOutSect;
	}
	dpStatSect = 0.;
	/*!
	* outlet conditions depend on outlet pressure
	*
	* outlet pressure = inlet pressure - dpStat - dpDyn
	*
	* dpStat as well as dpDyn depend on the inlet and outlet conditions\n
	* for simplicity the conditions at inlet and outlet are averaged and used for calculation of the pressure differences
	*
	* Iteration for outlet pressure (and conditions, consequently) starting with outletpressure = inletpressure
	*/
	volOut = volIn;
	PhaseChange = false;
	do {
		pPaOutPrev = pPaOutSect;
		// mean properties
		volMean = (volIn + volOut) / 2.;
		viscMean = (dynViscInSect + dynViscOutSect) / 2.;
		friction = FrictCoeff;
		if (region == TbRegion::Tube) { //tube
			Reynolds = Dia * MassVel / viscMean;
			friction += FrictFact(Reynolds, Base.Rough / Dia) * LengthSection / Dia;
			dpStatSect = 9.80665 * HeightSection / volMean;
		}
		else if (region == TbRegion::Orifice) { // orifice
			friction = ksiOrifice(LengthOrifice, viscMean);
		}
		else if (region == TbRegion::Bend) { //bend
			Reynolds = Dia * MassVel / viscMean;
			friction = ksiBend(Reynolds, beta);
		}
		/**
		 * Next step is calculation of dynamic pressure difference
		 *
		 * As the calculation methods are different the calculation is switched between the tube region like plain tube, orifice, bend.\n
		 * a single place is handled by the additional friction coefficient
		  */
		dpDynSect = ((volOut - volIn) + friction * volMean / 2.) * MassVel * MassVel;
		/**
		 * dpDynSect includes also acceleration pressure drop\n
		 * additional FrictionCoefficient is used with mean properties
		 */
		pPaOutSect = pPaInSect - (dpDynSect + dpStatSect);
		//      cout<<"\nsingle dpdyn "<<dpDynSect<<" dpstat "<<dpStatSect<<" pPaOut "<<pPaOutSect;
		pMPaOut = pPaOutSect * 1e-6;
		/// if outlet pressure is less than 10% of drum pressure or atmospheric pressure, there is something wrong-> calculation will be terminated  
		if (pMPaOut < fmax(0.1, 0.1 * Drum.pMPa)) {
			prot << "\n error: section exit pressure too low";
			cout << "\n error: section exit pressure too low";
			prot << "\n iteration step " << Base.iterg << " Branch " << NbBr <<
				" tube " << Number << " Flow " << Flow << "Mass Velocity " << MassVel;
			switch (region) {
			case TbRegion::Tube: //tube region
				prot << "\n tube region";
				break;
			case TbRegion::Orifice: //orifice region
				prot << "\n orifice region";
				break;
			case TbRegion::Bend: //tube bend region
				prot << "\n tube bend region";
			}
			prot << "\n friction = " << friction << " frictioncoeff " << FrictCoeff;
			prot << "\n volIn = " << volIn << " volOut = " << volOut << " volMean = " << volMean;
			prot << "\n xOutSect = " << xOutSect;

			prot << "\n pPaInSection = " << pPaInSect << " pPaOutSect = " << pPaOutSect;
			prot << "\n dpDynSection = " << dpDynSect << " dpStatSection = " << dpStatSect;
			exit(-1);
		}
		tempOut = H2O::temp(enthOutSect, pMPaOut);
		tSatOut = H2O::satTemp(pMPaOut);
		enthWSatOut = H2O::enth(tSatOut, pMPaOut, WATER);
		enthSSatOut = H2O::enth(tSatOut, pMPaOut, STEAM);
		/// checking for a phase change
		if ((index == WATER && enthOutSect > enthWSatOut) ||
			(index == STEAM && enthOutSect < enthSSatOut)) {
			// phase change
			tempOut = tSatOut;
			PhaseChange = true;
			xOutSect = (enthOutSect - enthWSatOut) / (enthSSatOut - enthWSatOut);
		}
		volOut = H2O::specVol(tempOut, pMPaOut, index);
		dynViscOutSect = H2O::dynVisc(tempOut, volOut);
		/// if after 50 iteration steps no solution is found, there is something wrong and the program terminates
		if (++i > 50) {
			prot << "\n tube " << Number << " dP2phase iteration not converging ";
			prot << "\n flow iteration step " << Base.iterg << "\n Branch " << NbBr << " Flow " << Flow;
			cout << "\n tube " << Number << " dPTwoPhase iteration not converging ";
			exit(-1);
		}
	} while (fabs(pPaOutPrev / pPaInSect - pPaOutSect / pPaInSect) > 1.e-6);
	rhoOutSect = 1. / volOut;
	rhoMeanSect = 1. / volMean;
	VelSect = MassVel * volMean;

	return 0;
}
