//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

// It is assumed that the tube sections are small enough that the properties
// are mean value between inlet and outlet

int _tube::dPTwoPhase(
	TbRegion region, // single place is default but FrictCoeffAdd must be given, orifice, bend, tube
	//        double pPaInSect, // absolute pressure at inlet [Pa]
	double FrictCoeffAdd, //friction coefficient [-]
	// if single place: friction coefficient of this particular place
	// if tube: additional friction coefficient
	double enthIn, // inlet enthalpy [kJ/kg]
	double enthOut, // outlet enthalpy [kJ/kg]
	double LengthOrifice, // length (thickness) of orifice [m]
	bool& PhaseChange,
	double& tSatOut, //saturation temperature at outlet [degC]
	double& volWSatOut, //spec. volume of saturated water at outlet  [m3/kg]
	double& volSSatOut, //spec. volume of saturated steam at outlet [m3/kg]
	double& SurfTensOut, //surface tension at outlet [N/m]
	double& dynVisSSatOut, // dynamic viscosity of saturated steam at outlet [Pa s]
	double& dynVisWSatOut, //dynamic viscosity of saturated water at outlet [Pa s]
	double& enthWSatOut, //spec. enthalpy of saturated water at outlet [kJ/kg]
	double& enthSSatOut, //spec. enthalpy of saturated steam at outlet [kJ/kg]
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
	double& rhoMeanSect, // mean density [kg/m3]
	double& VelSect) // mean velocity [m/s]
{
	/*    Input values */
	/*     All units in SI-units */
	/*     g in kg/s */
	/*     l,h,d in m */
	/*     w,w0 in m/s */
	/*     p, dpdyn, dpstat in Pa */

	/* Local variables */
	double xMean, Reynolds = 0.;
	double pMPaIn, pMPaOut, pPrev, zeta, rhoW, rhoS, B;
	double dynVisIn, dynVisOut;
	double  tSatIn, //saturation temperature at inlet (from previous section) [deg C]
		volSSatIn, //spec. Volume of saturated steam at inlet (from previous section) [m3/kg]
		volWSatIn, //spec. Volume of saturated water at inlet (from previous section) [m3/kg]
		rhoWSatIn, //density  of saturated water at inlet (from previous section) [kg/m3]
		rhoSSatIn, //density  of saturated steam at inlet (from previous section) [kg/m3]
		SurfTensIn, //surface tension at inlet (from previous section) [N/m]
		enthWSatIn, //enthalpy of saturated water at inlet (from previous section) [kJ/kg]
		enthSSatIn, //enthalpy of saturated steam at inlet (from previous section) [kJ/kg]
		dynVisSSatIn, //dynamic viscosity of saturated steam at inlet (from previous section) [Pa s]
		dynVisWSatIn, //dynamic viscosity of saturated water at inlet (from previous section) [Pa s]
		pPaInSect; // absolute pressure at inlet [Pa]
	int i = 0;
	//   prot << " h "<<h<<" l "<<l<<endl;
	//   cout << "\n entIn " << enthIn << " enthOut " << enthOut << "p " << pPaOutSect;
	// conditions at inlet
	//   cout << "\n 63 pPaOut " << pPaOutSect << " rhoOutSect " << rhoOutSect;
/**
 * Most of the function parameters are properties at outlet of tube section.
 *
 * to avoid double calculation of the values at inlet the outlet values of preceding section are used.
 *
 * only for the first section in a tube the properties have be calculated indicated by rhoOutSect = 0
 */
	pPaInSect = pPaOutSect;
	pMPaIn = pPaInSect * 1e-6;
	if (rhoOutSect < 1e-3) { // calculate conditions at inlet
		tSatIn = H2O::satTemp(pMPaIn);
		tSatOut = tSatIn;
		volSSatIn = H2O::specVol(tSatIn, pMPaIn, STEAM);
		volSSatOut = volSSatIn;
		volWSatIn = H2O::specVol(tSatIn, pMPaIn, WATER);
		volWSatOut = volWSatIn;
		rhoWSatIn = 1. / volWSatIn;
		rhoSSatIn = 1. / volSSatIn;
		SurfTensIn = H2O::sigma(tSatIn) / 1000.;
		SurfTensOut = SurfTensIn;
		dynVisSSatIn = H2O::dynVisc(tSatIn, volSSatIn);
		dynVisSSatOut = dynVisSSatIn;
		dynVisWSatIn = H2O::dynVisc(tSatIn, volWSatIn);
		dynVisWSatOut = dynVisWSatIn;
		enthWSatIn = H2O::enth(tSatIn, pMPaIn, WATER);
		enthWSatOut = enthWSatIn;
		enthSSatIn = H2O::enth(tSatIn, pMPaIn, STEAM);
		enthSSatOut = enthSSatIn;
		xInSect = (enthIn - enthWSatIn) / (enthSSatIn - enthWSatIn);
		if (Base.showDPTubeDetail) {
			prot << "\n row 96 xInSect " << xInSect << endl;
		}
		//     cout << "\n110 xin " << xInSect << " pMPa " << pMPaIn << " region " << region;
		if (xInSect <= 0.) {
			xInSect = 0.;
			VoidInSect = 0.;
			rhoInSect = rhoWSatIn;
		}
		else if (xInSect >= 1.) {
			xInSect = 1.;
			VoidInSect = 1.;
			rhoInSect = rhoSSatIn;
		}
		else {
			if (region == TbRegion::Orifice || xInSect < 1e-7) {
				rhoInSect = 1. / (xInSect * volSSatIn + (1. - xInSect) * volWSatIn);
				VoidInSect = xInSect / volWSatIn / (xInSect / volWSatIn + (1. - xInSect) / volSSatIn);
			}
			else {
				if (Base.showDPTubeDetail) {
					prot << "\n row 113 xInSect " << xInSect << endl;
				}
				switch (Base.Method) {
				case 'J':
					rhoInSect = Density_Jirous(xInSect, volWSatIn, volSSatIn, VoidInSect);
					break;
				case 'W':
					rhoInSect = Density_Welzer(xInSect, volWSatIn, volSSatIn, VoidInSect);
					break;
				case 'R':
					rhoInSect = Density_Rouhani(xInSect, rhoWSatIn, rhoSSatIn, SurfTensIn, VoidInSect);
					break;
				case 'C':
					//                rhoOutSect = density_chawla();
					break;
				case 'E':
					rhoInSect = Density_chexal(xInSect, pMPaIn, volWSatIn, volSSatIn,
						dynVisWSatIn, dynVisSSatIn, SurfTensIn, VoidInSect);
					break;
				case 'G':
					rhoOutSect = Density_Woldesemayat(xInSect, rhoWSatIn, rhoSSatIn, SurfTensIn, pMPaIn, VoidInSect);
				}
				if (rhoInSect < rhoSSatIn || rhoInSect > rhoWSatIn) {
					return 1;
				}
			}
		}
	}
	else { //from previous section or tube
		rhoInSect = rhoOutSect;
		xInSect = xOutSect;
		VoidInSect = VoidOutSect;
		tSatIn = tSatOut;
		volSSatIn = volSSatOut;
		volWSatIn = volWSatOut;
		rhoWSatIn = 1. / volWSatIn;
		rhoSSatIn = 1. / volSSatIn;
		SurfTensIn = SurfTensOut;
		dynVisSSatIn = dynVisSSatOut;
		dynVisWSatIn = dynVisWSatOut;
		enthWSatIn = enthWSatOut;
		enthSSatIn = enthSSatOut;
		if (Base.showDPTubeDetail) {
			prot << "\n row 155 xInSect " << xInSect;
		}
	}
	//
	pMPaOut = pMPaIn;
	xOutSect = (enthOut - enthWSatIn) / (enthSSatIn - enthWSatIn); //starting value
	if (Base.showDPTubeDetail) {
		prot << "\nrow 162 xOutSect " << xOutSect << " enthOut " << enthOut << " enthWsatIn " << enthWSatIn << " enthSSatIn " << enthSSatIn;
	}
	dynVisIn = H2O::dynVisc(tSatIn, 1. / rhoInSect); 
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
	do {/*         do */
		pPrev = pPaOutSect;
		//      cout << " \n 172 xout " << xOutSect;
/**
 * In iteration loop first density at outlet is calculated .
  */
		if (xOutSect < 0.) {
			xOutSect = 0.;
			rhoOutSect = 1. / volWSatOut;
			VoidOutSect = 0.;
		}
		else if (xOutSect > 1.) {
			xOutSect = 1.;
			rhoOutSect = 1. / volSSatOut;
			VoidOutSect = 1.;
		}
		else {
			if (region == TbRegion::Orifice || xOutSect < 1e-7) { // homogeneous model
				VoidOutSect = xOutSect / volWSatOut / (xOutSect / volWSatOut + (1. - xOutSect) / volSSatOut);
				rhoOutSect = 1. / (xOutSect * volSSatOut + (1. - xOutSect) * volWSatOut);
			}
			else {
				//            cout << "\n192 xout " << xOutSect << " pMPa " << pMPaOut;
				switch (Base.Method) {
				case 'J':
					rhoOutSect = Density_Jirous(xOutSect, volWSatOut, volSSatOut, VoidOutSect);
					break;
				case 'W':
					rhoOutSect = Density_Welzer(xOutSect, volWSatOut, volSSatOut, VoidOutSect);
					break;
				case 'R':
					rhoOutSect = Density_Rouhani(xOutSect, 1. / volWSatOut, 1. / volSSatOut, SurfTensOut, VoidOutSect);
					break;
				case 'C':
					//               rhoOutSect = density_chawla();
					break;
				case 'E':
					rhoOutSect = Density_chexal(xOutSect, pMPaOut, volWSatOut, volSSatOut,
						dynVisWSatOut, dynVisSSatOut, SurfTensOut, VoidOutSect);
					break;
				case 'G':
					rhoOutSect = Density_Woldesemayat(xOutSect, 1. / volWSatOut, 1. / volSSatOut, SurfTensOut, pMPaOut, VoidOutSect);
				}
			}
		}
		xMean = (xInSect + xOutSect) / 2.;
		dpDynSect = 0.;
		//      cout << "\n region " << region;
/**
 * Next step is calculation of dynamic pressure difference
 *
 * As the calculation methods are different the calculation is switched between the tube region like plain tube, orifice, bend.\n
 * a single place is handled by the additional friction coefficient
  */

		switch (region) {
		case TbRegion::Tube: //tube region
			rhoMeanSect = (rhoInSect + rhoOutSect) / 2.;
			dpStatSect = rhoMeanSect * 9.80665 * HeightSection;
			Reynolds = Dia * MassVel / (dynVisWSatIn + dynVisWSatOut) * 2.;
			zeta = FrictFact(Reynolds, Base.Rough / Dia);
			switch (Base.Method) {
			case 'J':
				dpDynSect = dpdyn_Jirous(zeta, (pMPaIn + pMPaOut) / 2., volSSatIn, volSSatOut, volWSatIn, volWSatOut, xMean);
				break;
			case 'W':
				dpDynSect = dpdyn_Welzer(zeta, (volWSatIn + volWSatOut) / 2., (volSSatIn + volSSatOut) / 2.,
					(enthWSatIn + enthWSatOut) / 2., (enthSSatIn + enthSSatOut) / 2.,
					(dynVisWSatIn + dynVisWSatOut) / 2., xInSect, xOutSect);
				break;
			case 'R':
				dpDynSect = dpdyn_Becker(zeta, xInSect, xOutSect, pPaInSect, pPaOutSect,
					volWSatIn, volWSatOut);
				break;
			case 'C':
				//    dpDynSect =      rhoOutSect = density_chawla();
				break;
			case 'E':
			case 'G':
				dpDynSect = dpdyn_HeatAtlas(pPaInSect, pPaOutSect, volSSatIn,
					volSSatOut, volWSatIn, volWSatOut, dynVisWSatIn,
					dynVisWSatOut, dynVisSSatIn, dynVisSSatOut, xMean);
			}
			if (Base.showDPTubeDetail) {
				prot << "\n frictCoeffAdd " << FrictCoeffAdd << " rhoIn " << rhoInSect << " rhoOut " << rhoOutSect;
				prot << "\n 1/volWin " << 1. / volWSatIn << " 1/volWout " << 1. / volWSatOut << " MassVel " << MassVel;
				prot << "\n head " << (volWSatIn + volWSatOut) / 4. * MassVel * MassVel << endl;
			}
			break;
		case TbRegion::Orifice: //orifice region
			dynVisOut = H2O::dynVisc(tSatOut, 1. / rhoOutSect);
			FrictCoeffAdd = ksiOrifice(LengthOrifice, (dynVisIn + dynVisOut) / 2.);
			if (Base.showDPTubeDetail) {
				prot << "\n orifice dynVisIn " << dynVisIn << " dynVisOut " << dynVisOut;
				prot << "\n friction orifice " << FrictCoeffAdd << endl;
			}
			break;
		case TbRegion::Bend: //tube bend region
			Reynolds = Dia * MassVel / (dynVisWSatIn + dynVisWSatOut) * 2.;
			//           prot << " \n dynVisWSatIn " << dynVisWSatIn << " dynVisWSatOut " << dynVisWSatOut << " Reynolds " << Reynolds;
			xMean = (xInSect + xOutSect) / 2.;
			rhoW = (1. / volWSatIn + 1. / volWSatOut) / 2.;
			rhoS = (1. / volSSatIn + 1. / volSSatOut) / 2.;
			double ksiTubeBend = ksiBend(Reynolds, beta);
			// including 2-phase multiplier VDI Heat Atlas L2.2.3
			if (beta <= 90.) {
				double ksiTubeBend90 = ksiBend(Reynolds, 90.);
				//               prot << "\n ksitube90 " << ksiTubeBend90;
				B = 1. + 2.2 / (ksiTubeBend90 * (2. + RadiusBend / Dia));
			}
			else {
				B = 1. + 2.2 / (ksiTubeBend * (2. + RadiusBend / Dia));
			}
			double phi2 = 1. + (rhoW / rhoS - 1.) * (B * xMean * (1. - xMean) + xMean * xMean);
			//                prot << "\n beta " <<beta <<" radius "<<RadiusBend<<" B " << B << " phi2 " << phi2 << " ksiTubeBend " << ksiTubeBend;
			dpDynSect = ksiTubeBend * phi2 * MassVel * MassVel / 2. / rhoW;
			if (Base.showDPTubeDetail) {
				prot << "\n friction bend " << FrictCoeffAdd;
			}
		}
		//		prot << "\n dp " << dpDynSect << " acceleration " << (1. / rhoOutSect - 1. / rhoInSect)*MassVel*MassVel<<" rhoin "<<rhoInSect<<
		//			" rhoout "<<rhoOutSect <<" dpstatsection " << dpStatSect;
		 /**
		  * dpDynSect includes also acceleration pressure drop\n
		  * additional FrictionCoefficient is distributed to inlet and outlet
		  */
		dpDynSect += ((FrictCoeffAdd / 2. - 2.) / rhoInSect + (FrictCoeffAdd / 2. + 2.) / rhoOutSect)
			/ 2. * MassVel * MassVel;
		/// a new outlet pressure is found and used in next iteration step 
		pPaOutSect = pPaInSect - dpDynSect - dpStatSect;
		pMPaOut = pPaOutSect / 1e6;
		if (pMPaOut < max(0.1, 0.1 * Drum.pMPa)) {
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
			prot << "\n xInSect = " << xInSect << " xOutSect = " << xOutSect;
			prot << "\n pPaInSection = " << pPaInSect << " pPaOutSect = " << pPaOutSect;
			prot << "\n dpDynSection = " << dpDynSect << " dpStatSection = " << dpStatSect;

			exit(-1);
		}
		//      cout << "\n region " << region << " dpdyn " << dpDynSect << " dpstat " << dpStatSect << " frict " << FrictCoeffAdd;
		if (Base.showDPTubeDetail) {
			prot << "\n dpdyn " << dpDynSect << " dpstat " << dpStatSect;
			prot << "\n pPaIn " << pPaInSect << " pPaOutSect " << pPaOutSect;
		}
		tSatOut = H2O::satTemp(pMPaOut);
		volWSatOut = H2O::specVol(tSatOut, pMPaOut, WATER);
		volSSatOut = H2O::specVol(tSatOut, pMPaOut, STEAM);
		SurfTensOut = H2O::sigma(tSatOut) / 1000.;
		dynVisSSatOut = H2O::dynVisc(tSatOut, volSSatOut);
		dynVisWSatOut = H2O::dynVisc(tSatOut, volWSatOut);
		enthWSatOut = H2O::enth(tSatOut, pMPaOut, WATER);
		enthSSatOut = H2O::enth(tSatOut, pMPaOut, STEAM);
		//      cout << "\n dpdyn " << dpDynSect << " dpstat " << dpStatSect << " pPaOut " << pPaOutSect << " tsat " << tSatOut;
		xOutSect = (enthOut - enthWSatOut) / (enthSSatOut - enthWSatOut);
		if (Base.showDPTubeDetail) {
			prot << "\n 344 xOutSect " << xOutSect;
		}
		/// if after 50 iteration steps no solution is found, there is something wrong and the program terminates
		if (++i > 50) {
			cout << "\n Branch " << this->NbBr << " iterg " << Base.iterg << " Flow " << Flow << " xInSect " << xInSect << " xOutSection " << xOutSect;
			prot << "\n tube " << Number << " dP2phase iteration not converging ";
			prot << "\n Branch " << this->NbBr << " iterg " << Base.iterg << " Flow " << Flow << " xInSect " << xInSect << " xOutSection " << xOutSect;
			switch (region) {
			case TbRegion::Tube: //tube region
				prot << "\n tube region"<<endl;
				break;
			case TbRegion::Orifice: //orifice region
				prot << "\n orifice region"<<endl;
				break;
			case TbRegion::Bend: //tube bend region
				prot << "\n tube bend region";
				prot << "\n one possible reason can be a sharp edged bend (radius = 0) and a Reynolds number close to 40 000 Reynold = " << Reynolds<<endl ;
			}
			cout << "\n tube " << Number << " dPTwoPhase iteration not converging ";
			exit(-1);
		}
		//     double bla = fabs(pPaOutSect / pPaInSect - pPrev / pPaInSect);
	} while (fabs(pPaOutSect / pPaInSect - pPrev / pPaInSect) > 1e-8);

	rhoMeanSect = (rhoInSect + rhoOutSect) / 2.;
	VelSect = MassVel / rhoMeanSect;
	//	prot << "\nreturn from dp2phase ";
/// checking for a phase change, i.e. xOutSect is less than 0 -> condensation or xOutSect is greater than 1 -> only steam leaving
	PhaseChange = false;
	if (xOutSect < 0. || xOutSect > 1.) {
		PhaseChange = true;
		return 1;
	}
	return 0;
} /* dp2Phase */


