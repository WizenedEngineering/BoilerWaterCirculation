/*     ----------------------------------------------------- */
/*     function to for safety criteria of tube               */
/*     ----------------------------------------------------- */
//#include "stdafx.h"

#include "CommonHeader.h"

/**
* @brief function to determine safety factor against different criteria
*
*/
void _tube::Criteria() {
	double  k, k1, k2, hf, hf1,
		xCrit;
	// checking for safety criteria
	double pMPa = pPaOut / 1e6;
	double tSat = H2O::satTemp(pMPa);
	double rhoWSat = 1. / H2O::specVol(tSat, pMPa, WATER);
	double rhoSSat = 1. / H2O::specVol(tSat, pMPa, STEAM);
	double enthWSat = H2O::enth(tSat, pMPa, WATER);
	double enthSSat = H2O::enth(tSat, pMPa, STEAM);
	double deltaEnthEvap = enthSSat - enthWSat;
	SafetyFactor = 100.;
	//   prot <<"pMPa "<<pMPa<< " q " << q << " HeatFlux " << HeatFlux  << " xOut " << xOut << endl;
	/// only heated tubes are checked 
	if (q > 1e-3 && xOut > 1.e-6 && xOut < 1.) {
		/// 1) criterium of Korneev for horizontal tubes
		if (fabs(Height / Length) <= 0.5) { //  horizontal tubes up to 30 degree inclination
			double VelW = MassVel / rhoWSat * (1. - xOut);
			double b1;
			if (pMPa < 4.4129925) {
				b1 = 1.02e-4 * pow(pMPa * 9.80665, 0.25);
			}
			else {
				b1 = 0.22e-4 * pow(pMPa * 9.80665, 0.65);
			}
			double VelWmin = b1 * pow(HeatFlux * 859.8, 0.42) * pow(Dia * 1e3, 0.76);
			//         prot << " VelW " << VelW << " VelWMin " << VelWmin << endl;

			if (fabs(Height / Length) <= 0.1737) {// horizontal tubes up to 10 degree inclination
				if (pMPa > 6.0) {
					VelWmin = VelWmin * 0.55;
				}
				else {
					VelWmin = VelWmin * (1. - 7.5e-3 * pMPa * 10.);
				}
			}
			SafetyFactor = fmin(VelW / VelWmin, SafetyFactor);
		}
		else { // vertical tubes
/**
 * 2) for critical steam quality (x) the criteria are only valid in a certain range of pressure and mass velocity.
 *
 */
			if (pMPa > 0.49 && MassVel > 200. &&
				(Dia >= 4e-3 && Dia <= 32e-3)) {
				/** for mass velocity > 200 and pressure > 4.9 MPa -> Dry-out according Kon'kov */
				if (pMPa <= 2.94) {
					xCrit = pow(HeatFlux, -.125) * 10.795 * pow(MassVel, -1. / 3.) *
						pow(Dia * 1e3, -.07) * exp(pMPa * 10. * .01715);
				}
				else if (pMPa <= 9.8) {
					xCrit = pow(HeatFlux, -.125) * 19.398 * pow(MassVel, -1. / 3.) *
						pow(Dia * 1e3, -.07) * exp(pMPa * 10. * -.00255);
				}
				else {
					xCrit = pow(HeatFlux, -.125) * 32.302 * pow(MassVel, -1. / 3.) *
						pow(Dia * 1e3, -.07) * exp(pMPa * 10. * -.00795);
				}
				SafetyFactor = fmin(SafetyFactor, xCrit / xOut);
			}
			/** for mass velocity > 500 and pressure > 2.9 MPa -> Film boiling according Doroshchuk */
			if (MassVel >= 500. && (pMPa >= 2.9 && pMPa <= 20.) && (Dia >= .004
				&& Dia <= .025)) {
				double pred = pMPa / 22.064;
				double c = (10.3 - (17.5 + pred * 8.) * pred) * 1e3 * sqrt(.008 / Dia);
				xCrit = (log(MassVel / 1e3) * (pred * .68 - .3) - log(HeatFlux * 1e3) + log(c)) /
					(log(MassVel / 1e3) * 1.2 + 1.5);
				if (xCrit > 0. && xCrit < 1.) {
					SafetyFactor = fmin(SafetyFactor, xCrit / xOut);
				}
			}

			/**
			 * in the remaining regions Katto-Ohno and Groeneveld should be used\n
			 * VDI Heatatlas doesn't recommend using Katto-Ohno if xIn > 0\n
			 * in this case Groeneveld should be used
			 */
			double xInlet = (EnthIn - enthWSat) / deltaEnthEvap;
			if (xInlet <= 0.) {
				double relLength = Length / Dia;
				double relSigma = H2O::sigma(tSat) * rhoWSat / (MassVel * MassVel * Length);
				double relRho = rhoSSat / rhoWSat;
				double c = .25;
				if (relLength > 150.) {
					c = .34;
				}
				else if (relLength > 50.) {
					c += (relLength - 50.) * 9e-4;
				}
				hf1 = c * pow(relSigma, .043) / relLength;
				k1 = 1.043 / (c * 4. * pow(relSigma, 0.043));
				k2 = (1. / relLength + .0124) * 5. / (pow(relRho, 0.133) * 6. *
					pow(relSigma, 1. / 3.));
				if (relRho <= .15) {
					double hf2 = pow(relRho, 0.133) * .1 * pow(relSigma, 1. / 3.) /
						(relLength * .0031 + 1.);
					if (hf1 <= hf2) {
						hf = hf1;
					}
					else {
						double hf3 = pow(relRho, 0.133) * .098 * pow(relSigma, 0.433) *
							pow(relLength, 0.27) / (relLength * .0031 + 1.);
						hf = fmin(hf2, hf3);
					}
					k = fmax(k1, k2);
				}
				else {/* rhostern > 0.15 */
					double hf4 = pow(relRho, 0.513) * .234 * pow(relSigma, 0.433) *
						pow(relLength, 0.27) / (relLength * .0031 + 1.);
					if (hf1 <= hf4) {
						hf = hf1;
					}
					else {
						double hf5 = pow(relRho, 0.6) * .0384 * pow(relSigma, 0.173) /
							(pow(relSigma, 0.233) * .28 * relLength + 1.);
						hf = fmax(hf4, hf5);
					}

					if (k1 >= k2) {
						k = k1;
					}
					else {
						double k3 = (pow(relSigma, 0.233) * 1.52 + 1. / relLength) * 1.12 /
							(pow(relRho, 0.6) * pow(relSigma, 0.173));
						k = fmin(k2, k3);
					}
				}
				xCrit = hf * (1. - k * xIn) * 4. * relLength + xIn;
				SafetyFactor = fmin(SafetyFactor, xCrit / xOut);
			}
			/// 3) comparing to critical heat flux according to Groeneveld table
			if (HeatFlux > 1.) {
				//            prot << " chf " << chftable(pMPa, rhoW, rhoS) << " HeatFlux " << HeatFlux  << endl;
				SafetyFactor = fmin(chftable(pMPa, rhoWSat, rhoSSat) / HeatFlux, SafetyFactor);
			}
		}
	}
}
//