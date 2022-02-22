/*
* function to determine safety factor against different criteria
*
* VDI Heat Atlas divides the range of mass velociy and pressure into 4 regions and recommends a method for each region
* Those methods are:
* Korneev for horizontal/inclined heated tubes
* Taitel-Dukler and Steiner for flow separation (only critical in heated tubes)
* Kon'kov for critical steam quality
* Doroshchuk for critical steam quality
* Katto-Ohno for critical steam quality
* Groeneveld for critical heat flux
* Kon'kov has following limitation: Dia < 33mm, MassVel > 200 kg/m2s; outside the limitation it will be disregarded
* Katto-Ohno is only valid for xIn <= 0
*
* the minimum safety factor is used
* Author: Rainer_Jordan@<very, very warm>mail.com
* Licence: 
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* date: September 2021
*/

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
	SafetyFactor = 10.;
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
			SafetyKorneev = VelW / VelWmin;
			SafetyFactor = fmin(SafetyKorneev, SafetyFactor);

// in horizontal heated tubes (angle to horizontal < 10 deg) the flow should not be stratified/wavy or mist
			// Taitel-Dukler
			if (fabs(Height / Length) < 0.17365) {//horizontal (< 10 deg)
				SafetyTaitel_Dukler = Taitel_Dukler(pPaOut, xOut);
				SafetyFactor = fmin(SafetyTaitel_Dukler, SafetyFactor);

				// Steiner
				SafetySteiner = Steiner(pPaOut, xOut);
				SafetyFactor = fmin(SafetySteiner, SafetyFactor);
			}
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
				SafetyKonkov = xCrit / xOut;
				SafetyFactor = fmin(SafetyFactor, SafetyKonkov);
			}
			/** for mass velocity > 500 and pressure > 2.9 MPa -> Film boiling according Doroshchuk */
			if (MassVel >= 500. && (pMPa >= 2.9 && pMPa <= 20.) && (Dia >= .004
				&& Dia <= .025)) {
				double pred = pMPa / 22.064;
				double c = (10.3 - (17.5 + pred * 8.) * pred) * 1e3 * sqrt(.008 / Dia);
				xCrit = (log(MassVel / 1e3) * (pred * .68 - .3) - log(HeatFlux * 1e3) + log(c)) /
					(log(MassVel / 1e3) * 1.2 + 1.5);
				if (xCrit > 0. && xCrit < 1.) {
					SafetyDoroshchuk = xCrit / xOut;
					SafetyFactor = fmin(SafetyFactor, SafetyDoroshchuk);
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
				SafetyKatto_Ohno = xCrit / xOut;
				SafetyFactor = fmin(SafetyFactor,SafetyKatto_Ohno);
			}
			/// 3) comparing to critical heat flux according to Groeneveld table
			if (HeatFlux > 1.) {
				//            prot << " chf " << chftable(pMPa, rhoW, rhoS) << " HeatFlux " << HeatFlux  << endl;
				SafetyGroeneveld = chftable(pMPa, rhoWSat, rhoSSat) / HeatFlux;
				SafetyFactor = fmin(SafetyGroeneveld, SafetyFactor);
			}
		}
	}
}
//
