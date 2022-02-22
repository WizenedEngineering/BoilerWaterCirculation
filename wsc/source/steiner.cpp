/* Determines the flow pattern in horizontal tubes
* according the flow chart by Steiner
*
* Source: VDI Heat Atlas, Second Edition, H3.1
*
*param[in] pPa Pressure[Pa]
*param[in] quality steam mass content, steam quality[-]
*Author: Rainer_Jordan@<very, very warm>mail.com
*Licence:
* Licensed under the European Union Public Licence(EUPL), Version 1.2 or later
* date: September 2021
*/

//#include "stdafx.h"

#include "CommonHeader.h"

double _tube::Steiner(//return safety factor against flow separation FrGm / FrGmlim 
	double pPa, // Pressure in Pa
	double quality // steam quality, steam mass content
) {
	// checking for flow pattern
	// according to Heat Atlas, H3.1 MassVel Patterns in Evaporator Tubes
	/// keeping nomenclature close to Literature\n
	/// index G=gas -> steam\n
	/// index L=liquid -> water\n

	//sin and cos of angle to horizontal
	double sinTheta = Height / Length;
	if (fabs(sinTheta) > 0.17365 ||//not horizontal (> 10 deg)
		quality < 1e-6 ||//only liquid
      RadiusBend > 1e-3) { // by secondary flow no separation in bends/elbows likely
		steiner = FlowPattern::undetermined;
	}
	else {
		double cosTheta = sqrt(1. - sinTheta * sinTheta);

		double pMPa = pPa / 1e6;
		double tSat = H2O::satTemp(pMPa);
		double rhoL = 1. / H2O::specVol(tSat, pMPa, WATER);
		double rhoG = 1. / H2O::specVol(tSat, pMPa, STEAM);
		double dynVisG = H2O::dynVisc(tSat, 1. / rhoG);
		double dynVisL = H2O::dynVisc(tSat, 1. / rhoL);
		double SurfTens = H2O::sigma(tSat) / 1000.;

		//calculation of parameters
		//MArtinelli parameter
//		prot << "32 steiner  quality " << quality << endl;

		double X = pow((1. - quality) / quality, 0.875) * sqrt(rhoG / rhoL) * pow(dynVisL / dynVisG, 0.125);

		double ReLFrG = ((MassVel * MassVel * MassVel * quality * quality * (1. - quality)) / (rhoG * (rhoL - rhoG) * dynVisL * 9.8066 * cosTheta));
		double FrGm = (MassVel * MassVel * quality * quality / (9.8066 * Dia * rhoL * rhoG));
		double ReL = MassVel * (1. - quality) * Dia / dynVisL;
		double ksiL = 0.3164 / pow(ReL, 0.25);
		double FrEuL = (ksiL * MassVel * MassVel * (1. - quality) * (1. - quality) / (2. * Dia * rhoL * (rhoL - rhoG) * 9.8066 * cosTheta));
		double WeFrL = 9.8066 * Dia * Dia * rhoL / SurfTens;

		//        prot<<"X "<<X<<"ReLFrg "<<ReLFrG<<" FrGm "<<FrGm<<" FrEuL "<<FrEuL<<" WeFrL "<<WeFrL<<endl;
				// determination of limit curves
				// this part not needed because of Bisection with starting value 0.5
				// un-wetted angle of tube cross section for stratified flow phi in radians
		double phi, phiu, phil, diff;
		double Ui, UL, UG, hL, hLu, hLl, AG, AL, psi, ReG, ksiG, psiu, psil,
			FrEuG_inv;//inverse of FrEuG used, if sinTheta is 0 FrEuG is inf , inverse -> 0

		hLu = 1.;
		hLl = 0.;
		ReG = MassVel * quality * Dia / dynVisG;
		ksiG = 0.3164 / pow(ReG, 0.25);
		FrEuG_inv = (2. * Dia * 9.8066 * rhoG * (rhoL - rhoG) * sinTheta) / (ksiG * MassVel * MassVel * quality * quality);
		//      prot << "\n quality " << quality << " frEuG_inv " << FrEuG_inv <<" X "<<X;
		bool fail = false;
		for (int i = 1; i <= 50; i++) {
			hL = (hLl + hLu) / 2.;
			Ui = 2. * sqrt(hL * (1. - hL));
			if (hL <= 0.5) {
				psi = 2. * asin(Ui);
				UL = psi / 2.;
				UG = M_PI - UL;
				AL = (psi - sin(psi)) / 8.;
				AG = M_PI_4 - AL;
			}
			else {
				phi = 2. * asin(Ui);
				UG = phi / 2;
				UL = M_PI - UG;
				AG = (phi - sin(phi)) / 8.;
				AL = M_PI_4 - AG;
			}
			double arg = (pow((UG + Ui) / M_PI, 0.25) *
				M_PI * M_PI / (64. * AG * AG) *
				((UG + Ui) / AG + Ui / AL) - FrEuG_inv);
			//        prot << "\n arg" << arg;
			if (arg < 0.) {
				fail = true;
				break;
			}

			diff = X - sqrt(arg * pow(M_PI / UL, 0.25) *
				64. * AL * AL * AL / (M_PI * M_PI * UL));
			//         prot << "\n hl " << hL << " AL " << AL << " AG " << AG << " diff " << diff << endl;
			//         prot << "\n UG " << UG << " Ui " << Ui <<" UL "<<UL;
			if (fabs(diff) < 1e-6) {
				break;
			}
			if (diff > 0.) {
				hLl = hL;
			}
			else {
				hLu = hL;
			}
			if (i == 50) {
				prot << "\n  Steiner iteration not converging" << endl;
				fail = true;
			}
		}

		if (fail) { // iteration by Martinelli parameter failed probably caused by low pressure  
		//  Steiner based the steam volume content on Rouhani
			double epsilon = quality / rhoG / ((1. + 0.12 * (1. - quality)) * (quality / rhoG + (1. - quality) / rhoL) + 1.18 * (1. - quality) * pow(9.8066 * SurfTens * (rhoL - rhoG), 0.25) / (MassVel * sqrt(rhoL))); //void fraction
			if (epsilon == 0.5) {
				hL = 0.5;
				Ui = 1.;
				UL = M_PI / 2.;
				UG = UL;
				AL = M_PI / 8.;
				AG = AL;
			}
			else if (epsilon < 0.5) {
				//     hl above centerline
				phiu = M_PI;
				phil = 0.;
				for (int i = 1; i <= 50; i++) {
					phi = (phiu + phil) / 2.;
					diff = 2. * M_PI * epsilon + sin(phi) - phi;
					if (fabs(diff) < 1e-6) {
						break;
					}
					if (diff < 0.) {
						phiu = phi;
					}
					else {
						phil = phi;
					}
				}
				hL = 15. * M_PI * (1. - epsilon) / (8. * (3 * sin(phi / 2.) + 4. * sin(phi / 4.)));
				Ui = 2. * sqrt(hL * (1. - hL));
				phi = 2. * asin(Ui);
				UG = phi / 2;
				UL = M_PI - UG;
				AG = (phi - sin(phi)) / 8.;
				AL = M_PI_4 - AG;
			}
			else { //hl below centerline
				psiu = 2. * M_PI;
				psil = 0.;
				for (int i = 1; i <= 50; i++) {
					psi = (psiu + psil) / 2.;
					diff = 2. * M_PI * epsilon + sin(psi) - psi;
					if (fabs(diff) < 1e-6) {
						break;
					}
					if (diff < 0.) {
						psiu = psi;
					}
					else {
						psil = psi;
					}
				}
				hL = 15. * M_PI * (1. - epsilon) / (8. * (3 * sin(psi / 2.) + 4. * sin(psi / 4.)));
				Ui = 2. * sqrt(hL * (1. - hL));
				psi = 2. * asin(Ui);
				UL = psi / 2.;
				UG = M_PI - UL;
				AL = (psi - sin(psi)) / 8.;
				AG = M_PI_4 - AL;
			}
		}
		//		prot << " hl " << hL << " AL "<< AL << " AG "<<AG <<     endl;
				// limits based on Martinelli parameter
		double ReLFrGlim = (226.3 * 226.3) / M_PI / M_PI / M_PI * AL * AG * AG;
		//		prot << " ReLFrG " << ReLFrG << " ReLFrGlim " << ReLFrGlim << endl;
		if (ReLFrG <= ReLFrGlim) {
			//             prot << "stratified" << endl;
			steiner = FlowPattern::stratified; //stratified flow
			return 0.;
		}
		double FrGmlim = 16. * AG * AG * AG / (M_PI * M_PI * sqrt(1. - (2. * hL - 1.) * (2. * hL - 1.))) *
			(M_PI * M_PI / (25. * hL * hL) / WeFrL + 1. / cosTheta);
		//	prot << " FrGm " << FrGm << " FrGmlim " << FrGmlim << endl;
		double Safety = FrGm / FrGmlim;
		if (Safety <= 1.) {
			//          prot << "wavy" << endl;
			steiner = FlowPattern::wavy; //wavy flow
			return Safety;
		}
		double FrEuLlim = 128. * AG * AL * AL / (M_PI * M_PI * Ui);
		//		prot << " FrEuL " << FrEuL << " FrEuLlim " << FrEuLlim << endl;
		if (FrEuL >= FrEuLlim) {
			//               prot << "bubble" << endl;
			steiner = FlowPattern::bubble; // bubble flow
			return Safety;
		}
		if ((X >= 0.34 && ReL >= 1187. && ReG >= 1187. && FrGm > FrGmlim) ||
			(X >= 0.51 && ReL < 1187. && ReG >= 1187. && FrGm >= FrGmlim)) {
			//                   prot << "plug" << endl;
			steiner = FlowPattern::plugSlug; //plug or slug flow
			return Safety;
		}
		double zetaPh = (1.138 + 2. * log(M_PI / (1.5 * AL)));
		zetaPh = 1. / zetaPh / zetaPh;
		double FrGmlim2 = 7680. * AG * AG / (M_PI * M_PI * zetaPh) / WeFrL;
		//		prot << " FrGm " << FrGm << " FrGmlim2 " << FrGmlim2 << endl;
		if (X < 0.51 && FrGm >= FrGmlim2) {
			//              prot << "mist" << endl;
			steiner = FlowPattern::mist; //mist flow
			return 0.;
		}
		if (X < 0.51 && FrGm > FrGmlim && FrGm < FrGmlim2) {
			//               prot << "annular" << endl;
			steiner = FlowPattern::annular; //annular flow
			return Safety;
		}
	}
	return 10.;
}
