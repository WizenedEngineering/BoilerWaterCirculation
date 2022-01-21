/* Determines the safety factor against wavy / stratified flow in horizontal tubes
* according the criterion by Taitel - Dukler
*
*Source:  A Model for Predicting Flow Regime Transitions in Horizontaland Near Horizontal Gas - Liquid Flow, Taitel, Y, Dukler, A.E., AIChE Journal(Vol.22, No 1)
*
* param[in] pPa Pressure[Pa]
* param[in] quality steam mass content, steam quality[-]
* author Rainer_Jordan@<very, very warm>mail.com
* Licence
* Licensed under the European Union Public Licence(EUPL), Version 1.2 or later
* date January 2022
*/
#include "stdafx.h"

#include "CommonHeader.h"

double _tube::Taitel_Dukler(
	double pPa, // Pressure in Pa
	double quality // steam quality, steam mass content
) {
	// checking for flow pattern, wavy or stratified flow
	// according to: A Model for Predicting Flow Regime Transitions in Horizontal and
	// Near Horizontal Gas-Liquid Flow, Taitel,Y,Dukler,A.E. , AIChE Journal (Vol.22, No 1)
	/// keeping nomenclature close to Literature\n
	/// index G=gas -> steam\n
	/// index L=liquid -> water\n

	//sin and cos of angle to horizontal
	//alpha positive for downward inclination
	double sinAlpha = -Height / Length;
	if (fabs(sinAlpha) > 0.17365 ||//not horizontal (> 10 deg)
		quality < 1e-6 ||//only liquid
		RadiusBend > 1e-3) { // by secondary flow no separation in bends/elbows likely
		return 100.;
	}

	double cosAlpha = sqrt(1. - sinAlpha * sinAlpha);

	double pMPa = pPa / 1e6;
	double tSat = H2O::satTemp(pMPa);
	double rhoL = 1. / H2O::specVol(tSat, pMPa, WATER);
	double rhoG = 1. / H2O::specVol(tSat, pMPa, STEAM);
	double dynViscG = H2O::dynVisc(tSat, 1. / rhoG);
	double dynViscL = H2O::dynVisc(tSat, 1. / rhoL);
	//			double SurfTens = H2O::sigma(tSat) / 1000.;

	double diff,
		CL = 0.046, //turbulent 
		CG = 0.046, //turbulent
		exp_m = 0.2, //turbulent
		exp_n = 0.2; //turbulent
	double Xsquare;
	double Y;
	double DhydL, DhydG, DLtilde, DGtilde;
	double uLtilde, uGtilde, hLtilde, hLu, hLl, AGtilde, ALtilde, ReG, ReL, fact, SLtilde, SGtilde, Sitilde;
	double uLs = (1. - quality) * MassVel / rhoL;
	double uGs = quality * MassVel / rhoG;
	hLu = 1.;
	hLl = 0.;
	for (int i = 0; i <= 20; i++) {
		hLtilde = (hLl + hLu) / 2.;
		fact = (2. * hLtilde - 1.);
		SGtilde = acos(fact);
		SLtilde = M_PI - SGtilde;
		Sitilde = sqrt(1. - fact * fact);
		ALtilde = 0.25 * (SLtilde + fact * Sitilde);
		AGtilde = 0.25 * (SGtilde - fact * Sitilde);
		uLtilde = area / Dia / Dia / ALtilde;
		uGtilde = area / Dia / Dia / AGtilde;
		DLtilde = 4. * ALtilde / SLtilde;
		DGtilde = 4. * AGtilde / (SLtilde + Sitilde);
		DhydL = DLtilde * Dia;
		DhydG = DGtilde * Dia;
		ReL = MassVel * (1. - quality) * DhydL / dynViscL;
		if (ReL > 2300.) {
			CL = 0.046, //turbulent 
			exp_n = 0.2; //turbulent
		}
		else {
			CL = 16.; //laminar
			exp_n = 1.; //laminar
		}
		ReG = MassVel * quality * DhydG / dynViscG;
		if (ReG > 2300.) {
			CG = 0.046; //turbulent
			exp_m = 0.2; //turbulent
		}
		else {
			CG = 16.; //laminar
			exp_m = 1.;//laminar
		}

		Xsquare = (CL * pow(uLs * Dia * rhoL / dynViscL, -exp_n) * rhoL * uLs * uLs) /
			(CG * pow(uGs * Dia * rhoG / dynViscG, -exp_m) * rhoG * uGs * uGs);

		Y = (rhoL - rhoG) * 9.8066 * sinAlpha /
			((4. * CG / Dia * pow(uGs * Dia * rhoG / dynViscG, -exp_m) * rhoG * uGs * uGs / 2.));

		diff = Xsquare * (pow(uLtilde * DLtilde, -exp_n) * uLtilde * uLtilde * SLtilde / ALtilde) -
			(pow(uGtilde * DGtilde, -exp_m) * uGtilde * uGtilde * (SGtilde / AGtilde + Sitilde / ALtilde + Sitilde / AGtilde))
			- 4. * Y;
		//				 cout << "\n hL " << hL << " diff " << diff << " CL " << CL << " CG " << CG<<" AL "<<AL<<" AG "<<AG<<" sum "<< (AL+AG) / (M_PI/4.);
		if (fabs(diff) < 1e-6) {
			break;
		}
		if (diff > 0.) {
			hLl = hLtilde;
		}
		else {
			hLu = hLtilde;
		}
		if (i == 50) {
			cout << "\n  taitel iteration not converging" << endl;
		}
	}
	double C2 = 1. - hLtilde;
	double FroudeSquare = rhoG / (rhoL - rhoG) * uGs * uGs / (Dia * 9.8066 * cosAlpha);
	return FroudeSquare * (uGtilde * Sitilde / AGtilde / C2 / C2); // should be above 1 to avoid separation
}

