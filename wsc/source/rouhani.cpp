//#include "stdafx.h"
#undef MAINFUNCTION
// It is assumed that the tube sections are small enough that the properties
// are mean value between inlet and outlet
#include "CommonHeader.h"

//dpdyn according Becker

double _tube::dpdyn_Becker(double zeta, double xInSection, double xOutSection, double pPaInSection, double pPaOutSection,
	double volWSatIn, double volWSatOut) {
	/*  dynamic pressure drop */
	/*  turbulent flow */
	double dp, r1, xpIn, xpOut;

	if (fabs(xOutSection - xInSection) <= 1e-6) {
		r1 = 1. + 2400. * pow((xOutSection + xInSection) / ((pPaInSection + pPaOutSection) / 98066.5), 0.96);
	}
	else {
		xpIn = xInSection * 98066.5 / pPaInSection;
		xpOut = xOutSection * 98066.5 / pPaOutSection;
		r1 = 1. + 1224.49 * (pow(xpOut, 1.96) - pow(xpIn, 1.96))
			/ (xpOut - xpIn);
	}
	dp = r1 * zeta * LengthSection / Dia * (volWSatIn + volWSatOut) / 4. * MassVel * MassVel;
	return dp;
}

double _tube::Gomez(double x, double rhoW, double rhoS, double SurfTens, double VoidFractionInput,
	double HeightRatio, double VoidFractionHomogeneous, double UsG) {
	double C0 = 1.15;
	double vr = 1.53 * sqrt(sqrt(SurfTens * 9.80665 * (rhoW - rhoS)) / rhoW) * (1. - VoidFractionInput) * HeightRatio;
	return 1. / (C0 / VoidFractionHomogeneous + vr / UsG);//VoidFraction
}

double _tube::Density_Rouhani(double x, double rhoW, double rhoS, double SurfTens, double& VoidFraction) {
	/* Local variables */
	double C0, vr, AngleFactor, HeightRatio, UsG;
	/** the equations from Rouhani are used for upward flow /
	 *    to allow for tube inclination other than vertical following assumptions are used
	 *    1. Above 30deg to horizontal the density is like in vertical tubes
	 *    2. in horizontal tubes Steiner (Heat Atlas, H3.1, eq.26) recommends using eq. 4 of Rouhani's report
	 *    3. for continuity between 30deg and horizontal a factor 2*height/length is used and applied to the factor in C0\n
	 *       for horizontal tubes the factor is 0.12 and for vertical tubes 0.2
	 *    4. downward flown tubes formula from Gomez */
	HeightRatio = HeightSection / LengthSection;
	UsG = x * MassVel / rhoS; //superficial gas (steam) velocity
	double VoidFractionHomogeneous = rhoW * x / (rhoW * x + rhoS * (1. - x)); // homogeneous

	if (HeightSection >= 0.) {
		AngleFactor = fmin(1., 2. * HeightRatio);
		C0 = sqrt(sqrt(Dia * 9.80665) * rhoW / MassVel) * (AngleFactor * 0.08 + 0.12) * (1. - x) + 1.;
		vr = (1. - x) * 1.18 * sqrt(sqrt(SurfTens * 9.80665 * (rhoW - rhoS)) / rhoW);
		VoidFraction = 1. / (C0 / VoidFractionHomogeneous + vr / UsG);
		/**
		 * Upward flow: the range for void fraction is between homogeneous void fraction and x\n
		 * in some cases like low mass velocity or high x content the result can lay outside this range \n
		 * that means we are outside the validity of this formula\n
		 * for high x content a homogeneous solution is sensible (the flow pattern can change to mist flow which is more like homogeneous ) \n
		 * for low mass velocity the steam velocity can be significantly higher than water velocity. The void fraction approaching x\n
		 */
		if (VoidFraction >= .25) {
			C0 = (1. - x) * 0.2 + 1.;
			VoidFraction = 1. / (C0 / VoidFractionHomogeneous + vr / UsG);
		}
		VoidFraction = fmax(VoidFraction, x);
		VoidFraction = fmin(VoidFraction, VoidFractionHomogeneous);
	}
	else {
		/**
		 * Downward flow: the range for void fraction is between x and 1.\n
		 * in some cases like low mass velocity or high x content the result can lay outside this range \n
		 * that means we are outside the validity of this formula\n
		 * for high x content a homogeneous solution is sensible (the flow pattern can change to mist flow and mist flow is more like homogeneous ) \n
		 * for low mass velocity there is a high chance for finely dispersed bubbles that again lead to no velocity difference of the phases -> homogeneous\n
		 */
		double epsLower = x;
		double VoidFractionInput;
		double epsHigher = 0.9999999;
		for (int i = 0; i < 50; i++) {
			VoidFractionInput = (epsLower + epsHigher) / 2.;
			VoidFraction = Gomez(x, rhoW, rhoS, SurfTens, VoidFractionInput,
				HeightRatio, VoidFractionHomogeneous, UsG);
			//            prot << "\n i " << i << " start " << VoidFractionInput << " void " << VoidFraction<<" diff "<<VoidFractionInput - VoidFraction;
			if (fabs(VoidFraction - VoidFractionInput) < 1e-6) break;
			if ((VoidFractionInput - VoidFraction) > 0.) {
				epsHigher = VoidFractionInput;
			}
			else {
				epsLower = VoidFractionInput;
			}
		}
	}
	if (VoidFraction < x || VoidFraction > 0.9999999) {
		VoidFraction = VoidFractionHomogeneous;
	}

	return rhoW * (1. - VoidFraction) + rhoS * VoidFraction;
} /* Density_Rouhani */

