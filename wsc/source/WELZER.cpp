//#include "stdafx.h"
#include "CommonHeader.h"

double _tube::dpdyn_Welzer(double zeta, double volW, double volS, double enthW, double enthS,
	double visW, double xInlet, double xOutlet) {
	/* System generated locals */
	double d__1, d__2;

	/* Local variables */
	double phiIn, phiOut, phiMean, beta1, beta2;
	double ksiTubeW = 0.;
	if (Base.showDPTube) {
		prot << "\n welzer massvel " << MassVel << " xIn " << xIn << " xOut " << xOut << " zeta " << zeta;
	}

	double RatioDensity = volS / volW;
	double xStar = sqrt(RatioDensity) / (sqrt(RatioDensity) + 1.); // transition point from slip to homogeneous
	d__1 = (RatioDensity - 1.) / (1000. * volS);
	double bStar = sqrt(d__1 * d__1 * d__1);
	/* Computing 2nd power */
	d__1 = 1. - bStar * 6. / (Dia * 1e3);
	beta1 = 1. / (d__1 * d__1);
	/* Computing 5th power */
	d__2 = d__1;
	d__1 *= d__1;
	beta2 = 1. / (d__2 * (d__1 * d__1));
	double AngleFactor = fmin(1., 2. * HeightSection / LengthSection);
	AngleFactor = fmax(0., AngleFactor);
	if (HeatSection > .001) {
		double enthEvap = enthS - enthW;
		d__1 = volS * volS * 14200. * Dia * MassVel * (HeatFlux * HeatFlux) /
			(bStar * 9.80665 * visW * (enthEvap * enthEvap));
		//14200 = 14.2*1e3; velW = MassVel*volW, kinematic Viscosity = dyn.Viscosity(i.e.visW)*volW, volW can be crossed out
		double lw = log10(d__1) * 2.;
		lw = 1. / (lw * lw);
		/* Computing 2nd power */
		d__1 = 1. - (xIn + xOut) / 2.;
		double lf = zeta + d__1 * d__1 * (lw * beta2 - zeta);
		if (HeightSection / LengthSection < .001) { //horizontal or down
			phiMean = (xInlet + xOutlet) / 2. * (RatioDensity - 1.) + 1.;
			ksiTubeW = lf * phiMean * LengthSection / Dia;
		}
		else {
			if (xInlet < xStar && xOutlet < xStar) {
				phiIn = xInlet * (RatioDensity - 1.) / (AngleFactor *
					xInlet * (RatioDensity - 1.) / sqrt(RatioDensity) + 1.);
				phiOut = xOutlet * (RatioDensity - 1.) / (AngleFactor *
					xOutlet * (RatioDensity - 1.) / sqrt(RatioDensity) + 1.);
				phiMean = (phiOut + phiIn) / 2. + 1.;
				ksiTubeW = phiMean * lw * beta2 * LengthSection / Dia;
			}
			else if (xInlet >= xStar && xOutlet >= xStar) {
				phiMean = (xInlet + xOutlet) / 2. * (RatioDensity - 1.) + 1.;
				ksiTubeW = zeta / RatioDensity * (phiMean * phiMean) * LengthSection / Dia;
			}
			else if (xInlet < xStar && xOutlet > xStar) {
				double lstar = LengthSection * (xStar - xInlet) / (xOutlet - xInlet);
				phiIn = xInlet * (RatioDensity - 1.) /
					(AngleFactor * xInlet * (RatioDensity - 1.) / sqrt(RatioDensity) + 1.);
				double phistar1 = xStar * (RatioDensity - 1.) /
					(AngleFactor * xStar * (RatioDensity - 1.) / sqrt(RatioDensity) + 1.);
				double phiMean1 = (phistar1 + phiIn) / 2. + 1.;
				double phiMean2 = (xStar + xOutlet) / 2. * (RatioDensity - 1.) + 1.;
				ksiTubeW = phiMean1 * lw * beta2 * lstar / Dia + zeta /
					RatioDensity * phiMean2 * phiMean2 * (LengthSection - lstar) / Dia;
			}
		}
	}
	else {
		if (HeightSection / LengthSection < .001) { //horizontal or down
			phiMean = xInlet * (RatioDensity - 1.) + 1.;
			ksiTubeW = zeta * phiMean * LengthSection / Dia;
		}
		else {
			if (xInlet < xStar && xOutlet < xStar) {
				phiMean = xInlet * (RatioDensity - 1.) / (AngleFactor * xInlet * (RatioDensity - 1.) / sqrt(RatioDensity) + 1.) + 1.;
				ksiTubeW = zeta * phiMean * phiMean * LengthSection / Dia;
			}
			else {
				phiMean = xInlet * (RatioDensity - 1.) + 1.;
				ksiTubeW = zeta / RatioDensity * phiMean * phiMean * LengthSection / Dia;
			}
			//         beta1 = 1.;
		}
	}
	double dp = ksiTubeW * volW / 2. * MassVel * MassVel;
	if (Base.showDPTube) {
		prot << "\n dp " << dp;
	}
	return dp;
} /* welzer_ */


double _tube::Density_Welzer(double x, double volW, double volS, double& VoidFraction) {
	/*    this function is based on Welzer */
	/*    to allow for tube inclination other than vertical following assumptions are used */
	/*    1. Above 30deg to horizontal the density is like in vertical tubes */
	/*    2. in horizontal tubes no slip (homogeneous) */
	/*    3. for continuity between 30deg and horizontal a factor 2* height/length is used */
	/*    4. downward flown tubes no slip (homogeneous)  */
	double phi;
	double RatioDensity = volS / volW;

	double xStar = sqrt(RatioDensity) / (sqrt(RatioDensity) + 1.); // transition point from slip to homogeneous
	if (x < xStar) {
		double AngleFactor = fmin(1., 2. * HeightSection / LengthSection);
		AngleFactor = fmax(0., AngleFactor);
		phi = (x * (RatioDensity - 1.) + 1.) /
			(AngleFactor * x * (RatioDensity - 1.) / sqrt(RatioDensity) + 1.);
	}
	else {
		phi = x * (RatioDensity - 1.) + 1.;
	}

	VoidFraction = RatioDensity * (1. - 1. / phi) / (RatioDensity - 1);
	if (Base.showDPTube) {
		prot << "\n x " << x << " xstar " << xStar << " Phi " << phi
			<< "VoidFraction " << VoidFraction
			<< " rho " << 1. / (volW * phi) <<
			" rho hom " << 1. / (x * volS + (1. - x) * volW) <<
			" rhoW " << 1. / volW << " rhoS " << 1. / volS;
	}
	return 1. / (volW * phi);
} /* welzer_ */


