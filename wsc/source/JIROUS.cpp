// source: F.Jirous, Analytische Methode der Berechnung des Naturumlaufs
// bei Dampferzeugern, VGB Kraftwerkstechnik 58, Heft 5 


#include "CommonHeader.h"

double _tube::dpdyn_Jirous(double zeta, double pMPa, double volSIn, double volSOut, double volWIn, double volWOut, double x) {
	double RatioDensity, pi1, tpf, R1, dp;

	RatioDensity = (volSIn + volSOut) / (volWIn + volWOut);
	pi1 = RatioDensity - 1.;
	tpf = (1.58 - 0.0213 * pMPa - 2.25956161e-4 * pMPa * pMPa) * pi1;

	if (x < 0.8) {
		R1 = 1. + tpf * x;
	}
	else {
		double d__1 = ((x - .8) * (1. + tpf - RatioDensity) - .02 * tpf);
		R1 = 1. - d__1 * d__1 /
			(.04 * (1.0 + tpf - RatioDensity)) +
			tpf * ((0.81 * tpf - 0.8 * (RatioDensity - 1)) /
				(1. + tpf - RatioDensity));
	}
	dp = R1 * zeta * LengthSection / Dia * (volWIn + volWOut) / 4. * MassVel * MassVel;
	//   prot << "dpdyn jiroush " << dp << " x "<<x << " R1 " <<R1<<endl;
	return dp;
}

//   jirousdensity

double _tube::Density_Jirous(double x, double volW, double volS, double& VoidFraction) {
	/*    this function is based on Jirous article*/
	/*    to allow for tube inclination other than vertical following assumptions are used */
	/*    1. Above 30deg to horizontal the density is like in vertical tubes */
	/*    2. in horizontal tubes no slip (homogeneous) */
	/*    3. for continuity between 30deg and horizontal a factor 2*height/length is used */
	/*    4. downward flown tubes no slip (homogeneous)  */

	double RatioDensity = volS / volW;

	double slip = pow(10., (0.031636 + (0.07773 + .051306 * log10(RatioDensity)) *
		log10(RatioDensity)));

	//static pressure difference
	if (HeightSection < 1e-6) {
		slip = 1.;
	}
	else if (HeightSection / LengthSection < 0.5) {
		slip = (slip - 1.) * 2. * HeightSection / LengthSection + 1.;
	}
	VoidFraction = x / (slip / RatioDensity + x * (1. - slip / RatioDensity));

	double Density = (1. - VoidFraction) / volW + VoidFraction / volS;
	//   prot << "dens jiroush " << Dens << " x " << x << " pi " << pi << " slip " << slip << endl;

	return Density;
}

