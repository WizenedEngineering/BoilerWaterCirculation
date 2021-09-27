//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

// It is assumed that the tube sections are small enough that the properties
// are mean value between inlet and outlet

double _tube::dpdyn_HeatAtlas(double pPaInSect, double pPaOutSect, double volSSatIn,
	double volSSatOut, double volWSatIn, double volWSatOut, double dynVisWSatIn,
	double dynVisWSatOut, double dynVisSSatIn, double dynVisSSatOut, double xMean
) {
	double Reynolds,
		zeta,
		FrictCoeffTube = 0.,
		dpdynSect = 0.,
		r1;

	if (Base.showDPTubeDetail) {
		prot << "\n pdyn according Vdi Heat Atlas  xMean " << xMean;
	}
	/**
	 * For horizontal tubes the calculation according Garcia (VDI Heat Atlas, H2.2. equations 15 -24) is used.
	 */
	double sinTheta = HeightSection / LengthSection;
	if (fabs(sinTheta) <= 0.17365) {// horizontal (< 10 deg)
		/**
		* the flow pattern is determined according Steiner
		*/
		Steiner((pPaInSect + pPaOutSect) / 2., xMean);
		/**
		 * mist flow is calculated as steam flow (only steam portion regarded).
		 */
		if (steiner == FlowPattern::mist) { //mist flow (calculated as steam flow, only steam portion regarded)
			Reynolds = Dia * xMean * MassVel / (dynVisWSatIn + dynVisWSatOut) * 2.;
			zeta = FrictFact(Reynolds, Base.Rough / Dia);
			FrictCoeffTube = zeta * LengthSection / Dia;
			if (Base.showDPTubeDetail) {
				prot << "\n steiner = mist  reynolds " << Reynolds << " zeta " << zeta << " frictcoeff " << FrictCoeffTube;
			}
			dpdynSect = FrictCoeffTube * (volSSatIn + volSSatOut) / 4. * MassVel * MassVel;
		}
		else {
			double a1, a2, b1, b2, c, d, t, wL, wG, wM, rhoM, lambda, fm;
			switch (steiner) {
			case FlowPattern::stratified:
				//				Stratified flow
			case FlowPattern::wavy:
				// wavy flow (regarded as stratified)
				if (Base.showDPTubeDetail) {
					prot << "\n steiner: stratified or wavy ";
				}

				a1 = 13.98;
				b1 = -0.9501;
				a2 = 0.0445;
				b2 = -0.1874;
				c = 9.275;
				d = 0.0324;
				t = 300.;
				break;
			case FlowPattern::bubble:
				//				Disperse bubble flow
				if (Base.showDPTubeDetail) {
					prot << "\n steiner: Disperse bubble flow";
				}
				a1 = 13.98;
				b1 = -0.9501;
				a2 = 0.1067;
				b2 = -0.2629;
				c = 2.948;
				d = 0.2236;
				t = 304.;
				break;
			case FlowPattern::plugSlug:
				//			Slug flow
				if (Base.showDPTubeDetail) {
					prot << "\n steiner: Slug flow";
				}
				a1 = 13.98;
				b1 = -0.9501;
				a2 = 0.1067;
				b2 = -0.2629;
				c = 3.577;
				d = 0.2029;
				t = 293.;
				break;
				/**
				 * for annular flow the set of parameters result in a factor that is by 10e6 higher than for unknown flow pattern\n
				 * this difference is too high and the parameter set for unknown flow pattern is used
				 */
				 //                  case FlowPattern::annular:
				 //                     //				Annular flow
				 //                     a1 = 3.671;      fm is by factor 10e6 higher than others
				 //                     b1 = 0.6257;
				 //                     a2 = 0.0270;
				 //                     b2 = -0.1225;
				 //                     c = 2.191;
				 //                     d = 0.2072;
				 //                     t = 10000.;
				 //                     break;
			default:
				if (Base.showDPTubeDetail) {
					prot << "\n case: no distinct flow pattern ";
				}
				a1 = 13.98;
				b1 = -0.9501;
				a2 = 0.0925;
				b2 = -0.2534;
				c = 4.864;
				d = 0.1972;
				t = 293.;
			}

			wL = (1. - xMean) * MassVel * (volWSatIn + volWSatOut) / 2.;
			wG = xMean * MassVel * (volSSatIn + volSSatOut) / 2.;
			wM = wL + wG;

			lambda = wL / wM;
			rhoM = 2. / (volWSatIn + volWSatOut) * lambda + 2. / (volSSatIn + volSSatOut) * (1. - lambda);
			if (Base.showDPTubeDetail) {
				prot << "\n 120 wL " << wL << " wg " << wG << " wm " << wM << " lambda " << lambda << " RhoM " << rhoM;
			}

			Reynolds = wM * Dia * 4. / (volWSatIn + volWSatOut) / (dynVisWSatIn + dynVisWSatOut);
			fm = a2 * pow(Reynolds, b2) + (a1 * pow(Reynolds, b1) - a2 * pow(Reynolds, b2)) / pow(1. + pow(Reynolds / t, c), d);

			dpdynSect = fm * 2. * rhoM * wM * wM * LengthSection / Dia;
			if (Base.showDPTubeDetail) {
				prot << "\n 128 Reynolds " << Reynolds << " fm " << fm << " dPdyn " << dpdynSect;
			}
		}
	}
	else {
		/**
		 * dynamic pressure drop of vertical tubes according Chisholm (VDI Heat Atlas, H2.2. equations 5 - 8).
		 */
		double MassVelW = MassVel * (1. - xMean);
		Reynolds = Dia * MassVelW / (dynVisWSatIn + dynVisWSatOut) * 2.;
		zeta = FrictFact(Reynolds, Base.Rough / Dia);
		//				pBar = (pMPaIn + pMPaOut) * 5.;
		double gamma = sqrt((1. / volWSatIn + 1. / volWSatOut) / (1. / volSSatIn + 1. / volSSatOut)) * pow((dynVisSSatIn + dynVisSSatOut) / (dynVisWSatIn + dynVisWSatOut), 0.1);
		r1 = 1. + (gamma * gamma - 1.) * (21. / gamma * pow(xMean * (1. - xMean), 0.9) + pow(xMean, 1.8));
		FrictCoeffTube = r1 * zeta * LengthSection / Dia;
		dpdynSect = FrictCoeffTube * (volWSatIn + volWSatOut) / 4. * MassVelW * MassVelW;
	}

	if (Base.showDPTubeDetail) {
		prot << "\n 146 FrictCoeffTube " << FrictCoeffTube;
		prot << "\n 1/volWin " << 1. / volWSatIn << " 1/volWout " << 1. / volWSatOut << " MassVel " << MassVel;
		prot << "\n head " << (volWSatIn + volWSatOut) / 4. * MassVel * MassVel << endl;
	}

	return dpdynSect;
} /* dpdyn_HeatAtlas */

