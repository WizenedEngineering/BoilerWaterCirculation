//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"
// the number of sections in the tube can vary with flow (different steam quality or void fraction)
// it is only calculated once for the smallest flow giving the higher steam quality / void fraction

void _tube::DetermineNoSections(double g) {
	int ns;
	/**
	 * criteria
	 *
	 * 1) tubes longer than 2*diameter or bends with r/D > 5 are split into sections
	 *
	 */
	NoSections = 1;
	if ((RadiusBend < 1e-3 && Length / Dia > 2.) ||
		(RadiusBend > 1e-3 && RadiusBend / Dia > 5.) 
		) {
		//      if (Number == 363) {
		//         prot << "\n determine Length " << Length << " dia " << Dia;
		//         prot << "\n q " << q;
		//      }
		if (q < 1e-3) {
			ns = static_cast<int> (Length); // number of tube section calculated in this step
			if (ns > 1)NoSections = ns;

		}
		else {
			ns = static_cast<int> (2. * Length); // number of tube section calculated in this step
			if (ns > 1)NoSections = ns;
			double deltaEnth = q / g; //enthalpy difference in tube
			if (deltaEnth / 10. >= 2.) {
				ns = static_cast<int> (deltaEnth / 10. + 0.5); ///3.) heated tubes: in each subsection enthalpy difference should not be more than 10 kJ/kg
			}
			if (ns > NoSections) NoSections = ns;
			//         if (Number == 363) {
			//            prot << "\n from deltaH " << NoSections << " deltaH " << deltaEnth;
			//         }
						/// 4.) heated tubes: additional limitation of change of void fraction (for simplicity homogeneous at inlet condition), should not be higher than 5% 
			double pMPaIn = pPaIn * 1e-6;
			double tSat = H2O::satTemp(pMPaIn);
			double rhoSSat = 1. / H2O::specVol(tSat, pMPaIn, STEAM);
			double rhoWSat = 1. / H2O::specVol(tSat, pMPaIn, WATER);
			double enthWSat = H2O::enth(tSat, pMPaIn, WATER);
			double enthSSat = H2O::enth(tSat, pMPaIn, STEAM);
			double xTbIn = (EnthIn - enthWSat) / (enthSSat - enthWSat);
			double VoidFractIn;
			if (xTbIn < 0.) {
				VoidFractIn = 0.;
			}
			else if (xTbIn > 1.) {
				VoidFractIn = 1.;
			}
			else {
				VoidFractIn = rhoWSat * xTbIn / (rhoWSat * xTbIn + rhoSSat * (1. - xTbIn));
			}
			double xTbOut = (EnthIn + deltaEnth - enthWSat) / (enthSSat - enthWSat);
			double VoidFractOut;
			if (xTbOut < 0.) {
				VoidFractOut = 0.;
			}
			else if (xTbOut > 1.) {
				VoidFractOut = 1.;
			}
			else {
				VoidFractOut = rhoWSat * xTbOut / (rhoWSat * xTbOut + rhoSSat * (1. - xTbOut));
			}

			ns = static_cast<int> (fabs((VoidFractOut - VoidFractIn) / 0.05 + 0.5));
			//         if (Number == 363) {
			//            prot << "\n xTbIn " << xTbIn << " voidIn " << VoidFractIn << " xTbOut " << xTbOut << " voidOut " << VoidFractOut;
			//            prot << "\n deltaEnth " << deltaEnth << " (VoidFractOut - VoidFractIn) " << (VoidFractOut - VoidFractIn) << " ns " << ns << " sections " << NoSections << endl;
			//         }
			if (ns > NoSections) NoSections = ns;
			if (Base.showDPTubeDetail) {
				prot << "\n deltah" << deltaEnth << " NoSection " << NoSections << endl;
			}
		}
	}
	///5.) after determine again additional check: the length of a section should not be smaller than diameter (except total length is smaller than diameter)
	if (Length / NoSections < Dia) {
		if (Base.showDPTubeDetail) {
			prot << "\n l<dia  Length/dia " << Length / Dia << endl;
		}
		if (Length / Dia >= 2.) {
			NoSections = static_cast<int> (Length / Dia);
			if (Base.showDPTubeDetail) {
				prot << " l<dia NoSection" << NoSections << endl;
			}
		}
		else {
			NoSections = 1;
		}
	}
	LengthSection = Length / NoSections;
	HeightSection = Height / NoSections;
	HeatSection = q / (NoSections * NoParallel);

	if (Base.showDPTube) {
		//   if(Number == 363){
		prot << "\n NoSections " << NoSections << " h " << HeightSection << " l " << LengthSection << " dia " << Dia << " q " << q << " gcalc " << g;
		prot << "\n pIn " << pPaIn << " enthIn " << EnthIn << " enthOut " << EnthOut << " enthdrum " << Drum.enthW << endl;
	}
}
