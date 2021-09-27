//#include "stdafx.h"
/*     ----------------------------------------------------- */
/*     function to determine resistance factor of bends      */
/*     ----------------------------------------------------- */

#include "CommonHeader.h"

double _tube::ksiBend(double Reynolds, double AngleDeg) {
	/* Local variables */
	int j;
	double  us[7], ksiBend = 0., d__2;
	double eta = 0.;
	double fh = 2.;
	double fm = 1.;
	double fa = 1.;
	double fr = 1.;
	double RelRough = Base.Rough / Dia;
	double RelRadius = RadiusBend / Dia;
	/*!
	 * Different calculation is done whether it is a smooth or sharp edged bend\n
	 * In sharp edged bends the bending radius (at tube centerline) is small\n
	 * The criterion is Radius/Diameter less than 0.1. A sharp inner edge is possible below Radius/Diameter less than 0.5.\n
	 *
	 */
	 //   double zeta = FrictFact(Reynolds, Dia);
	 //   cout << "\nreyn " << Reynolds << " zeta " << zeta << " radius " << RadiusBend << " AngleDeg " << AngleDeg;
	if (RelRadius > 0.1) { //smooth bend 
		if (AngleDeg <= 70.) {
			fa = sin(AngleDeg * M_PI / 180.) * .9;
		}
		else if (AngleDeg >= 100.) {
			fa = 0.7 + 0.35 * AngleDeg / 90.;
		}
		if (RelRadius > 1.) {
			fm = fa * 0.21 / sqrt(RelRadius);
		}
		else {
			fm = fa * 0.21 / pow(RelRadius, 2.5);
		}
		if (RelRadius < 0.55) {
			if (Reynolds < 4e4) {
				fr = FrictFact(Reynolds, RelRough) / FrictFact(4e4, RelRough);
			}
		}
		else {
			if (Reynolds < 2e5) {
				fr = FrictFact(Reynolds, RelRough) / FrictFact(2e5, RelRough);
			}
		}

		//! for smooth bends:\n
		//! at Reynolds numbers close to 4e4 factor fh can jump between 1 and 2 (in worst case) depending on relative roughness\n
		//! numeric instability at flows that lead to Reynolds close to 4e4 can happen\n
		//! minimum flow difference but one time above=> high pressure drop, one time below => significantly lower pressure drop\n  
		//! Reynolds is set to minimum 4e4 for factor fh\n 
		//! this might lead to higher ksiBend at low Reynolds numbers\n 

		Reynolds = fmax(Reynolds, 4e4);
		if (RelRadius < 0.55) {
			fh = fmin(1.5, 1. + 500. * RelRough);
		}
		else if (RelRadius < 1.5) {
			if (RelRough < 1e-3) {
				if (Reynolds > 2e5) {
					fh = 1. + 1e3 * RelRough;
				}
				else {
					fh = FrictFact(Reynolds, RelRough) / FrictFact(Reynolds, 0.);
				}
			}
		}
		else { // relRadius > 1.5
			if (RelRough < 1e-3) {
				fh = 1. + RelRough * RelRough * 1e6;
			}
		}
		ksiBend = fh * fr * fm;

		/*     ------------------------------------------ */
		/*     correction for bends close to each other in
		 *     "S" or "U" arrangement                     */
		 /*     ------------------------------------------ */
		if (UorS != USArrangement::No) {
			double dist = Length - 2. * RadiusBend;
			if (dist < 1e-3) {
				if (UorS == USArrangement::S) {
					eta = 1.1964706; // is  ‬ 1.0982353 *2. -1.;
				}
				else {
					eta = 0.21844428; // is  .60922214* 2.-1.;
				}
			}
			else {
				if (UorS == USArrangement::S) {
					us[0] = 1.0982353;
					us[1] = -.20151263;
					us[2] = .05792232;
					us[3] = -.0080823559;
					us[4] = 6.0783983e-4;
					us[5] = -2.3790387e-5;
					us[6] = 3.8269694e-6;
				}
				else {
					us[0] = .60922214;
					us[1] = .060370803;
					us[2] = -.0072293581;
					us[3] = 7.8684409e-4;
					us[4] = -8.0123244e-5;
					us[5] = 4.8920378e-6;
					us[6] = -1.1609907e-7;
				}
				for (j = 0; j < 7; ++j) {
					eta += us[j] * pow(dist / Dia, j);
				}

				//! For smooth bends the arrangement of consecutive bends in U- or S-arrangement is taken into account.\n 
				//! a correction factor "eta" is calculated\n 
				//! for both bends the total ksiBend = 2 * eta * ksiBend(single)\n
				//! eta is only applied to the first bend,\n
				//! assuming ~ equal ksiBend for both bends the first bend takes full advantage\n 
				//! ksiBend = ksiBend *(2*eta - 1) and the second uses eta = 1\n
				if (eta > 0.5) {
					ksiBend *= (2. * eta - 1.);
					if (Base.showDPTubeDetail) {
						prot << "\neta < 0  Length " << Length << "  Dia " << Dia << " UorS ";
					}
				}
			}
		}
	}
	else { // sharp edged bend
	/* Computing 4th power */
		d__2 = sin(AngleDeg * M_PI / 360.);
		d__2 *= d__2;
		fm = d__2 * (.95 + d__2 * 2.05);
		fa = 1.2;
		if (AngleDeg <= 87.1) { //fitting of data from diagram 6-7b and table (Idel'chik)
			fa = (0.00019 * AngleDeg - 0.04) * AngleDeg + 3.24265;
		}
		fh = 1.;
		if (Reynolds >= 4e4) {
			fh = fmin(1.5, 1. + 500. * RelRough);
			if (RelRough <= 1e-6) {
				fr = 1.1; //for smooth tube fh still is ~1 but fr should be 1.1
			}
		}
		else {
			fr = FrictFact(Reynolds, RelRough) / FrictFact(4e4, RelRough);
		}
		//      cout << "\n fh " << fh << " fr " << fr << " fa " << fa << " fm" << fm;
		ksiBend = fh * fr * fa * fm;
		//      cout << " ksibend " << ksiBend;
	}
	//   prot<<" AngleDeg "<<AngleDeg<< " zeta "<<zeta<<endl;
	//   prot << " 1phase ksiBend " << ksiBend << endl;

	if (Base.showDPTubeDetail) {
		if (Base.showDPTube) {
			prot << "\n  Bend";
			prot << " \n  reynolds   AngleDeg      eta        r     dhyd    fh       fr    "
				"   fm       fa    ksiBend ";
			prot << "\n" << setw(8) << Reynolds << " " << setw(5) << AngleDeg << " " << setw(8)
				<< eta << " " << setw(8) << RadiusBend << " " << setw(8) << Dia << " "
				<< setw(8) << fh << " " << setw(8) << fr << " " << setw(8) << fm << " "
				<< setw(8) << fa << " " << " " << setw(8) << ksiBend;
		}
	}
	return ksiBend;
} /* Bend */
