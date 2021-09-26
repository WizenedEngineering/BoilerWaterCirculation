/*****************************************************************//**
 * \file   KsiTee.cpp
  ******************************************************************** */
//#include "stdafx.h"
#undef MAINFUNCTION

#include "CommonHeader.h"
using namespace std;
double _branch::KsiTee(TFlow TCase, double Flowz, double Flowa, double Areaz, double Areaa,
	double Fillet, double angle, double dhyd) {

	/* Local variables */
	double FlowRatio, /* ratio of branching flow to main flow */
		AreaRatio,
		ksi, C, D, E, y;
	double d__1, d__2, fa, fm, alpha;

	/*!
	###Nomenclature :
	*  Flow : Flow kg/s\n
	*  Area : cross section m2\n
	*  ksi  : friction factor\n
		###Index:
		####a) Flow separation
	*      z arriving total flow\n
	*      a branching-off flow\n
	*      d flow in same direction as arriving flow\n
	*      Areaz = Aread\n
	*/
	//               |  ^  |
	//               |  |  |
	//               |     |
	//               |  a  |
	//     ----------       -----------------
	//       ->> z           -> d
	//     ----------------------------------
	/**
		####case TFlow::StraightSeparationInletStraight:
		####case TFlow::StraightSeparationOutletStraight:
		####case TFlow::StraightSeparationInletOff:
	* Source: F.Brandt, Dampferzeuger, FDBR Fachbuch Band 3\n
	![](@ref Tseparation.png)
			*\n
		*   if total flow is coming from branch-in tube\n
		*      z arriving total flow\n
		*      a, d flow in both inline tubes\n
		*      Areaa = Aread\n
	 */
	 //        ---------------------------
	 //            <- a           -> d
	 //        -----------     -----------
	 //                  |     |
	 //                  |  ^  |
	 //                  |  ^  |
	 //                  |  |  |
	 //                  |  z  |
	 //                  |     |
	  /**
	 ####case TFlow::OffSeparationInletStraight:
	 ####case TFlow::OffSeparationOutletOff:
	 * Source: IE Idel'chik, Handbook of hydraulic resistance\n
	 ![](@ref Tseparation1.png)

			 Idel'chik gives a formula in Diagram 7-29 for angle of 90 deg

			 ksi = 1 + k (velocity_a/velocity_z)^2

			 the pressure drop is ksi * velocity head of z
			 this allows the formula to be split up

			 dp = 1 * velocity head of z + k * velocity head of a (or d)

			 ksiOut of OffSeparationOutletOff is the first term = 1
			 ksiIn of OffSeparationInletStraight is k

			 Idel'chik gives k as 0.3 for welded tees and 1.5 for standard threaded malleable-iron tees

			 Comparing to sharp-edged smooth tube deflection k would be 1.32. As the velocity head for z is already taken 0.32 remains.

			 To handle other angles than 90 deg a sharp-edged deflection with the corresponding angle is assumed
	 ####b) Flow union
		 *      z leaving total flow\n
		 *      a flow from branching-in tube\n
		 *      d flow in same direction as leaving flow\n
		 *      Aread = Areaz\n
	 ####case TFlow::StraightUnionOutletOff:
	 ####case TFlow::StraightUnionOutletStraight:
	 ####case TFlow::StraightUnionInletStraight:
	 * Source: F.Brandt, Dampferzeuger, FDBR Fachbuch Band 3
	  */
	  //        ---------------------------\n
	  //            -> d           ->> z\n
	  //        -----------     -----------\n
	  //                  |     |\n
	  //                  |  ^  |\n
	  //                  |  |  |\n
	  //                  |  a  |\n
	  //                  |     |\n
	  /**
		 ![](TUnion.png)
		 * \n
		  * if total flow is leaving to branch-off tube\n
		  *      z leaving total flow\n
		  *      a, d flow in both inline tubes\n
		  *      Areaa = Aread\n
			  ####case TFlow::OffUnionOutletStraight:
			  ####case TFlow::OffUnionInletOff:
	  * Source: IE Idel'chik, Handbook of hydraulic resistance\n
		![](TUnion1.png)
		*/
		//                 |  ^  |
		//                 |  ^  |
		//                 |  |  |
		//                 |  z  |
		//       ----------       -----------------
		//            a ->         <- d
		//       ----------------------------------

		/**
		particular:

				Idel'chik gives a formula in Diagram 7-29 for angle of 90 deg

				the diagram only show cases with areaa<= areaz

				his formula is used in the case of 90deg and equal cross section

				the pressure drop is ksi * velocity head of z
				this allows the formula to be split up

				dp = 1 * velocity head of z + factor * velocity head of a(or d)

				ksiIn of OffUnionInletOff is the first term = 1
				ksiOut of OffUnionOutletStraight is factor

				To handle other angles than 90 deg or different cross section a sharp - edged deflection with the corresponding
				angle is assumed for OffUnionOutletStraight

		*/
		// ----------------------------------------------------------
	//            prot << " FlowA " << Flowa << " Flowz " << Flowz << endl;
	//            prot << " AreaA " << Areaa << " Areaz " << Areaz << endl;

	if (Flowa == 0. || Flowz == 0. || Areaz == 0. || Areaa == 0.) {
		cout << "\n Tee : Flow or Area = 0 ";
		prot << "\n Tee : Flow or Area = 0 ";

		return 0.;
	}
	switch (TCase) {
	case TFlow::StraightSeparationInletOff:
		/*     ---------------------------------------------------------- */
		/*     Flow separation inlet of branch-off tube                   */
		/*     ---------------------------------------------------------- */
		AreaRatio = Areaz / Areaa;
		FlowRatio = Flowa / Flowz;
		y = Fillet / dhyd;
		C = 0.87 * (AreaRatio - 0.82);
		D = 0.0685 * (1.36 + AreaRatio);
		E = 0.1 + 1.32 * exp(-0.4 * (AreaRatio - 1.));
		d__1 = (1. - FlowRatio);
		d__2 = AreaRatio * FlowRatio;
		ksi = 2. * (0.98 * d__1 * d__1 - 1. + d__2 * (FlowRatio * (C + D * exp(-40. * y)) + E - 0.6 * cos(angle))) - (d__2 * d__2 - 1.);
		//      prot << "case 1 ksiz " << ksi << " AreaRatio " << AreaRatio << " FlowRatio " << FlowRatio << endl;
				  // correction to inlet of a
		ksi = ksi / d__2 / d__2;
		//         prot<<" Ksia "<<ksi<<endl;
		break;
	case TFlow::StraightSeparationInletStraight:
		/*     ---------------------------------------------------------- */
		/*     Flow separation inlet of leaving straight tube             */
		/*     ---------------------------------------------------------- */
		AreaRatio = Areaz / Areaa;
		FlowRatio = Flowa / Flowz;
		C = 0.04 * (1. - 1.25 * exp(-0.68 * AreaRatio));
		d__1 = (1. - FlowRatio);
		//         d__2 = AreaRatio * FlowRatio;
		ksi = 2. * (1.02 * d__1 * d__1 - 1. + 0.8 * FlowRatio * (1. - 0.0139 * FlowRatio + C * cos(angle))) - (d__1 * d__1 - 1.);
		//         prot << "case 2 ksiz " << ksi << " AreaRatio " << AreaRatio << " FlowRatio " << FlowRatio << endl;
					// correction to inlet of d
		ksi = ksi / d__1 / d__1;
		//                 prot<<" Ksia "<<ksi<<endl;

		break;
	case TFlow::StraightSeparationOutletStraight:
		/*     ---------------------------------------------------------- */
		/*     Flow separation outlet of arriving straight                */
		/*     ---------------------------------------------------------- */
		ksi = 0.;
		//                  prot<<" Ksia "<<ksi<<endl;
		break;
	case TFlow::StraightUnionOutletOff:
		/*     ---------------------------------------------------------- */
		/*     Flow union outlet of branch-in tube                        */
		/*     ---------------------------------------------------------- */
		y = Fillet / dhyd;
		AreaRatio = Areaz / Areaa;
		FlowRatio = Flowa / Flowz;
		d__1 = 1. - FlowRatio;
		if (AreaRatio < 0.04) {
			C = 0.;
		}
		else {
			C = 0.25 * pow((AreaRatio - 0.04), 0.26);
		}
		E = 1.29 * pow(AreaRatio, -0.17) - 0.94;
		d__2 = AreaRatio * FlowRatio;
		//                 prot << " d__1 " << d__1 << " C " << C << " E " << E << " d__2 " << d__2 << " y " << y << " angle " << angle << endl;
		ksi = 2. * (1. - 0.97 * d__1 * d__1 - d__2 * (FlowRatio * (C - 0.1 * exp(-40. * y) + 0.4 * cos(angle)) + E)) - (1. - d__2 * d__2);
		//                  prot << "case 3 ksiz " << ksi << " AreaRatio " << AreaRatio << " FlowRatio " << FlowRatio << endl;
		ksi = ksi / d__2 / d__2;
		//          prot<<" Ksia "<<ksi<<endl;
		break;
	case TFlow::StraightUnionOutletStraight:
		/*     ---------------------------------------------------------- */
		/*     Flow union outlet of arriving straight tube                */
		/*     ---------------------------------------------------------- */
		AreaRatio = Areaz / Areaa;
		FlowRatio = Flowa / Flowz;
		d__1 = 1. - FlowRatio;
		C = 1.05 - 0.475 * exp(-2.26 * (AreaRatio - 1.));
		E = 0.468 * pow(AreaRatio, -2.17);
		d__2 = AreaRatio * FlowRatio;
		ksi = 2. * (1. - 0.981 * d__1 * d__1 - d__2 * (FlowRatio * (C * cos(angle) - 0.11) + E)) - (1. - d__1 * d__1);
		//                  prot << "case 4 ksiz " << ksi << " AreaRatio " << AreaRatio << " FlowRatio " << FlowRatio << endl;
		//                  prot<<"Flow a "<<Flowa<<" flow z "<<Flowz<< " Ksia "<<ksi<<endl;
		break;
	case TFlow::OffSeparationOutletOff:
		// Flow separation from branching tube
		ksi = 1.;
		break;
	case TFlow::OffSeparationInletStraight:
		// Flow separation from branching tube inlet of straight tube
	 //handled as sharp bend
//assuming a smooth tube as well as reynolds > 4e4
		d__2 = sin(angle / 2.);
		d__2 *= d__2;
		fm = d__2 * (.95 + d__2 * 2.05);
		fa = 1.2;
		alpha = angle / M_PI * 180.; // angle in deg
		if (alpha <= 87.1) { //fitting of data from diagram 6-7b and table (Idel'chik)
			fa = (0.00019 * alpha - 0.04) * alpha + 3.24265;
		}
		ksi = 1.1 * fm * fa - 1.; // 1 times velocity head is already handled in incoming branch as outlet resistance
		ksi = fmax(0., ksi);
		break;
	case TFlow::StraightUnionInletStraight:
		// Flow union from branching tube, inlet of straight tube
		ksi = 0.;
		break;
	case TFlow::OffUnionOutletStraight:
		// Flow union from straight tubes into branching tube, outlet of straight tube
		if (fabs(angle - M_PI / 2.) < 1e-3 && fabs(Areaa / Areaz - 1.) < 1e-6) {
			ksi = 3. + Flowz / Flowa * (Flowz / Flowa - 3.);
		}
		else {
			//assuming sharp-edged deflection with smooth tube as well as reynolds > 4e4
			d__2 = sin(angle / 2.);
			d__2 *= d__2;
			fm = d__2 * (.95 + d__2 * 2.05);
			fa = 1.2;
			alpha = angle / M_PI * 180.; // angle in deg
			if (alpha <= 87.1) { //fitting of data from diagram 6-7b and table (Idel'chik)
				fa = (0.00019 * alpha - 0.04) * alpha + 3.24265;
			}
			ksi = 1.1 * fm * fa;
		}
		break;
	case TFlow::OffUnionInletOff:
		// Flow union from straight tubes into branching tube, inlet of off tube
		 // handled as sharp edged bend
		ksi = 1.;
		break;
	default:
		cout << "T-Piece: case not available";
		ksi = 0.;
		break;
	}
	return ksi;
} /* tstueck_ */

