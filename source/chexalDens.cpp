//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

//common variables used in Chexal functions

static double jFstar,  ///<superficial fluid velocity at CCFL(ft/sec)
					jGstar, ///<superficial gas velocity at CCFL(ft/sec)
					alstar; ///<void fraction at CCFL 

static double rhoF,///< density of fluid (water) (lbm/ft^3)
					rhoG,///<density of gas (vapor,steam) (lbm/ft^3)
					sig,///<surface tension (lbf/ft) 
					rr,///<density ratio (-)
					cp,///<pressure factor (-)
					vj0xx; ///<non flow and alpha part of VGJ 
//Properties in British Units

static double dh, ///<hydraulic diameter (ft) 
					ReFxx,///<partial liquid Reynolds number (ReF/jF) 
					ReGxx,///<partial vapor Reynolds number (ReG/jG) 
					fra;///<tube orientation factor (-), tube orientation between 0 deg(vertical) and 90 deg (horizontal)  

/**
 * Comment in the original FORTRAN files:\n
 *     The  PCLV  subroutine  has  been  classified  by  EPRI  as  a\n
 *     Research Program. As such, it has not been developed and tested\n
 *     to the extent that a production program would be, and unforeseen\n
 *     results may occur when running the program.  EPRI does not make\n
 *     any warranty or representation whatsoever, expressed or implied\n
 *     including any warranty of  merchantability or fitness for any\n
 *     purpose with respect to the program.  Nor does EPRI assume any\n
 *     liability whatsoever with respect to any use of the program or\n
 *     any portion thereof or with  respect to any damages which may\n
 *     result from such use.

 *     This subroutine is based on the 1996 Revision of the EPRI
 *     Chexal-Lellouche Void Fraction model as documented in
 *     EPRI Report TR-106326, Void Fraction Technology for Design
 *     and Analysis, December 1996. */

 /**
  * @brief calculates the void fraction
  *
  * sets up jF and jG and calls chexal_glenns\n
  * All internal calculations are performed in Imperial units. Conversion is done in calling function
  * @param [in] x steam quality
  * @param [in] jF Superficial Liquid Velocity (ft/sec)
  * @param [in] jG Superficial Vapor Velocity (ft/sec)
  * @param [out] v1 Void Fraction
  * @param [out] c01 Co, Concentration Parameter
  * @param [out] vgj1 Vgj, Drift Velocity (same units as input)
  * @param [out] v2  2nd Void Fraction For Counter-Current Flows
  * @param [out] c02 Co, Concentration Parameter for 2nd void fraction
  * @param [out] vgj2 Vgj, Drift Velocity for 2nd void fraction (m/s or ft/sec)
  * @return int Error flag, set to nonzero value if trouble is encountered\n
  *         1 Cannot find 1st void fraction\n
  *         2 Cannot find 2nd void fraction, CC Flow only\n
  *         3 Cannot find either void fraction\n
  *         4 Bad CounterCurrent Flow Pair jF > jFSTAR\n
  *         5 CCFL Solution in FRESCO did not converge in allotted number of iterations
 */
int chexal_pclv(double x, double jF, double jG, double& v1, double& c01,
	double& vgj1, double& v2, double& c02, double& vgj2);

/**
 *
 * @brief calculates void fraction
 *
 * in counter flow cases (steam flowing opposite direction to water) we can have 2 solutions \n
 * both solutions are in in the arrays (gvoid, gc0 and gvgj)
 * @param [in] x steam quality
 * @param [in] jF liquid superficial velocity
 * @param [in] jG vapor superficial velocity
 * @param [out] gvoid array of void fraction
 * @param [out] gc0 array of c0
 * @param [out] gvgj array of drift velocity
 * @return int error code
 */
int chexal_glenns(double x, double jF, double jG, double* gvoid, double* gc0,
	double* gvgj);

/**
 * @brief chexal_c3func computes the Chexal-Lellouche C3 coefficient
 *
 * Chexal_c3func is called from ::chexal_glenns and ::chexal_fresco
 *
 * @param [in] ReF superficial liquid Reynolds number
 * @param [in] ReG superficial vapor Reynolds number
 * @param [in] jF liquid superficial velocity
 * @param [in] jF1 CCFL liquid superficial velocity (jFstar) or liquid superficial velocity (jF) depends on calling function
 * @return double C3
 */
double chexal_c3func(double ReF, double ReG, double jF, double jF1);

/**
 * @brief chexal_c0func_init calculates several constants for C0FUNC in the Chexal-Lellouche drift flux correlation.
 *
 * The parameters computed here are orientation dependent.
 * This function is called from ::chexal_glenns and ::chexal_fresco
 *
 * @param [in] ReF superficial liquid Reynolds number
 * @param [in] ReG superficial vapor Reynolds number
 * @param [out] c1x exponent used for (1-alpha) for C1 constant
 * @param [out] b1 B1 constant
 * @param [out] k0 parameter used in defining fluid parameter L
 * @param [out] R r parameter from the correlation
 * @return int error code
 */
int chexal_c0func_init(double ReF, double ReG,
	double& c1x, double& b1, double& k0, double& R);

/**
 * @brief  C0FUNC calculates several constants in the Chexal-Lellouche drift flux correlation.
 *
 * The parameters computed here are orientation dependent.
 *
 * This function is called from ::chexal_glenns and ::chexal_fresco
 *
 * @param [in] jSum sum of superficial velocity
 * @param [in] jG vapor superficial velocity
 * @param [in] vgj0 VGJ/C1 where C1 = ALM1**C1X
 * @param [in] al vapor volume fraction
 * @param [in] alm1 liquid volume fraction
 * @param [in] k0 parameter used in defining fluid parameter L
 * @param [in] R  r parameter from the correlation

 * @return double C0  distribution parameter
 */
double chexal_c0func(double jSum, double jG, double vgj0, double al,
	double alm1,  double k0, double R);

/**
 * @brief chexal-fresco is used to compute a CCFL (counter current flow limit) value of void fraction
 *     and jF (jFstar) on the CCFL line given a jG as input.
 *
 * Fresco is called from ::chexal_glenns
 * @param [in] x steam quality
 * @param [in] jF value of jF (MUST BE NEGATIVE)
 * @param [in] jG value of jG (MUST BE POSITIVE)
 * @param [in] CalcjFStar true if jFstar is computed, XjFIN IS THE INITIAL GUESS FOR THE ANSWER (jFSTAR)\n
 *                false if jGstar is computed, XjGIN IS THE INITIAL GUESS FOR THE ANSWER (jGSTAR)
 * @param [out] answer jFstar or jGsta
 * @param [out] Void void fraction on the CCFL line
 * @return int
 */
int chexal_fresco(double x, double jF, double jG, bool CalcjFStar,
	  double* answer, double* Void);

/**
 * @brief iteration of void fraction
 *
 * chexal_iterat solves the Chexal-Lellouche correlation for jF or jG and void fraction\n
 * It uses the regula-falsi method (Anderson-Bjoerk) for co-current flow\n
 * For counter-current flow the stepped method from original code is kept to find the 2 solutions.
 *
 * chexal_iterat is called from subroutine ::chexal_fresco (iccfl = 1) and ::chexal_glenns (iccfl = 0).
 * @param [in] x steam quality
 * @param [in] jF superficial liquid velocity
 * @param [in] jG superficial vapor  velocity
 * @param [out] gvoid void fraction
 * @param [in] iccfl parameter what should be calculated
 * @param [in,out] c3 C3 constant
 * @param [out] gc0 distribution parameter
 * @param [in] vgj0 VGJ/C1 where C1 = ALM1**C1X
 * @param [in] c1x exponent used for (1-alpha) for C1 constant
 * @param [out] gvgj drift velocity
 * @param [in] k0 parameter used in defining fluid parameter L
 * @param [in] R  r parameter from the correlation
 * @return int
 */
int chexal_iterat(double x, double jF, double jG,
	double* gvoid, bool iccfl, double& c3, double* gc0,
	double vgj0, double c1x, double* gvgj,
	double k0, double R);

/**
 * @brief iteration of void fraction with given boundaries
 *
 * It uses the regula-falsi method (Anderson-Bjoerk) \n
 * The final result has to be between the 2 boundaries\n
 * chexal_regula_falsi is called from subroutine chexal_iterat.
 * @param [in] a void fraction at lower boundary
 * @param [in] fa function value at lower boundary
 * @param [in] b void fraction at upper boundary
 * @param [in] fb function value at upper boundary
 * @param [in] jSum sum of superficial velocities
 * @param [in] jG superficial vapor velocity
 * @param [in] vgj0 VGJ/C1 where C1 = ALM1**C1X
 * @param [in] c1x exponent used for (1-alpha) for C1 constant
 * @param [in] k0 parameter used in defining fluid parameter L
 * @param [in] R  parameter from the correlation
 * @param [out] vgj drift velocity
 * @param [out] c0 drift flux distribution parameter
 * @return double
 */

double chexal_regula_falsi(double a, double fa, double b, double fb, double jSum, double jG, double vgj0, double c1x,
	double k0, double R,  double& vgj, double& c0);

/**
 * \brief  Chexal_poly calculates the difference between the jG calculated from the void fraction and the actual jg
 *
 * this difference is used by the chexal_iterat function
 *
 * chexal_poly is called by function ::chexal_iterat
 *
 * @param [in] al void fraction
 * @param [in] jSum sum of superficial velocities
 * @param [in] jG superficial vapor velocity
 * @param [in] vgj0 VGJ/C1 where C1 = ALM1**C1X
 * @param [in] c1x exponent used for (1-alpha) for C1 constant
 * @param [in] k0 parameter used in defining fluid parameter L
 * @param [in] R r parameter from the correlation
 * @param [out] vgj drift velocity
 * @param [out] c0 drift flux distribution parameter
 * @return double
 */
double chexal_poly(double al, double jSum, double jG,
	double vgj0, double c1x, double k0, double R,
	 double& vgj, double& c0);

/**
 * \brief MARCH2 Checks convergence of a solution on the CCFL line
 *
 * and increments jG or jF for the next try at the solution of the Chexal-Lellouche drift flux.
 *
 * MARCH2 is called from function ::chexal_fresco
 *
 * \param [in] CalcjFStar = true if jFstar is computed\n
 *             = false if jGstar is computed
 * \param [in] key = 1, converged solution found\n
 *             = 2, converged solution not found
 * \param [in,out] jF superficial liquid velocity
 * \param [in,out] jG superficial vapor velocity
 * \param [in,out] jFold  old superficial liquid velocity
 * \param [in,out] jGold old superficial vapor velocity
 * \param [out] delta change in jF or jG for the next iteration
 * \param [out] ihist = 1, if solution found for last try\n
 *              = 2, if solution NOT found for the last try
 * \param [out] iboth = 0, last two values of jF or jG produced were on one side of the CCFL line\n
 *              = 2, values of jF or jG have been found on both sides of the CCFL line.
 * \return int
 */
int chexal_march2(bool CalcjFStar, int key, double& jF, double& jG,
	double& jFold, double& jGold, double& delta,
	int& ihist, int& iboth);

#define ISIGN(x) ((x) >= 0. ? 1 : -1)

// Properties are implicit at saturation (no ..Sat.. in variable name)

double _tube::Density_chexal(double x, double pMPa, double volW,
	double volS, double dynVisW, double dynVisS, double SurfTens, double& v1) {
	/* System generated locals */
	double dens; //mixture density [kg/m3]
	/* Local variables */
	static const double gravty = 32.17; // gravity constant, British Units
	//   double v1, // void fraction
	double v2, // second void fraction for counter-current flow
		jF, // superficial fluid (water) velocity [ft/s]
		jG; // superficial gas (steam) velocity [ft/s]
	int error = 0;
	double cz1, cz2;
	double vgj1, vgj2;

	double d__1;
	double c2, c4, c5, c7, c8, arg;
	jFstar = 0.;
	jGstar = 0.;
	alstar = 0.;

	/* if superficial velocity too low return with single phase */
	jF = (1. - x) * MassVel * volW / .3048; /* in ft/s */
	if (jF < 1e-7) { //no fluid, only steam
		v1 = 1.;
		cout << "\n jF < 1e-7 shortcut steam";
		return 1. / volS;
	}
	jG = x * MassVel * volS / .3048; /* in ft/s */
	if (jG < 1e-7) { // no steam, only fluid
		v1 = 0.;
		cout << "\njG < 1e-7 shortcut water";
		return 1. / volW;
	}

	/*  PROGRAM TO REDUCE ADIABATIC VOID FRACTION DATA USING THE CHEXAL- */
	/*  LELLOUCHE VOID FRACTION CORRELATION - WILL ALSO CALCULATE CCFL. */

	/*  GET PROPERTIES */
	// conversion to British Units */
	//    p = pMPa * 145.04; /* P   = pressure */
	rhoF = 1. / volW * 0.0624228; /* rhoF = liquid density */
	rhoG = 1. / volS * 0.0624228; /* rhoG = vapor density */
	sig = SurfTens * 0.22481 / 3.2808; /* SIG  = surface tension */
	cp = 22.129 * 22.129 * 4. / (pMPa * (22.129 - pMPa)); /* CP   = 4*Pcrit**2/P*(Pcrit-P) */
	rr = rhoF / rhoG; /* RR   = RHOF/RHOG */
	if (rr >= 18.) {
		c2 = 1.;
		//c5 = sqrt(150. / rr); 
		//if (c5 < 1.) {
		//	arg = c5 / (1. - c5);
		//	if (arg < 85.) {
		//		c2 = 1. - exp(-arg);
		//	}
		//}
		if (rr > 153.550173) { //arg = 85; c5 = 85/86; -> rr = 153.550173 
			c5 = sqrt(150. / rr); //if rr >150 -> c5 < 1
			arg = c5 / (1. - c5);
			c2 = 1. - exp(-arg);
		}
	}
	else {
		c2 = pow(log(rr), 0.7) * .4757;
	}
	dh = Dia / .3048; /* in ft *//*     DH    = hydraulic diameter */
	if (dh > 0.3) {
		c7 = pow(0.3 / dh, 0.6);
		//	if (c7 < 1.) {
		c8 = c7 / (1. - c7);
		c4 = 1. / (1. - exp(-c8));
	}
	else {
		c4 = 1.;
	}
	d__1 = (rr - 1.) * sig * gravty * gravty / (rr * rhoF);
	vj0xx = c2 * c4 * 1.41 * pow(d__1, 0.25); /* VJ0XX = non flow and alpha part of VGJ */

	//      prot<<"\n chexal dens x "<<x<<" MassVel "<<MassVel<<" jF "<<jF<<" jG "<<jG<<endl;
	if (HeightSection < 0.) {
		jF = -jF;
		jG = -jG;
	}
	double ort;
	if (fabs(HeightSection) < 1e-6) { // horizontal
		ort = 90.;
		fra = 0.;
	}
	else if (fabs(fabs(HeightSection / LengthSection) - 1.) < 1e-6) { // vertical
		ort = 0.;
		fra = 1.;
	}
	else {
		ort = 90. - asin(fabs(HeightSection / LengthSection)) *180. / 3.14159265;
		fra = pow(1. - ort / 90., 0.2);
	}
	/*     Subroutine BASIC calculates Chexal Lellouche drift flux */
	/*     correlation variables that are not flow or void dependent */
	/*     This subroutine is incorporated here */

	ReFxx = rhoF * dh * 3600. / (dynVisW * 2419.1); /*     REFXX = partial liquid Reynolds number (REf/jF) */
	ReGxx = rhoG * dh * 3600. / (dynVisS * 2419.1); /*     REGXX = partial vapor Reynolds number (REg/jG) */


 //  prot<<"\n psi"<<pMPa * 145.04<<" dh "<< dh<<" jF "<<jF<<" jG "<<jG<< " rhoF "<< rhoF << " rhoG "<< rhoG << " sigma "<< sig << " ReFxx "<< ReFxx<<" ReGxx "<<ReGxx;
	if (Base.showDPTubeDetail) {
	cout << "\n           -------Steam-Water Results-------\n" <<
		" System Pressure - " << pMPa << " MPa";
	cout << "\n  Hydraulic Diameter = " << Dia << " m" <<
		"\n Flow Direction from Vertical = " << ort << " degrees" <<
		"\n x = " << x <<
		"\n l =  " << LengthSection << " h " << HeightSection;
		}

		/*   CALL VOID FRACTION SUBROUTINE */
	error = chexal_pclv(x, jF, jG, v1, cz1, vgj1, v2, cz2, vgj2);
	if (error > 0) {
		cout << "\n ERROR IN VOID SUBROUTINE, IER = " << error;
		if (error == 4) {
			cout << " THE FOLLOWING RESULTS ARE INVALID" << endl;
		}
		else {
			cout << " THE FOLLOWING RESULTS MAY BE INVALID" << endl;
		}
		if (Base.showDPTube) {
			prot << "\n ERROR IN VOID SUBROUTINE, IER = " << error << endl;
			if (error == 4) {
				prot << "\n THE FOLLOWING RESULTS ARE INVALID";
			}
			else {
				prot << "\n THE FOLLOWING RESULTS MAY BE INVALID";
			}
		}
	}
	if (Base.showDPTubeDetail) {
	prot << "\n jF = " << jF * .3048 << " jG = " << jG * .3048;
	prot << "\nVoid Fraction = " << v1 << " Co = " << cz1 << " Vgj = " << vgj1 * 0.3048 << " m/s";
	if (jF * jG < 0.) {
		prot << "\nCounter Current Flow Second Root" <<
			"\nVoid Fraction = " << v2 << "  Co = " << cz2 << "  Vgj = " << vgj2 * 0.3048 << " m/s";
		}
	}
	dens = (1. - v1) / volW + v1 / volS; // in metric units
//   prot << "\n dia " << Dia << " height " << HeightSection << " length " << LengthSection << " ort " << ort <<
//      " pMpa" << pMPa << " x " << x << "  MassVel " << MassVel;
//     prot << "Jf = " << jF * 0.3048 << " m/s" << " Jg = " << jG * 0.3048 << " m/s\n";
// prot << " alpha " << v1 << " dens " << dens;
//   if(v1< x) prot << v1;
	return dens;
} /* density_chexal */

/* PCLV */
int chexal_pclv(double x, double jF, double jG, double& v1, double& c01,
	double& vgj1, double& v2, double& c02, double& vgj2) {

	double gc0[2], gvgj[2], sign, gvoid[2];
	int error;

	/*     The  PCLV  subroutine  has  been  classified  by  EPRI  as  a */
	/*     Research Program. As such, it has not been developed and tested */
	/*     to the extent that a production program would be, and unforeseen */
	/*     results may occur when running the program.  EPRI does not make */
	/*     any warranty or representation whatsoever, expressed or implied */
	/*     including any warranty of  merchantability or fitness for any */
	/*     purpose with respect to the program.  Nor does EPRI assume any */
	/*     liability whatsoever with respect to any use of the program or */
	/*     any portion thereof or with  respect to any damages which may */
	/*     result from such use. */

	/*     This subroutine is based on the 1996 Revision of the EPRI */
	/*     Chexal-Lellouche Void Fraction model as documented in */
	/*     EPRI Report TR-106326, Void Fraction Technology for Design */
	/*     and Analysis, December 1996. */

	/*     Subroutine calling arguments are: */
	/*     Input - */
	/*       DH   - Hydraulic Diameter (ft) */
	/*       jF  - Superficial Liquid Velocity (ft/sec) */
	/*       jG  - Superficial Vapor Velocity (ft/sec) */
	/*       ORT  - Flow Direction in degrees from vertical - (0 to 90 deg) */
	/*       V1   - Void Fraction */
	/*       C01  - Co, Concentration Parameter */
	/*       VGJ1 - Vgj, Drift Velocity (same units as input) */
	/*     For Counter-Current Flows */
	/*       V2   - 2nd Void Fraction */
	/*       C02  - Co, Concentration Parameter for 2nd void fraction */
	/*       VGJ2 - Vgj, Drift Velocity for 2nd void fraction (m/s or ft/sec) */
	/*       error flag ier is function return value   */
	/*       IER  - Error flag, set to nonzero value if trouble is encountered */
	/*          IER = 1 Cannot find 1st void fraction */
	/*          IER = 2 Cannot find 2nd void fraction, CC Flow only */
	/*          IER = 3 Cannot find either void fraction */
	/*          IER = 4 Bad CounterCurrent Flow Pair jF > jFSTAR */
	/*          IER = 5 CCFL Solution in FRESCO did not converge */
	/*                   in allotted number of iterations */
	// ier is used as return value

	/*     All calling argument parameters are double precision variables */
	/*     All internal calculations are performed */
	/*     in British units. Conversion is done in calling function  */
	/*     Property values of water and steam must be known */
	/*     before the subroutine is called.  These are */
	/*       MF   - liquid viscosity (lbm/ft-hr) */
	/*       MG   - vapor viscosity (lbm/ft-hr) */
	/*       P    - pressure (psia) */
	/*       PC   - critical pressure (psia) */
	/*       rhoF - liquid density (lbm/ft3) */
	/*       rhoG - vapor density (lbm/ft3) */
	/*       SIG  - surface tension (lbf/ft) */
	/*    Note, all properties are in British units. */

/* already covered in calling function
 *only values above 1e-7 allowed
	if (fabs(jF) < 1e-7) {
		if (jF < 0. || (jF == 0. && jG < 0.)) {
			jF = -1e-7;
		} else {
			jF = 1e-7;
		}
	}
	if (fabs(jG) < 1e-7) {
		if (jG < 0. || (jG == 0. && jF < 0.)) {
			jG = -1e-7;
		} else {
			jG = 1e-7;
		}
	}
*/
	sign = 1.;
	if (fra == 0.) { // horizontal
		if (jF * jG < 0. && jG < 0.) {
			sign = -1.;
			jF = -jF;
			jG = -jG;
		}
	}

	error = chexal_glenns(x, jF, jG, gvoid, gc0, gvgj);

	/*     SETUP RETURN ARGUMENTS
	 *  conversion to SI Units will be done in calling function */
	v1 = gvoid[0];
	c01 = gc0[0];
	v2 = gvoid[1];
	c02 = gc0[1];
	vgj1 = gvgj[0] * sign;
	vgj2 = gvgj[1] * sign;
	jFstar *= sign;
	return error;
} /* pclv_ */

/* GLENNS */
int chexal_glenns(double x, double jF, double jG, double* gvoid, double* gc0,
	double* gvgj) {
	/* Local variables */
	int i, error = 0;
	double R, b1, c3, k0, c1x, jSum,
		ReF, ReG, vgj0, vhomogen;

	/* Function Body */
	c3 = 0.;

	ReF = ReFxx * jF;
	ReG = ReGxx * jG;

	if (ReF * ReG < 0.) { /*  COUNTERCURRENT FLOW */
		//      if (ionce == 0) {
		error = chexal_fresco(x, jF, jG, true, &jFstar, &alstar);
		if (fabs(jF / jFstar) > 1.) {
			if (Base.showDPTube) {
				cout << "\nTHE INPUT jF AND jG PAIR VIOLATES CCFL";
				cout << "\njF = " << jF * .3048 << " m/s IS GREATER THAN jFSTAR = " << jFstar * .3048 << " m/s";
				prot << "\nTHE INPUT jF AND jG PAIR VIOLATES CCFL";
							prot << "\njF = " << jF << " FT/S IS GREATER THAN jFSTAR = " << jFstar << " FT/S";
			}
			for (i = 0; i <= 1; ++i) {
				gvoid[i] = 1e12;
				gc0[i] = 0.;
				gvgj[i] = 0.;
			}
			return 4;
		}
	}
	ReF = ReFxx * jF;
	ReG = ReGxx * jG;
	c3 = chexal_c3func(ReF, ReG, jF, jFstar);
	vgj0 = vj0xx * c3;
	jSum = jF + jG;

	chexal_c0func_init(ReF, ReG, c1x, b1, k0, R);

	error = chexal_iterat(x, jF, jG, gvoid, false, c3, gc0, vgj0, c1x,
		gvgj, k0, R);

	if (ReF * ReG > 0.) {
		vhomogen = jG / jSum;
		for (i = 0; i <= 1; ++i) {
			if (gvoid[i] > .999 && gvoid[i] < vhomogen) {
				gvoid[i] = vhomogen;
				gc0[i] = 1.;
				gvgj[i] = 0.;
			}
		}
	}
	return error;
} /* glenns_ */


double chexal_c3func(double ReF, double ReG, double jF, double jF2) {

	/* Local variables */
	double d1, temp, xReG, jFrx,
		term2, c3vert, c3horz;


	/*     C3FUNC COMPUTES THE CHEXAL LELLOUCHE C3 COEFFICIENT */
	/*     C3FUNC IS CALLED FROM SUBROUTINES EPRIDV AND FRESCO */

	/*      REF     = superficial liquid Reynolds number */
	/*      REG     = superficial vapor Reynolds number */
	/*      jF     = liquid superficial velocity */
	/*      jF2 = CCFL liquid superficial velocity (jFstar) or liquid superficial velocity (jF)
					  depends on calling function*/
	/*      DH      = hydraulic diameter (feet) */
	/*      FRA     =  mixing fraction based on orientation */

	/*     INPUT UNITS MUST BE CONSISTENT AND DOUBLE PRECISION */

	d1 = .125;
	c3horz = 0.;
	c3vert = 0.;

	/*     VERTICAL ------------------------------------------------ */
	if (fabs(fra) > 1e-7) {
		if (ReF >= 0. && ReG >= 0.) {
			/* COCURRENT UP */
			//temp = (ReF / 3e5);
			//if (temp < 85.) {
			//	c3vert = exp(-temp) * 2.;
			//}
			//if (c3vert < .5) {
			//	c3vert = .5;
			//}
			if (ReF < 415888.3083) {
				//0.5 = exp(-ReF/ 3e5) * 2. ; exp(-ReF/3e5) = 0.25 -> ReF = 3e5 * ln(0.25) = 415888.3083 
				c3vert = exp(-ReF / 3e5) * 2.;
			}
			else {
				c3vert = .5;
			}
		}
		else {
			/* COCURRENT DOWN AND COUNTERCURRENT */
			double xReF = -(ReF);
			double xd = d1 / dh;
			if (ReG <= 0.) {
				xReG = 0.;
				jFrx = 1.;
			}
			else {
				xReG = ReG;
				if (jFstar < 0.) {
					temp = fabs(jF / jF2);
					double z = .8;
					if (temp > .3) {
						z -= (temp - .3);
					}
					jFrx = pow(1. - temp, z);
					if (jFrx > 1.) {
						jFrx = 1.;
					}
					else if (jFrx < 0.) {
						jFrx = 0.;
					}
				}
				else {
					jFrx = 1.;
				}
			}
			double ajFrx = 1. - jFrx;

			if (fabs(xReG) < 10. / 85) {
				temp = 0.;
			}
			else {
				temp = exp(-10. / fabs(xReG));
			}
			double term1 = exp(pow((xReF + xReG * jFrx * temp * pow(8., 0.5 * xd)) / 3.5e5, 0.4)) * 2.;
			temp = xReF / (jFrx * 3.5e4 + 2.5e4) * xd * xd;
			if (temp < 85.) {
				term2 = pow(xReF, 0.035) * -1.7 * exp(-temp);
			}
			else {
				term2 = 0.;
			}
			double term3 = (jFrx * .26 + ajFrx * .85) * pow(xd, 0.1) * pow(xReF, 0.001);
			double c10 = term1 + term2 + term3;

			if (c10 > 0.) {
				double b2 = pow(1. / (xReF * .05 / 3.5e5 + 1.), 0.4);
				c3vert = pow(c10 / 2., b2) * 2.;
			}
			if (ReG > 0.) { /* THE LARGE DH MODEL IS BASED */
				if (dh > d1) { /* ON THE CSA SIMULATIONS OF */
					if (dh < 1.) { /* THE GE 1 AND 4 FT TESTS */
						c3vert = ((dh - d1) * .6 / (1. - d1) + c3vert *
							(1. - (dh - d1) / (1. - d1))) * jFrx + c3vert * ajFrx;
					}
					else if (dh < 3.) {
						c3vert = (.6 - (dh - 1.) * .27) * jFrx + c3vert * ajFrx;
					}
					else {
						c3vert = jFrx * .06 + c3vert * ajFrx;
					}
				}
			}
		}
		c3vert *= fra;
	}
	/*     HORIZONTAL ------------------------------------------------ */
	if (fra != 1.) {
		if (ReF * ReG >= 0.) { /*  COCURRENT */
			//temp = fabs(ReF / 3e5);
			//if (temp < 85.) {
			//	c3horz = exp(-temp) * .5;
			//}
			//if (c3horz < .125) {
			//	c3horz = .125;
			//}
			if (ReF < 415888.3083) {
				//0.125 = exp(-ReF/ 3e5) * 0.5 ; exp(-ReF/3e5) = 0.25 -> re = 3e5 * ln(0.25) = 415888.3083 
				c3horz = exp(-ReF / 3e5) * 0.5;
			}
			else {
				c3horz = .125;
			}
		}
		else { /* COUNTERCURRENT */
			temp = fabs(ReF / 600.);
			if (temp < 85.) {
				c3horz = pow(fabs(ReF / 4.5e5), 0.19) * 4.5 + exp(-temp) * .5;
			}
			else {
				c3horz = pow(fabs(ReF / 4.5e5), 0.19) * 4.5;
			}
		}

		if (ReG < 0.) {
			c3horz *= (fra - 1.);
		}
		else {
			c3horz *= (1. - fra);
		}
	}

	return fmin(10., c3vert + c3horz);
} /* c3func_ */

int chexal_c0func_init(double ReF, double ReG,
	double& c1x, double& b1, double& k0, double& R) {

	/* Local variables */
	double re, b1horz, b1vert, c1xhorz, c1xvert;


	/*     Subroutine C0FUNC_init calculates several constants for C0FUNC in the  */
	/*     Chexal-Lellouche drift flux correlation.  The parameters computed */
	/*     here are orientation dependent. */
	/*     This subroutine is called from EPRIDV and FRESCO */

	/*    Input Parameters */
	/*      ReF     = superficial liquid Reynolds number */
	/*      ReG     = superficial vapor Reynolds number */

	/*     Output Parameters */
	/*      C0      = distribution parameter */
	/*      C1X     = exponent used for (1-alpha) for C1 constant */
	/*      B1      = B1 constant */
	/*      K0      = parameter used in defining fluid parameter L */
	/*      R       = r parameter from the correlation */


	/*  EXPONENTIAL UNDERFLOW LIMIT is set to -85.*/
	/*       VERTICAL -------------------------------------------------- */
	if (fabs(fra) < 1e-7) {
		b1vert = 0.;
		c1xvert = 0.;
	}
	else {
		re = ReG;
		if (ReF >= ReG && ReG > 0.) {
			re = ReF;
		}
		if (re < 0.) {
			b1vert = 1e-10;
		}
		else {
			//b1vert = 1.;
			//if (fabs(re / 6e4) < 85.) {
			//	b1vert = 1. / (exp(-re / 6e4) + 1.);
			//}
			//if (b1vert > .8) {
			//	b1vert = .8;
			//}
			if (re < 83177.66166) {
				//0.8 = 1./ (0.25 + 1) ; exp(-re/6e4) = 0.25 -> re = 6e4 * ln(0.25) = 83177.66166 
				b1vert = 1. / (exp(-re / 6e4) + 1.);
			}
			else {
				b1vert = .8;
			}
		}
		if (ReG >= 0.) {
			c1xvert = b1vert;
		}
		else {
			c1xvert = .5;
		}
		b1vert *= fra;
		c1xvert *= fra;
	}
	/*       HORIZONTAL ------------------------------------------------ */
	if (fabs(fra - 1.) < 1e-7) {
		b1horz = 0.;
		c1xhorz = 0.;
	}
	else {
		re = fabs(ReG);
		if (ReF * ReG >= 0.) {
			if (fabs(ReF) >= fabs(ReG)) {
				re = fabs(ReF);
			}
		}
		//b1horz = 1.;
		//if (re / 6e4 < 85.) {
		//	b1horz = 1. / (exp(-re / 6e4) + 1.);
		//}
		//if (b1horz > .8) {
		//	b1horz = .8;
		//}
		if (re < 83177.66166) {
			//0.8 = 1./ (0.25 + 1) ; exp(-re/6e4) = 0.25 -> re = 6e4 * ln(0.25) = 83177.66166 
			b1horz = 1. / (exp(-re / 6e4) + 1.);
		}
		else {
			b1horz = .8;
		}
		b1horz *= (1. - fra);
		c1xhorz = b1horz;
	}
	c1x = c1xvert + c1xhorz;
	b1 = b1vert + b1horz;
	k0 = b1 + (1. - b1) / pow(rr, 0.25);
	R = (1.57 / rr + 1.) / (1. - b1);

	return 0;
} /* c0func_init */

double chexal_c0func(double jSum, double jG, double vgj0, double al,
	double alm1,  double k0, double R) {

	/* Local variables */
	double l, lh, lv, c0a, c0, denom;


	/*     Subroutine C0FUNC calculates several constants in the Chexal- */
	/*     Lellouche drift flux correlation.  The parameters computed */
	/*     here are orientation dependent. */
	/*     This subroutine is called from EPRIDV and FRESCO */

	/*    Input Parameters */
	/*      jSum     = sum of superficial velocities */
	/*      jG     = vapor superficial velocity */
	/*      VGJ0    = VGJ/C1 where C1 = ALM1**C1X */
	/*      AL      = vapor  volume fraction */
	/*      ALM1    = liquid volume fraction */

	/*     return value */
	/*      C0      = distribution parameter */


	/*  EXPONENTIAL UNDERFLOW LIMIT is set to -85.*/
	lv = 1.;
	denom = 1.;
	if (cp < 85.) {
		denom = 1. - exp(-(cp));
	}
	if (cp * al < 85.) {
		lv = (1. - exp(-(cp)*al)) / denom;
	}
	lh = lv;

	lh *= pow(al, 0.05) * (alm1 * alm1) + 1.;
	l = lv * fra + lh * (1. - fra);
	c0 = l / (k0 + (1. - k0) * pow(al, R));

	if (jG < 0.) {
		c0a = fmin(6., -(vgj0)*pow(alm1, 0.2) / jSum);
		c0 = fmax(c0, c0a);
	}

	return c0;
} /* c0func_ */

int chexal_fresco(double x, double jF, double jG, bool isCalcjFstar, 
	 double* answer, double* Void) {
	/* Initialized data */

	double frac = .1; /*  INITIAL DELTA FRACTION OF jF OR jG */

	/* System generated locals */
	int ier;

	/* Local variables */
	double  R, b1, c3 = 0., k0, gc0[2], c1x, ReF = 0., ReG = 0.;
	int key, inum, iboth, ihist;
	double  gvgj[2], gvoid[2]; //, b1ccfl, c0ccfl, c3ccfl, vgjccfl;
	double delta, jFold, jGold;

	/*     Subroutine FRESCO is used to compute a ccfl value of void */
	/*     and jF (jFstar) on the CCFL line given a jG as input.  Fresco */
	/*     is called from EPRIDV. */

	/*    Input Parameters */
	/*      jF = value of jF (MUST BE NEGATIVE) */
	/*      jG = value of jG (MUST BE POSITIVE) */
	/*      isCalcjFstar true if jFstar is computed */
	/*                false if jGstar is computed */
	/*      isCalcjFstar = true,  XjFIN IS THE INITIAL GUESS FOR THE ANSWER (jFstar) */
	/*      isCalcjFstar = false, XjGIN IS THE INITIAL GUESS FOR THE ANSWER (jGstar) */

	/*      REFXX   = partial liquid Reynolds number (REf/jF) */
	/*      REGXX   = partial vapor Reynolds number (REg/jG) */
	/*      VJ0XX   = non flow and alpha part of VGJ */
	/*      DH      = hydraulic diameter */
	/*      FRA     =  mixing fraction based on orientation */
	/*      RR      = RHOF/RHOG */
	/*      CP      = 4*Pcrit**2/P*(Pcrit-P) */
	/*      C0OLD   = prior time step value of C0 */
	/*      C0TAU   = C0 smoothing constant */
	/*      TSTEP   = time step */

	/*   Output parameters */
	/*      ANSWER = jFstar or jGstar */
	/*      VOID   = void fraction on the CCFL line */

	* Void = 0.;
	ihist = 0;
	iboth = 0;
	jFold = jF;
	jGold = jG;
	if (isCalcjFstar) {
		if (jF >= 0.) {
			jF = -1e4;
		}
		delta = -frac * jF;
		ReG = ReGxx * jG;
	}
	else {
		delta = frac * jG;
		ReF = ReFxx * jF;
		ReG = ReGxx * jG;
		c3 = chexal_c3func(ReF, ReG, jF, jF);
	}

	for (inum = 1; inum <= 200; ++inum) {
		if (isCalcjFstar) {
			ReF = ReFxx * jF;
			c3 = chexal_c3func(ReF, ReG, jF, jF);
		}
		else {
			ReG = ReGxx * jG;
		}

		chexal_c0func_init(ReF, ReG, c1x, b1, k0, R);

		ier = chexal_iterat(x, jF, jG, gvoid, true, 
			c3, gc0, vj0xx * c3, c1x, gvgj, k0, R);

		if (gvoid[0] >= 2.) {
			key = 2;
		}
		else {
			key = 1;
			*Void = gvoid[0];
			//         b1ccfl = b1; /*  CCFL */
			//         c3ccfl = c3; /*  CCFL */
			//         c0ccfl = gc0[0]; /*  CCFL */
			//         vgjccfl = gvgj[0]; /*  CCFL */
		}
		//		prot<< "\n inum " << inum << " jF " << jF << " jG " << jG << " gvoid " << gvoid[0] << " jfOld " << jFold << " jGold " << jGold;
		chexal_march2(isCalcjFstar, key, jF, jG, jFold, jGold, delta, ihist, iboth);

		if (delta > 1e10) {
			break;
		}
	}

	if (inum >= 200) {
		cout << "BAD NEWS: DID 200 ITERATIONS AND COULD NOT CONVERGE ON A CCFL SOLUTION" << endl;
		cout << "SUBSEQUENT RESULTS FOR THIS PAIR MAY BE CLOSE OR THEY MAY BE TRASH" << endl;
		//		prot << "\nBAD NEWS: DID 200 ITERATIONS AND COULD NOT CONVERGE ON A CCFL SOLUTION";
		//		prot << "\nSUBSEQUENT RESULTS FOR THIS PAIR MAY BE CLOSE OR THEY MAY BE TRASH" << endl;

		ier = 5;
	}
	if (isCalcjFstar) {
		*answer = (jF + jFold) * .5;
		prot << "\n jFstar " << *answer << " x " << x;
	}
	else {
		*answer = (jG + jGold) * .5;
	}

	return ier;
} /* fresco_ */

int chexal_iterat(double x, double jF, double jG,
	double* gvoid, bool iccfl, double& c3, double* gc0,
	double vgj0, double c1x, double* gvgj,
	double k0, double R) {
	/* Local variables */
	int i, ier = 0;
	double C0;
	int np;
	double a, b, fa, fb, alpha, vgj, tol;
	double ffmin, zlmin, direct, vhomogen;
	double fx, fh, f1;
	/*     ITERAT USES THE REGULA-FALSI METHOD (Anderson-Bjoerk) for co-current flow */
	/*     TO SOLVE THE CHEXAL-LELLOUCHE CORRELATION FOR jF OR jG AND VOID */
	/*     FRACTION. For counter-current flow the stepped method to establish to boundaries
	 *     is kept to find the 2 solutions.
	 *
	 *  ITERAT is called from subroutine FRESCO (iccfl = 1) and Glenns (iccfl = 0). */

	 /*      jF    = superficial liquid velocity   (INPUT) */
	 /*      jG    = superficial vapor  velocity   (INPUT) */
	 /*      GVOID = void fraction                        (OUTPUT) */
	 /*      C3    = C3 constant                   (INPUT/output) */
	 /*      GC0      = distribution parameter            (OUTPUT) */
	 /*      VGJ0    = VGJ/C1 where C1 = ALM1**C1X (INPUT) */
	 /*      CP      = 4*Pcrit**2/P*(Pcrit-P)      (INPUT) */
	 /*      C1X     = exponent used for (1-alpha) for C1 constant  (INPUT) */
	 /*      GVGJ     = drift velocity                    (OUTPUT) */
	 /*      FRA     =  mixing fraction based on orientation (INPUT) */
	 /*      K0      = parameter used in defining fluid parameter L (INPUT) */
	 /*      R       = r parameter from the correlation   (INPUT) */
	 /*      REF     = superficial liquid Reynolds number (INPUT) */
	 /*      REG     = superficial vapor Reynolds number   (INPUT) */

	 /* Function Body */
 //   ier1 = 0;

	if (jF == 0. && jG < 0.) {
		gvoid[0] = 1.;
		gvgj[0] = 0.;
		gc0[0] = 1.;
		c3 = 2.;
		return 0;
	}
	double jSum = jF + jG;

	if ((jF * jG > 0.) && (!iccfl)) { //co-current flow
		// establish range 
		vhomogen = jG / jSum;

		if (jF > 0.) {
			// upward or horizontal flow
			// solution is only possible between x and homogeneous void fraction
			a = x;
			fa = chexal_poly(a, jSum, jG, vgj0, c1x, k0, R, vgj, C0);
			if (fa > 0.) { // solution void will be lower than x this happens on low mass velocities and high x -> most probably annular (churn) flow 
				gvoid[0] =vhomogen; 
				gvgj[0] = 0.;
				gc0[0] = 1.;
				return 0;
         }
			fx = chexal_poly(a, jSum, jG, vgj0, c1x, k0, R, vgj, C0);
			b = vhomogen;
			fb = chexal_poly(b, jSum, jG, vgj0, c1x, k0, R, vgj, C0);
			fh = chexal_poly(b, jSum, jG, vgj0, c1x, k0, R,  vgj, C0);
		}
		else {
			//downward flow
			// solution between homogeneous void fraction and 1
			a = x;
			fa = chexal_poly(a, jSum, jG, vgj0, c1x, k0, R,  vgj, C0);
			fh = chexal_poly(a, jSum, jG, vgj0, c1x, k0, R,  vgj, C0);
			b = .9999999; // close to 1
			f1 = chexal_poly(b, jSum, jG, vgj0, c1x, k0, R,  vgj, C0);
			fb = chexal_poly(b, jSum, jG, vgj0, c1x, k0, R,  vgj, C0);
		}
		/*     RANGE ESTABLISHED: NOW USE REGULA FALSI TO FIND SOLUTION */

		alpha = chexal_regula_falsi(a, fa, b, fb, jSum, jG, vgj0, c1x,
			k0, R,  vgj, C0);

		gvoid[0] = alpha;
		gvgj[0] = vgj;
		gc0[0] = C0;

	}
	else { // counter-current flow
		np = 0;
		if (!iccfl) {
			if (jF * jG < 0.) {
				np = 1;
			}
		}
		for (i = 0; i <= np; ++i) {
			a = 0.;
			direct = 1.;
			b = a;

			if (i == 1) {
				b = 1.;
				a = b;
				direct = -1.;
			}

			fa = chexal_poly(a, jSum, jG, vgj0, c1x, k0, R, vgj, C0);

			ffmin = fabs(fa);
			zlmin = a;

			while (true) {
				/*    COARSE OUTER LOOP TO ESTABLISH RANGE ----------- */
				if (b >= .015 && b <= .985) {
					tol = direct * .01;
				}
				else if (b < .001 || b > .999) {
					tol = direct * 1e-4;
				}
				else { // is only remaining if (zl < .015 ||zl > .985) {
					tol = direct * .001;
				}

				b = a + tol;
				fb = chexal_poly(b, jSum, jG, vgj0, c1x, k0, R, vgj, C0);
				//				prot << "\n i " << i << " np "<<np<<" zll " << zll << " ff2 " << ff2;
				if (fa * fb < 0.) {
					break;
				}
				if (fabs(fb) < fabs(ffmin)) {
					ffmin = fb;
					zlmin = b;
				}
				fa = fb;
				a = b;

				if (i == 0) {
					if (b > 1e-7 && b < .9999999) {
						continue;
					}
					else {
						prot << "\nb " << b;
						if (iccfl) {
							gvoid[0] = 10.;
							return 0;
						}
						if (ier != -1) {
							cout << "\n** SORRY: NO SOLUTION FOUND FOR ROOT " << i;
							cout << "\nTHE MINIMUM ERROR OF " << ffmin << " OCCURRED AT A VOID FRACTION OF " << zlmin <<
								"\nEXAMINE THE RESULTS CAREFULLY" << endl;
							//							prot << "\n** SORRY: NO SOLUTION FOUND FOR ROOT " << i;
							//							prot << "\nTHE MINIMUM ERROR OF " << ffmin << " OCCURRED AT A VOID FRACTION OF " << zlmin <<
							//								"\nEXAMINE THE RESULTS CAREFULLY";
						}
						//						ier = 1;
						if (np == 1) {
							return 1;
						}
					}
				}
				else {
					if (b >= .9999999) {
						if (ier != -1) {
							cout << "\n** SORRY: NO SOLUTION FOUND FOR ROOT " << i;
							cout << "\nTHE MINIMUM ERROR OF " << ffmin << " OCCURRED AT A VOID FRACTION OF " << zlmin <<
								"\nEXAMINE THE RESULTS CAREFULLY";
						}
						ier = 2;
					}
					else if (b > 1e-7) {
						continue;
					}
					else {
						b = 0.;
						a = b;
						direct = 1.;
						continue;
					}
				}
			}
			/*     RANGE ESTABLISHED: NOW USE REGULA FALSI TO FIND SOLUTION */
	//		zr = zl + tol;

			if ((a > 1e-7 && a < .9999999) && (b > 1e-7 && b < .9999999)) {
				alpha = chexal_regula_falsi(a, fa, b, fb, jSum, jG, vgj0, c1x,
					k0, R, vgj, C0);
			}
			else {
				a = b;
				fa = chexal_poly(a, jSum, jG, vgj0, c1x, k0, R,  vgj, C0);
			}

			gvoid[i] = (a + b) * .5;
			gvgj[i] = vgj;
			gc0[i] = C0;
		}
	}
	return ier;
} /* iterat_ */


double chexal_regula_falsi(double a, double fa, double b, double fb, double jSum, double jG, double vgj0, double c1x,
	double k0, double R,  double& vgj, double& C0) {
	double c = 0.;
	double fc, abj;
	for (int i = 1; i <= 100; ++i) {
		if (fabs(b - a) < 1e-10) {
			c = b;
			break;
		}
		c = a - fa * (b - a) / (fb - fa);

		fc = chexal_poly(c, jSum, jG, vgj0, c1x, k0, R,  vgj, C0);
		if (fabs(fc) < 1e-14) {
			break;
		}
		if (ISIGN(fc) * ISIGN(fb) < 0) {
			a = b;
			fa = fb;
		}
		else {
			abj = 1. - fc / fb;
			if (abj <= 1e-6) {
				abj = .5;
			}
			fa *= abj;
		}
		b = c;
		fb = fc;
	}
	return c;
}


double chexal_poly(double alpha, double jSum, double jG,
	double vgj0, double c1x, double k0, double R,
	 double& vgj, double& C0) {
	/* Local variables */
	double diff, alm1;

	/*     POLY CALCULATES THE DIFFERENCE BETWEEN THE jG CALCULATED FROM */
	/*     THE VOID FRACTION AND THE ACTUAL jG, THIS DIFFERENCE IS USED */
	/*     BY THE ITERAT ROUTINE.   POLY IS CALLED BY SUBROUTINE ITERAT. */

	/*      ALPHA   = void fraction                  (INPUT) */
	/*      jSum    = sum of superficial velocities  (INPUT) */
	/*      jG      = superficial vapor velocity     (INPUT) */
	/*      VGJ     = drift velocity                         (OUTPUT) */
	/*      VGJ0    = VGJ/C1 where C1 = ALM1**C1X    (INPUT) */
	/*      C1X     = exponent used for (1-alpha) for C1 constant  (INPUT) */
	/*      CP      = 4*Pcrit**2/P*(Pcrit-P)         (INPUT) */
	/*      C0      = drift flux distribution parameter      (OUTPUT) */
	/*      K0      = parameter used in defining fluid parameter L (INPUT) */
	/*      R       = r parameter from the correlation      (INPUT) */
	/*      FRA     =  mixing fraction based on orientation (INPUT) */


	if (alpha >= .9999999) alpha = .9999999;
	if (alpha <= 1e-7) alpha = 1e-7;
	alm1 = 1. - alpha;
	vgj = vgj0 * pow(alm1, c1x);

	C0 = chexal_c0func(jSum, jG, vgj0, alpha, alm1,  k0, R);

	diff = alpha * (C0 * jSum + vgj) - jG; /*  CALCULATE ERROR AND RETURN */
	return diff;
} /* poly_ */

int chexal_march2(bool isCalcjFstar, int key, double& jF, double& jG,
	double& jFold, double& jGold, double& delta,
	int& ihist, int& iboth) {
	/*      MARCH2 Checks convergence of a solution on the CCFL line */
	/*      and increments jG or jF for the next try at the solution */
	/*      of the Chexal-Lellouche drift flux.  MARCH2 is called from */
	/*      subroutine FRESCO. */

	/*      isCalcjFstar = true if jFstar is computed */
	/*                  false if jGstar is computed */
	/*      KEY   = 1, converged solution found      (input)*/
	/*              2, converged solution not found */
	/* input or output (c++ references)*/
	/*      jF    = superficial liquid velocity      (INPUT/OUTPUT) */
	/*      jG    = superficial vapor velocity       (INPUT/OUTPUT) */
	/*      jFOLD = old superficial liquid velocity  (INPUT/OUTPUT) */
	/*      jGOLD = old superficial vapor velocity   (INPUT/OUTPUT) */
	/*      DELTA = change in jF or jG for the next iteration (OUTPUT) */
	/*      IHIST = 1, if solution found for last try         (OUTPUT) */
	/*              2, if solution NOT found for the last try */
	/*      IBOTH = 0, last two values of jF or jG produced   (OUTPUT) */
	/*                 were on one side of the CCFL line */
	/*              2, values of jF or jG have been found */
	/*                 on both sides of the CCFL line. */
	/* Initialized data */

	double convrg = 1e12; /*  RETURN VALUE INDICATING CONVERgence*/
	double faster = 1.2; /*  ACCELERATION FACTOR ON DELTA */
	double delmax = 1e4; /*  MAXIMUM ALLOWABLE DELTA */

	/* Local variables */
	static double jFbad, jGbad, jFgood, jGgood;
	double enough = 1e-6;

	if (key == 1) { /* ---- LAST TRY GOOD */
		if (ihist == 2) {
			iboth = 2;
		}
		if (!isCalcjFstar) { // jGStar is calculated
			jGgood = jG;
			if (ihist == 2 || iboth >= 2) {
				if (fabs(jG - jGold) / jGold < enough) {
					delta = convrg;
					return 0;
				}
				jGold = jG;
				jG = (jG + jGbad) * .5;
			}
			else {
				jGold = jG;
				//				delta = faster * delta;
				//				if (delta > delmax) {
				//					delta = delmax;
				//				}
				delta = fmin(faster * delta, delmax);
				jG += delta;
			}
		}
		else { // jFstar is calculated
			jFgood = jF;
			if (ihist == 2 || iboth >= 2) {
				if (fabs(jF - jFold) / fabs(jFold) < enough) {
					delta = convrg;
					return 0;
				}
				jFold = jF;
				jF = (jF + jFbad) * .5;
			}
			else {
				jFold = jF;
				//				delta = faster * delta;
				//				if (delta > delmax) {
				//					delta = delmax;
				//				}
				delta = fmin(faster * delta, delmax);
				jF -= delta;
			}
		}
		ihist = 1;
	}
	else { /* ---- LAST TRY BAD */
		if (ihist == 1) {
			iboth = 2;
		}
		if (!isCalcjFstar) { // jGStar is calculated
			jGbad = jG;
			if (ihist == 1 || iboth >= 2) {
				if (fabs(jG - jGold) / jGold < enough) {
					delta = convrg;
					return 0;
				}
				jGold = jG;
				jG = (jG + jGgood) * .5;
			}
			else {
				jGold = jG;
				//				delta = faster * delta;
				//				if (delta > delmax) {
				//					delta = delmax;
				//				}
				delta = fmin(faster * delta, delmax);
				if (jG - delta <= 0.) {
					delta = 0.99 * jG;
				}
				jG -= delta;
			}
		}
		else { // jFStar is calculated
			jFbad = jF;
			if (ihist == 1 || iboth >= 2) {
				if (fabs((jF - jFold) / jFold) < enough) {
					delta = convrg;
					return 0;
				}
				jFold = jF;
				jF = (jF + jFgood) * .5;
			}
			else {
				jFold = jF;
				//				delta = faster * delta;
				//				if (delta > delmax) {
				//					delta = delmax;
				//				}
				delta = fmin(faster * delta, delmax);
				if (jF + delta >= 0.) {
					delta = -0.99 * jF;
				}
				jF += delta;
			}
		}
		ihist = 2;
	}
	return 0;
} /* march2 */

