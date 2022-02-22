//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

//common variables used in Chexal functions

//static double jFstar,  ///<superficial fluid velocity at CCFL(ft/sec)
//					jGstar, ///<superficial gas velocity at CCFL(ft/sec)
//					alstar; ///<void fraction at CCFL 

//Properties in British Units
static double rhoF,///< density of fluid (water) (lbm/ft^3)
rhoG,///<density of gas (vapor,steam) (lbm/ft^3)
sig,///<surface tension (lbf/ft) 
rr,///<density ratio (-)
cp;///<pressure factor (-)

// constant values set once
static double dh, ///<hydraulic diameter (ft) 
ReF,///<liquid Reynolds number  
ReG,///<vapor Reynolds number  
fra,///<tube orientation factor (-), tube orientation between 0 deg(vertical) and 90 deg (horizontal)  
x,///<steam quality [-]
jF,///<superficial water (Fluid) velocity [ft/s]
jG,///<superficial steam (Gas) velocity [ft/s]
jSum,///<sum of jF and jG [ft/s]
VoidHomogeneous,///<homogeneous void fraction [-]
c1x,///<constant 
k0,///<constant 
R, ///<constant 
vgj0; ///<non alpha part of VGJ 

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

 *     This subroutine is based on the 1996 Revision of the EPRI\n
 *     Chexal-Lellouche Void Fraction model as documented in\n
 *     EPRI Report TR-106326, Void Fraction Technology for Design\n
 *     and Analysis, December 1996.\n
 *     -end quote-\n

 *     as counter-current flows, i.e. steam against water and CCFL situations are unlikely\n
 *     within natural circulation boilers those parts and references are removed from code
  */

  /**
	* @brief chexal_c3func computes the Chexal-Lellouche C3 coefficient
	*
	* Chexal_c3func is called once
	*
	* @return double C3
	*/
double chexal_c3func();

/**
 * @brief chexal_c0func_init calculates several constants for C0FUNC in the Chexal-Lellouche drift flux correlation.
 *
 * The parameters computed here are orientation dependent.
 * This function is called once
 *
 * @return int error code
 */
int chexal_c0func_init();

/**
 * @brief  C0FUNC calculates several constants in the Chexal-Lellouche drift flux correlation.
 *
 * The parameters computed here are orientation dependent.
 *
 * This function is called from ::chexal_poly
 *
 * @param [in] al vapor volume fraction
 * @param [in] alm1 liquid volume fraction

 * @return double C0  distribution parameter
 */
double chexal_c0func(double al, double alm1);

/**
 * @brief iteration of void fraction
 *
 * chexal_iterat solves the Chexal-Lellouche correlation for jF or jG and void fraction\n
 * It uses the regula-falsi method (Anderson-Bjoerk) for co-current flow\n
 *
 * chexal_iterat is called once
 * @param [out] Void void fraction
 * @param [out] C0 distribution parameter
 * @param [out] vgj drift velocity
 * @return int
 */
int chexal_iterat(double& Void, double& C0, double& vgj);


/**
 * \brief  Chexal_poly calculates the difference between the jG calculated from the void fraction and the actual jg
 *
 * this difference is used by the chexal_iterat function
 *
 * chexal_poly is called by function ::chexal_iterat
 *
 * @param [in] al void fraction
 * @param [out] vgj drift velocity
 * @param [out] c0 drift flux distribution parameter
 * @return double
 */
double chexal_poly(double al, double& c0, double& vgj);

#define ISIGN(x) ((x) >= 0. ? 1 : -1)

// Properties are implicit at saturation (no ..Sat.. in variable name)

double _tube::Density_chexal(double xIn, double pMPa, double volW,
	double volS, double dynVisW, double dynVisS, double SurfTens, double& Void) {
	/* System generated locals */
	double dens; //mixture density [kg/m3]
	/* Local variables */
	static const double gravty = 32.17; // gravity constant, British Units
	int error = 0;
	double C0, vgj;

	double d__1;
	double c2, c4, c5, c7, c8, arg, vgj0xx, c3;
	x = xIn;
	/* if superficial velocity too low return with single phase */
	jF = (1. - x) * MassVel * volW / .3048; /* in ft/s */
	if (jF < 1e-7) { //no fluid, only steam
		Void = 1.;
		cout << "\n jF < 1e-7 shortcut steam";
		return 1. / volS;
	}
	jG = x * MassVel * volS / .3048; /* in ft/s */
	if (jG < 1e-7) { // no steam, only fluid
		Void = 0.;
		cout << "\njG < 1e-7 shortcut water";
		return 1. / volW;
	}

	/*  PROGRAM TO REDUCE ADIABATIC VOID FRACTION DATA USING THE CHEXAL- */
	/*  LELLOUCHE VOID FRACTION CORRELATION */

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
		if (rr > 153.550173) { //arg = 85; c5 = 85/86; -> rr = 153.550173 
			c5 = sqrt(150. / rr); //if rr >150 -> c5 < 1
			arg = c5 / (1. - c5);
			c2 = 1./(1. - exp(-arg));
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
	vgj0xx = c2 * c4 * 1.41 * pow(d__1, 0.25); /* vgj0xx = non flow and alpha part of VGJ */
	if (HeightSection < 0.) {
		jF = -jF;
		jG = -jG;
	}
	jSum = jF + jG;
	VoidHomogeneous = jG / jSum;
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
		ort = 90. - asin(fabs(HeightSection / LengthSection)) * 180. / 3.14159265;
		fra = pow(1. - ort / 90., 0.2);
	}
	double sign = 1.;
	if (fra == 0.) { // horizontal
		if (jF * jG < 0. && jG < 0.) {
			sign = -1.;
			jF = -jF;
			jG = -jG;
		}
	}

	ReF = rhoF * dh * 3600. / (dynVisW * 2419.1) * jF;
	ReG = rhoG * dh * 3600. / (dynVisS * 2419.1) * jG;
	if (Base.showDPTubeDetail) {
		prot << "\n           -------Steam-Water Results-------\n" <<
			" System Pressure - " << pMPa << " MPa";
		prot << "\n  Hydraulic Diameter = " << Dia << " m" <<
			"\n Flow Direction from Vertical = " << ort << " degrees" <<
			"\n x = " << x <<
			"\n l =  " << LengthSection << " h " << HeightSection;
	}

	c3 = chexal_c3func();
	vgj0 = vgj0xx * c3;

	chexal_c0func_init();

	error = chexal_iterat(Void, C0, vgj);

	if (Void > .999 && Void < VoidHomogeneous) {
		Void = VoidHomogeneous;
		C0 = 1.;
		vgj = 0.;
	}
	if (error > 0) {
		cout << "\n ERROR IN VOID SUBROUTINE, IER = " << error;
		if (Base.showDPTube) {
			prot << "\n ERROR IN VOID SUBROUTINE, IER = " << error << endl;
		}
	}

	dens = (1. - Void) / volW + Void / volS; // in metric units
	return dens;
} /* density_chexal */

double chexal_c3func() {

	/* Local variables */
	double d1, temp, xReG, jFrx,
		term2, c3vert, c3horz;

	/*     C3FUNC COMPUTES THE CHEXAL LELLOUCHE C3 COEFFICIENT */

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
			/* COCURRENT DOWN  */
			double xReF = -(ReF);
			double xd = d1 / dh;
			xReG = 0.; // cocurrent!!!
			jFrx = 1.; // cocurrent!!!
//			double ajFrx = 0.;

			if (fabs(xReF) < 10. / 85) {
				temp = 0.;
			}
			else {
				temp = exp(-10. / fabs(xReF));
			}
//			double term1 = exp(pow((xReF + xReG * jFrx * temp * pow(8., 0.5 * xd)) / 3.5e5, 0.4)) * 2.;
			double term1 = exp(pow(xReF / 3.5e5, 0.4)) * 2.;
//			temp = xReF / (jFrx * 3.5e4 + 2.5e4) * xd * xd;
			temp = xReF / 6.0e4 * xd * xd;
			if (temp < 85.) {
				term2 = pow(xReF, 0.035) * -1.7 * exp(-temp);
			}
			else {
				term2 = 0.;
			}
//			double term3 = (jFrx * .26 + ajFrx * .85) * pow(xd, 0.1) * pow(xReF, 0.001);
			double term3 = .26 * pow(xd, 0.1) * pow(xReF, 0.001);
			double c10 = term1 + term2 + term3;

			if (c10 > 0.) {
				double b2 = pow(1. / (xReF * .05 / 3.5e5 + 1.), 0.4);
				c3vert = pow(c10 / 2., b2) * 2.;
			}
// it only applies to ReF<0. and ReG>0., i.e. countercurrent
			//if (ReG > 0.) { /* THE LARGE DH MODEL IS BASED */
			//	if (dh > d1) { /* ON THE CSA SIMULATIONS OF */
			//		if (dh < 1.) { /* THE GE 1 AND 4 FT TESTS */
			//			c3vert = ((dh - d1) * .6 / (1. - d1) + c3vert *
			//				(1. - (dh - d1) / (1. - d1))) * jFrx + c3vert * ajFrx;
			//		}
			//		else if (dh < 3.) {
			//			c3vert = (.6 - (dh - 1.) * .27) * jFrx + c3vert * ajFrx;
			//		}
			//		else {
			//			c3vert = jFrx * .06 + c3vert * ajFrx;
			//		}
			//	}
			//}
		}
		c3vert *= fra;
	}
	/*     HORIZONTAL ------------------------------------------------ */
	if (fra != 1.) {
//		if (ReF * ReG >= 0.) { /*  COCURRENT */
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
//		}
//		else { /* COUNTERCURRENT */
//			temp = fabs(ReF / 600.);
//			if (temp < 85.) {
//				c3horz = pow(fabs(ReF / 4.5e5), 0.19) * 4.5 + exp(-temp) * .5;
//			}
//			else {
//				c3horz = pow(fabs(ReF / 4.5e5), 0.19) * 4.5;
//			}
//		}

		if (ReG < 0.) {
			c3horz *= (fra - 1.);
		}
		else {
			c3horz *= (1. - fra);
		}
	}

	return fmin(10., c3vert + c3horz);
} /* c3func_ */

int chexal_c0func_init() {

	/* Local variables */
	double re, b1, b1horz, b1vert, c1xhorz, c1xvert;


	/*     Subroutine C0FUNC_init calculates several constants for C0FUNC in the  */
	/*     Chexal-Lellouche drift flux correlation.  The parameters computed */
	/*     here are orientation dependent. */

	/*    Input as global */
	/*      ReF     = superficial liquid Reynolds number */
	/*      ReG     = superficial vapor Reynolds number */

	/*     Output as global */
	/*      C1X     = exponent used for (1-alpha) for C1 constant */
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
		b1vert = 1.;
		if (re < 0.) {
			b1vert = 1e-10;
		}
		if (fabs(re / 6e4) < 85.) {
			b1vert = 1. / (exp(-re / 6e4) + 1.);
		}
		if (b1vert > .8) {
			b1vert = .8;
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

double chexal_c0func(double al, double alm1) {

	/* Local variables */
	double l, lh, lv, c0a, c0, denom;


	/*     Subroutine C0FUNC calculates several constants in the Chexal- */
	/*     Lellouche drift flux correlation.  The parameters computed */
	/*     here are orientation dependent. */

	/*    Input Parameters */
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


int chexal_iterat(double& alpha, double& C0, double& vgj) {
	/* Local variables */
	int ier = 0;
	double a, b, fa, fb, fc, abj;
	/*     ITERAT USES THE REGULA-FALSI METHOD (Anderson-Bjoerk) for co-current flow */
	/*     TO SOLVE THE CHEXAL-LELLOUCHE CORRELATION FOR jF OR jG AND VOID */
	/*     FRACTION. */

	 /*      alpha = void fraction                     (OUTPUT) */
	 /*      C0   = distribution parameter            (OUTPUT) */
	 /*      VGJ  = drift velocity                    (OUTPUT) */

	//co-current flow only
	// establish range 
	// upward or horizontal flow
	// solution is between x and homogeneous void fraction
	a = x;
	fa = chexal_poly(a, C0, vgj);
	if (jF > 0.) {
		if (fa > 0.) { // solution void will be lower than x 
							// to allow for continuity alpha is set to x
							// this happens on low mass velocities Rainer Jordan 2022-02-15
			alpha = x;
			C0 = 0.;
			vgj = 0.;
			return 0;
		}
		b = VoidHomogeneous;
		fb = chexal_poly(b, C0, vgj);
	}
	else {
		//downward flow, solution between x and 1
		b = .9999999; // close to 1
		fb = chexal_poly(b, C0, vgj);
	}
	/*     RANGE ESTABLISHED: NOW USE REGULA FALSI TO FIND SOLUTION */
	for (int i = 1; i <= 100; ++i) {
		if (fabs(b - a) < 1e-10) {
			alpha = b;
			break;
		}
		alpha = a - fa * (b - a) / (fb - fa);

		fc = chexal_poly(alpha, C0, vgj);
		if (fabs(fc) < 1e-14) {
			break;
		}
		if (ISIGN(fc) * ISIGN(fb) < 0) {//ISIGN to avoid underflow errors
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
		b = alpha;
		fb = fc;
	}

	return ier;
} /* iterat_ */

double chexal_poly(double alpha, double& C0, double& vgj) {
	/* Local variables */
	double diff, alm1;

	/*     POLY CALCULATES THE DIFFERENCE BETWEEN THE jG CALCULATED FROM */
	/*     THE VOID FRACTION AND THE ACTUAL jG, THIS DIFFERENCE IS USED */
	/*     BY THE ITERAT ROUTINE.   POLY IS CALLED BY SUBROUTINE ITERAT. */

	/*      ALPHA   = void fraction                  (INPUT) */
	/*      VGJ     = drift velocity                         (OUTPUT) */
	/*      C0      = drift flux distribution parameter      (OUTPUT) */


//	if (alpha >= .9999999) alpha = .9999999; limits already set in calling function
//	if (alpha <= 1e-7) alpha = 1e-7;
	alm1 = 1. - alpha;
	vgj = vgj0 * pow(alm1, c1x);
	C0 = chexal_c0func(alpha, alm1);
	diff = alpha * (C0 * jSum + vgj) - jG; /*  CALCULATE DIFFERENCE AND RETURN */

	return diff;
} /* poly_ */


