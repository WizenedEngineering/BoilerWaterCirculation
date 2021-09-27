/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*      Calculation of properties for  water and steam            */
/*      according to 1997 IFC formulations for industrial use     */
/*      Parameters :                                              */
/*      t1  Temperature in K                                      */
/*      p1  Pressure in MPa(abs)                                  */
/*      medium that should be calculated:                         */
/*              1  STEAM                                          */
/*              2  WATER                                          */
/*              3  Fluid (super critical)                         */
/*      Return values :                                           */
/*      Enthalpy in kJ/kg                                         */
/*      cp, cv in kJ/kgK                                          */
/*      spec. Volume in m3/kg                                     */
/*      dyn. Viscosity in Pas                                     */
/*      Conductivity in W/mK                                      */
/*      Surface tension in mN/m                                   */
/*      speed of sound (SoS) in m/s                               */
/*      temperature, saturation temperature in K                  */
/*      saturation pressure in MPa                                */
/*  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//#include "stdafx.h"
#include "propH2O.h"
#include "CommonHeader.h"

/* Table of constant values */
static const double GasConstH2O = 0.461526;
static double ngamma1[34] = { .14632971213167,-.84548187169114,
	-3.756360367204,3.3855169168385,-.95791963387872,.15772038513228,
	-.016616417199501,8.1214629983568e-4,2.8319080123804e-4,
	-6.0706301565874e-4,-.018990068218419,-.032529748770505,
	-.021841717175414,-5.283835796993e-5,-4.7184321073267e-4,
	-3.0001780793026e-4,4.7661393906987e-5,-4.4141845330846e-6,
	-7.2694996297594e-16,-3.1679644845054e-5,-2.8270797985312e-6,
	-8.5205128120103e-10,-2.2425281908e-6,-6.5171222895601e-7,
	-1.4341729937924e-13,-4.0516996860117e-7,-1.2734301741641e-9,
	-1.7424871230634e-10,-6.8762131295531e-19,1.4478307828521e-20,
	2.6335781662795e-23,-1.1947622640071e-23,1.8228094581404e-24,
	-9.3537087292458e-26 };

static double irgamma2[43] = { 1.,1.,1.,1.,1.,2.,2.,2.,2.,2.,3.,3.,3.,
	3.,3.,4.,4.,4.,5.,6.,6.,6.,7.,7.,7.,8.,8.,9.,10.,10.,10.,16.,16.,
	18.,20.,20.,20.,21.,22.,23.,24.,24.,24. };

static double jrgamma2[43] = { 0.,1.,2.,3.,6.,1.,2.,4.,7.,36.,0.,1.,
	3.,6.,35.,1.,2.,3.,7.,3.,16.,35.,0.,11.,25.,8.,36.,13.,4.,10.,14.,
	29.,50.,57.,20.,35.,48.,21.,53.,39.,26.,40.,58. };

static double nrgamma2[43] = { -.0017731742473213,-.017834862292358,
	-.045996013696365,-.057581259083432,-.05032527872793,
	-3.3032641670203e-5,-1.8948987516315e-4,-.0039392777243355,
	-.043797295650573,-2.6674547914087e-5,2.0481737692309e-8,
	4.3870667284435e-7,-3.227767723857e-5,-.0015033924542148,
	-.040668253562649,-7.8847309559367e-10,1.2790717852285e-8,
	4.8225372718507e-7,2.2922076337661e-6,-1.6714766451061e-11,
	-.0021171472321355,-23.895741934104,-5.905956432427e-18,
	-1.2621808899101e-6,-.038946842435739,1.1256211360459e-11,
	-8.2311340897998,1.9809712802088e-8,1.0406965210174e-19,
	-1.0234747095929e-13,-1.0018179379511e-9,-8.0882908646985e-11,
	.10693031879409,-.33662250574171,8.9185845355421e-25,
	3.0629316876232e-13,-4.2002467698208e-6,-5.9056029685639e-26,
	3.7826947613457e-6,-1.2768608934681e-15,7.3087610595061e-29,
	5.5414715350778e-17,-9.436970724121e-7 };

static double n0gamma2[9] = { -9.6927686500217,10.086655968018,
	-.005608791128302,.071452738081455,-.40710498223928,
	1.4240819171444,-4.383951131945,-.28408632460772,.021268463753307
};

static double n3[40] = { 1.0658070028513,-15.732845290239,
	20.944396974307,-7.6867707878716,2.6185947787954,-2.808078114862,
	1.2053369696517,-.0084566812812502,-1.2654315477714,
	-1.1524407806681,.88521043984318,-.64207765181607,.38493460186671,
	-.85214708824206,4.8972281541877,-3.0502617256965,
	.039420536879154,.12558408424308,-.2799932969871,1.389979956946,
	-2.018991502357,-.0082147637173963,-.47596035734923,
	.0439840744735,-.44476435428739,.90572070719733,.70522450087967,
	.10770512626332,-.32913623258954,-.50871062041158,
	-.022175400873096,.094260751665092,.16436278447961,
	-.013503372241348,-.014834345352472,5.7922953628084e-4,
	.0032308904703711,8.0964802996215e-5,-1.6557679795037e-4,
	-4.4923899061815e-5 };

/* gamma function in region 1(WATER T<350 degC)*/
double H2O::gamma1(double pNorm, double tNorm)
{
	double ret_val;
	static double pFact, tFact, pFact_2, pFact_4, tFact_2,
		tFact_3, tFact_4, tFact_5, tFact_6, tFact_7, tFact_8,
		tFact_9, pFact_8, tFact_10, tFact_11, tFact_31,
		tFact_40, tFact_41, tFact_17, tFact_29, tFact_38,
		tFact_39;

	pFact = 7.1 - pNorm;
	tFact = tNorm - 1.222;
	tFact_2 = tFact * tFact;
	tFact_3 = tFact * tFact_2;
	tFact_4 = tFact_2 * tFact_2;
	tFact_5 = tFact_4 * tFact;
	tFact_6 = tFact_5 * tFact;
	tFact_7 = tFact_6 * tFact;
	tFact_8 = tFact_7 * tFact;
	tFact_9 = tFact_8 * tFact;
	tFact_10 = tFact_9 * tFact;
	tFact_11 = tFact_10 * tFact;
	tFact_17 = tFact_11 * tFact_6;
	tFact_29 = tFact_17 * tFact_11 * tFact;
	tFact_31 = tFact_29 * tFact_2;
	tFact_38 = tFact_31 * tFact_7;
	tFact_39 = tFact_38 * tFact;
	tFact_40 = tFact_39 * tFact;
	tFact_41 = tFact_40 * tFact;
	pFact_2 = pFact * pFact;
	pFact_4 = pFact_2 * pFact_2;
	pFact_8 = pFact_4 * pFact_4;
	ret_val = ngamma1[0] / tFact_2 + ngamma1[1] / tFact + ngamma1[2] +
		ngamma1[3] * tFact + ngamma1[4] * tFact_2 + ngamma1[5] *
		tFact_3 + ngamma1[6] * tFact_4 + ngamma1[7] * tFact_5 +
		pFact * (ngamma1[8] / tFact_9 + ngamma1[9] / tFact_7 +
			ngamma1[10] / tFact + ngamma1[11] + ngamma1[12] * tFact + ngamma1[
				13] * tFact_3 + pFact * (ngamma1[14] / tFact_3 + ngamma1[15]
					+ ngamma1[16] * tFact + ngamma1[17] * tFact_3 + ngamma1[18] *
					tFact_17 + pFact * (ngamma1[19] / tFact_4 + ngamma1[20] +
						ngamma1[21] * tFact_6 + pFact * (ngamma1[22] / tFact_5 +
							ngamma1[23] / tFact_2 + ngamma1[24] * tFact_10 + pFact * (
								ngamma1[25] / tFact_8 + pFact_2 * pFact * (ngamma1[26] /
									tFact_11 + ngamma1[27] / tFact_6 + pFact_8 * pFact_4 *
									pFact * (ngamma1[28] / tFact_29 + pFact_2 * (ngamma1[29] /
										tFact_31 + pFact_4 * pFact_2 * (ngamma1[30] / tFact_38 +
											pFact * (ngamma1[31] / tFact_39 + pFact * (ngamma1[32] /
												tFact_40 + pFact * ngamma1[33] / tFact_41)))))))))));
	return ret_val;
} /* gamma1_ */

/* ------------------------------------------------------------------------------- */
/* derivative of gamma with respect to tau region 1 */
double H2O::gamma1Tau(double pNorm, double tNorm)
{
	double ret_val;
	static double pFact, tFact, pFact_2, pFact_4, tFact_2,
		tFact_3, tFact_4, tFact_5, tFact_6, tFact_7, tFact_8,
		tFact_9, pFact_8, tFact_10, tFact_12, tFact_30,
		tFact_32, tFact_40, tFact_16, tFact_41, tFact_42,
		tFact_39;

	pFact = 7.1 - pNorm;
	tFact = tNorm - 1.222;
	tFact_2 = tFact * tFact;
	tFact_3 = tFact * tFact_2;
	tFact_4 = tFact_2 * tFact_2;
	tFact_5 = tFact_4 * tFact;
	tFact_6 = tFact_5 * tFact;
	tFact_7 = tFact_6 * tFact;
	tFact_8 = tFact_7 * tFact;
	tFact_9 = tFact_8 * tFact;
	tFact_10 = tFact_9 * tFact;
	tFact_12 = tFact_10 * tFact_2;
	tFact_16 = tFact_12 * tFact_4;
	tFact_30 = tFact_16 * tFact_12 * tFact_2;
	tFact_32 = tFact_30 * tFact_2;
	tFact_39 = tFact_32 * tFact_7;
	tFact_40 = tFact_39 * tFact;
	tFact_41 = tFact_40 * tFact;
	tFact_42 = tFact_41 * tFact;
	pFact_2 = pFact * pFact;
	pFact_4 = pFact_2 * pFact_2;
	pFact_8 = pFact_4 * pFact_4;

	ret_val = ngamma1[0] * -2. / tFact_3 - ngamma1[1] / tFact_2 + ngamma1[3]
		 + ngamma1[4] * 2. * tFact + ngamma1[5] * 3. * tFact_2 +
			ngamma1[6] * 4. * tFact_3 + ngamma1[7] * 5. * tFact_4 + pFact
			* (ngamma1[8] * -9. / tFact_10 - ngamma1[9] * 7. / tFact_8 -
				ngamma1[10] / tFact_2 + ngamma1[12] + ngamma1[13] * 3. *
				tFact_2 + pFact * (ngamma1[14] * -3. / tFact_4 + ngamma1[16]
					+ ngamma1[17] * 3. * tFact_2 + ngamma1[18] * 17. * tFact_16 +
					pFact * (ngamma1[19] * -4. / tFact_5 + ngamma1[21] * 6. *
						tFact_5 + pFact * (ngamma1[22] * -5. / tFact_6 - ngamma1[23] *
							2. / tFact_3 + ngamma1[24] * 10. * tFact_9 + pFact * (
								ngamma1[25] * -8. / tFact_9 + pFact_2 * pFact * (ngamma1[26] *
									-11. / tFact_12 - ngamma1[27] * 6. / tFact_7 + pFact_8 *
									pFact_4 * pFact * (ngamma1[28] * -29. / tFact_30 + pFact_2 *
										(ngamma1[29] * -31. / tFact_32 + pFact_4 * pFact_2 * 
											(ngamma1[30] * -38. / tFact_39 + pFact * (ngamma1[31] * -39. /
												tFact_40 + pFact * (ngamma1[32] * -40. / tFact_41 + pFact *
													-41. * ngamma1[33] / tFact_42)))))))))));

		return ret_val;
} /* gamma1Tau_ */

/* ----------------------------------------------------------------------- */
/* derivative of gamma with respect to pi region 1 */
double H2O::gamma1Pi(double pNorm, double tNorm)
{
	double ret_val;
	static double pFact, tFact, pFact_2, pFact_4, tFact_2,
		tFact_3, tFact_4, tFact_5, tFact_6, tFact_7, tFact_8,
		tFact_9, pFact_8, tFact_10, tFact_11, tFact_31,
		tFact_40, tFact_41, tFact_17, tFact_29, tFact_38,
		tFact_39;

	pFact = 7.1 - pNorm;
	tFact = tNorm - 1.222;
	tFact_2 = tFact * tFact;
	tFact_3 = tFact * tFact_2;
	tFact_4 = tFact_2 * tFact_2;
	tFact_5 = tFact_4 * tFact;
	tFact_6 = tFact_5 * tFact;
	tFact_7 = tFact_6 * tFact;
	tFact_8 = tFact_7 * tFact;
	tFact_9 = tFact_8 * tFact;
	tFact_10 = tFact_9 * tFact;
	tFact_11 = tFact_10 * tFact;
	tFact_17 = tFact_11 * tFact_6;
	tFact_29 = tFact_17 * tFact_11 * tFact;
	tFact_31 = tFact_29 * tFact_2;
	tFact_38 = tFact_31 * tFact_7;
	tFact_39 = tFact_38 * tFact;
	tFact_40 = tFact_39 * tFact;
	tFact_41 = tFact_40 * tFact;
	pFact_2 = pFact * pFact;
	pFact_4 = pFact_2 * pFact_2;
	pFact_8 = pFact_4 * pFact_4;
	ret_val = -ngamma1[8] / tFact_9 - ngamma1[9] / tFact_7 - ngamma1[10] /
		tFact - ngamma1[11] - ngamma1[12] * tFact - ngamma1[13] *
		tFact_3 + pFact * ((-ngamma1[14] / tFact_3 - ngamma1[15] -
			ngamma1[16] * tFact - ngamma1[17] * tFact_3 - ngamma1[18] *
			tFact_17) * 2. + pFact * ((-ngamma1[19] / tFact_4 - ngamma1[20]
				 - ngamma1[21] * tFact_6) * 3. + pFact * ((-ngamma1[22] /
					tFact_5 - ngamma1[23] / tFact_2 - ngamma1[24] * tFact_10) *
					4. + pFact * (-ngamma1[25] * 5. / tFact_8 + pFact_2 * pFact *
						((-ngamma1[26] / tFact_11 - ngamma1[27] / tFact_6) * 8. +
							pFact_8 * pFact_4 * pFact * (ngamma1[28] * -21. / tFact_29
								+ pFact_2 * (ngamma1[29] * -23. / tFact_31 + pFact_4 *
									pFact_2 * (ngamma1[30] * -29. / tFact_38 + pFact * (ngamma1[31]
										 * -30. / tFact_39 + pFact * (ngamma1[32] * -31. /
											tFact_40 - pFact * 32. * ngamma1[33] / tFact_41))))))))));
	return ret_val;
} /* gamma1pi_ */

/* ---------------------------------------------------------------------- */
/*second derivative of gamma with respect to pi region 1 */
double H2O::gamma1PiPi(double pNorm, double tNorm)
{
	double ret_val;
	static double pFact, tFact, pFact_2, pFact_4, tFact_2,
		tFact_3, tFact_4, tFact_5, tFact_6, tFact_7, tFact_8,
		tFact_9, pFact_8, tFact_10, tFact_11, tFact_31,
		tFact_40, tFact_41, tFact_17, tFact_29, tFact_38,
		tFact_39;

	pFact = 7.1 - pNorm;
	tFact = tNorm - 1.222;
	tFact_2 = tFact * tFact;
	tFact_3 = tFact * tFact_2;
	tFact_4 = tFact_2 * tFact_2;
	tFact_5 = tFact_4 * tFact;
	tFact_6 = tFact_5 * tFact;
	tFact_7 = tFact_6 * tFact;
	tFact_8 = tFact_7 * tFact;
	tFact_9 = tFact_8 * tFact;
	tFact_10 = tFact_9 * tFact;
	tFact_11 = tFact_10 * tFact;
	tFact_17 = tFact_11 * tFact_6;
	tFact_29 = tFact_17 * tFact_11 * tFact;
	tFact_31 = tFact_29 * tFact_2;
	tFact_38 = tFact_31 * tFact_7;
	tFact_39 = tFact_38 * tFact;
	tFact_40 = tFact_39 * tFact;
	tFact_41 = tFact_40 * tFact;
	pFact_2 = pFact * pFact;
	pFact_4 = pFact_2 * pFact_2;
	pFact_8 = pFact_4 * pFact_4;
	ret_val = (ngamma1[14] / tFact_3 + ngamma1[15] + ngamma1[16] * tFact +
		ngamma1[17] * tFact_3 + ngamma1[18] * tFact_17) * 2. + pFact *
		((ngamma1[19] / tFact_4 + ngamma1[20] + ngamma1[21] *
			tFact_6) * 6. + pFact * ((ngamma1[22] / tFact_5 + ngamma1[23]
				/ tFact_2 + ngamma1[24] * tFact_10) * 12. + pFact * (ngamma1[25]
					 * 20. / tFact_8 + pFact_2 * pFact * ((ngamma1[26] /
						tFact_11 + ngamma1[27] / tFact_6) * 56. + pFact_8 *
						pFact_4 * pFact * (ngamma1[28] * 420. / tFact_29 + pFact_2 *
							(ngamma1[29] * 506. / tFact_31 + pFact_4 * pFact_2 * 
								(ngamma1[30] * 812. / tFact_38 + pFact * (ngamma1[31] * 870. /
									tFact_39 + pFact * (ngamma1[32] * 930. / tFact_40 + pFact *
										992. * ngamma1[33] / tFact_41)))))))));
	return ret_val;
} /* gamma1pipi_ */

/* --------------------------------------------------------------------------- */
/*second  derivative of gamma with respect to  tau region 1 */
double H2O::gamma1TauTau(double pNorm, double tNorm)
{
	double ret_val;
	static double pFact, tFact, pFact_2, pFact_4, tFact_2,
		tFact_3, tFact_4, tFact_5, tFact_6, tFact_7, tFact_8,
		tFact_9, pFact_8, tFact_10, tFact_11, tFact_13,
		tFact_31, tFact_15, tFact_33, tFact_40, tFact_41,
		tFact_42, tFact_43;

	pFact = 7.1 - pNorm;
	tFact = tNorm - 1.222;
	tFact_2 = tFact * tFact;
	tFact_3 = tFact * tFact_2;
	tFact_4 = tFact_2 * tFact_2;
	tFact_5 = tFact_4 * tFact;
	tFact_6 = tFact_5 * tFact;
	tFact_7 = tFact_6 * tFact;
	tFact_8 = tFact_7 * tFact;
	tFact_9 = tFact_8 * tFact;
	tFact_10 = tFact_9 * tFact;
	tFact_11 = tFact_10 * tFact;
	tFact_13 = tFact_11 * tFact_2;
	tFact_15 = tFact_13 * tFact_2;
	tFact_31 = tFact_15 * tFact_15 * tFact;
	tFact_33 = tFact_31 * tFact_2;
	tFact_40 = tFact_33 * tFact_7;
	tFact_41 = tFact_40 * tFact;
	tFact_42 = tFact_41 * tFact;
	tFact_43 = tFact_42 * tFact;
	pFact_2 = pFact * pFact;
	pFact_4 = pFact_2 * pFact_2;
	pFact_8 = pFact_4 * pFact_4;
	ret_val = ngamma1[0] * 6. / tFact_4 + ngamma1[1] * 2. / tFact_3 +
		ngamma1[4] * 2. + ngamma1[5] * 6. * tFact + ngamma1[6] * 12. *
		tFact_2 + ngamma1[7] * 20. * tFact_3 + pFact * (ngamma1[8] *
			90. / tFact_11 + ngamma1[9] * 56. / tFact_9 + ngamma1[10] *
			2. / tFact_3 + ngamma1[13] * 6. * tFact + pFact * (ngamma1[14] *
				12. / tFact_5 + ngamma1[17] * 6. * tFact + ngamma1[18] * 272. *
				tFact_15 + pFact * (ngamma1[19] * 20. / tFact_6 + ngamma1[21]
					* 30. * tFact_4 + pFact * (ngamma1[22] * 30. / tFact_7 +
						ngamma1[23] * 6. / tFact_4 + ngamma1[24] * 90. * tFact_8 +
						pFact * (ngamma1[25] * 72. / tFact_10 + pFact_2 * pFact * (
							ngamma1[26] * 132. / tFact_13 + ngamma1[27] * 42. / tFact_8 +
							pFact_8 * pFact_4 * pFact * (ngamma1[28] * 870. / tFact_31
								+ pFact_2 * (ngamma1[29] * 992. / tFact_33 + pFact_4 *
									pFact_2 * (ngamma1[30] * 1482. / tFact_40 + pFact * (ngamma1[31]
										 * 1560. / tFact_41 + pFact * (ngamma1[32] * 1640. /
											tFact_42 + pFact * 1722. * ngamma1[33] / tFact_43)))))))))));
	return ret_val;
} /* gamma1TauTau_ */

/* ------------------------------------------------------------------------- */
/*derivative of gamma with respect to pi and tau region 1 */
double H2O::gamma1PiTau(double pNorm, double tNorm)
{
	double ret_val;
	static double pFact, tFact, pFact_2, pFact_4, tFact_2,
		tFact_3, tFact_4, tFact_5, tFact_6, tFact_7, tFact_8,
		tFact_9, pFact_8, tFact_10, tFact_12, tFact_30,
		tFact_32, tFact_40, tFact_16, tFact_41, tFact_42,
		tFact_39;

	pFact = 7.1 - pNorm;
	tFact = tNorm - 1.222;
	tFact_2 = tFact * tFact;
	tFact_3 = tFact * tFact_2;
	tFact_4 = tFact_2 * tFact_2;
	tFact_5 = tFact_4 * tFact;
	tFact_6 = tFact_5 * tFact;
	tFact_7 = tFact_6 * tFact;
	tFact_8 = tFact_7 * tFact;
	tFact_9 = tFact_8 * tFact;
	tFact_10 = tFact_9 * tFact;
	tFact_12 = tFact_10 * tFact_2;
	tFact_16 = tFact_10 * tFact_6;
	tFact_30 = tFact_10 * tFact_10 * tFact_10;
	tFact_32 = tFact_30 * tFact_2;
	tFact_39 = tFact_32 * tFact_7;
	tFact_40 = tFact_39 * tFact;
	tFact_41 = tFact_40 * tFact;
	tFact_42 = tFact_41 * tFact;
	pFact_2 = pFact * pFact;
	pFact_4 = pFact_2 * pFact_2;
	pFact_8 = pFact_4 * pFact_4;
	ret_val = ngamma1[8] * 9. / tFact_10 + ngamma1[9] * 7. / tFact_8 +
		ngamma1[10] / tFact_2 - ngamma1[12] - ngamma1[13] * 3. *
		tFact_2 + pFact * (ngamma1[14] * 6. / tFact_4 - ngamma1[16] *
			2. - ngamma1[17] * 6. * tFact_2 - ngamma1[18] * 34. *
			tFact_16 + pFact * (ngamma1[19] * 12. / tFact_5 - ngamma1[21]
				* 18. * tFact_5 + pFact * (ngamma1[22] * 20. / tFact_6 +
					ngamma1[23] * 8. / tFact_3 - ngamma1[24] * 40. * tFact_9 +
					pFact * (ngamma1[25] * 40. / tFact_9 + pFact_2 * pFact * (
						ngamma1[26] * 88. / tFact_12 + ngamma1[27] * 48. / tFact_7 +
						pFact_8 * pFact_4 * pFact * (ngamma1[28] * 609. / tFact_30
							+ pFact_2 * (ngamma1[29] * 713. / tFact_32 + pFact_4 *
								pFact_2 * (ngamma1[30] * 1131. / tFact_39 + pFact * (ngamma1[31]
									 * 1170. / tFact_40 + pFact * (ngamma1[32] * 1240. /
										tFact_41 + pFact * 1312. * ngamma1[33] / tFact_42))))))))));
	return ret_val;
} /* gamma1pitau_ */

/* ---------------------------------------------------------------------------- */
 /* gamma function in region 2(STEAM)*/
double H2O::gamma2r(double pNorm, double tNorm)
{
	double ret_val;
	int i;
	double tFact;

	tFact = tNorm - .5;
	ret_val = 0.;
	for (i = 0; i <= 42; ++i) {
		ret_val += nrgamma2[i] * pow(pNorm, irgamma2[i]) * pow(tFact, jrgamma2[i]);
	}
	return ret_val;
} /* gamma2r */

/* ----------------------------------------------------------------------- */
/* gamma function in region 2(STEAM)*/
double H2O::gamma2(double pNorm, double tNorm)
{
	/* Initialized data */
	static double j0gamma2[9] = { 0.,1.,-5.,-4.,-3.,-2.,-1.,2.,3. };
	//static double n0gamma2[9] = { -9.6927686500217,10.086655968018,
	  // -.005608791128302,.071452738081455,-.40710498223928,
	  // 1.4240819171444,-4.383951131945,-.28408632460772,.021268463753307};

	double ret_val;
	int i;
	double gamma0;

	gamma0 = log(pNorm) + n0gamma2[0] + n0gamma2[1] * tNorm;
	for (i = 2; i <= 8; ++i) {
		gamma0 += n0gamma2[i] * pow(tNorm, j0gamma2[i]);
	}
	ret_val = gamma0 + gamma2r(pNorm, tNorm);
	return ret_val;
} /* gamma2 */

/* ---------------------------------------------------------------------- */
/* derivative of gamma with respect to pi region 2 */
double H2O::gamma2rPi(double pNorm, double tNorm)
{
	double ret_val;
	static double tFact, pFact_2, pFact_4, tFact_2, tFact_3,
		tFact_4, tFact_6, tFact_7, tFact_8, tFact_10,
		tFact_11, tFact_20, tFact_13, tFact_14, tFact_21,
		tFact_16, tFact_25, tFact_26, tFact_35, tFact_29,
		tFact_36, tFact_39, tFact_40, tFact_48, tFact_50,
		tFact_53, tFact_57, tFact_58;

	tFact = tNorm - .5;
	tFact_2 = tFact * tFact;
	tFact_3 = tFact * tFact_2;
	tFact_4 = tFact_2 * tFact_2;
	tFact_6 = tFact_4 * tFact_2;
	tFact_7 = tFact_6 * tFact;
	tFact_8 = tFact_7 * tFact;
	tFact_10 = tFact_8 * tFact_2;
	tFact_11 = tFact_10 * tFact;
	tFact_13 = tFact_11 * tFact_2;
	tFact_14 = tFact_13 * tFact;
	tFact_16 = tFact_10 * tFact_6;
	tFact_20 = tFact_10 * tFact_10;
	tFact_21 = tFact_20 * tFact;
	tFact_25 = tFact_21 * tFact_4;
	tFact_26 = tFact_25 * tFact;
	tFact_29 = tFact_26 * tFact_3;
	tFact_35 = tFact_29 * tFact_6;
	tFact_36 = tFact_35 * tFact;
	tFact_39 = tFact_36 * tFact_3;
	tFact_40 = tFact_39 * tFact;
	tFact_48 = tFact_40 * tFact_8;
	tFact_50 = tFact_48 * tFact_2;
	tFact_53 = tFact_50 * tFact_3;
	tFact_57 = tFact_53 * tFact_4;
	tFact_58 = tFact_57 * tFact;
	pFact_2 = pNorm * pNorm;
	pFact_4 = pFact_2 * pFact_2;
	ret_val = nrgamma2[0] + nrgamma2[1] * tFact + nrgamma2[2] * tFact_2 +
		nrgamma2[3] * tFact_3 + nrgamma2[4] * tFact_6 + pNorm * ((
			nrgamma2[5] * tFact + nrgamma2[6] * tFact_2 + nrgamma2[7] *
			tFact_4 + nrgamma2[8] * tFact_7 + nrgamma2[9] * tFact_36) *
			2. + pNorm * ((nrgamma2[10] + nrgamma2[11] * tFact + nrgamma2[12]
				* tFact_3 + nrgamma2[13] * tFact_6 + nrgamma2[14] *
				tFact_35) * 3. + pNorm * ((nrgamma2[15] * tFact + nrgamma2[16]
					* tFact_2 + nrgamma2[17] * tFact_3) * 4. + pNorm * (nrgamma2[
						18] * 5. * tFact_7 + pNorm * ((nrgamma2[19] * tFact_3 +
							nrgamma2[20] * tFact_16 + nrgamma2[21] * tFact_35) * 6. +
							pNorm * ((nrgamma2[22] + nrgamma2[23] * tFact_11 + nrgamma2[24]
								* tFact_25) * 7. + pNorm * ((nrgamma2[25] * tFact_8 +
									nrgamma2[26] * tFact_36) * 8. + pNorm * (nrgamma2[27] * 9. *
										tFact_13 + pNorm * ((nrgamma2[28] * tFact_4 + nrgamma2[29] *
											tFact_10 + nrgamma2[30] * tFact_14) * 10. + pFact_2 *
											pFact_4 * ((nrgamma2[31] * tFact_29 + nrgamma2[32] *
												tFact_50) * 16. + pFact_2 * (nrgamma2[33] * 18. * tFact_57
													+ pFact_2 * ((nrgamma2[34] * tFact_20 + nrgamma2[35] *
														tFact_35 + nrgamma2[36] * tFact_48) * 20. + pNorm * (
															nrgamma2[37] * 21. * tFact_21 + pNorm * (nrgamma2[38] * 22. *
																tFact_53 + pNorm * (nrgamma2[39] * 23. * tFact_39 + pNorm *
																	24. * (nrgamma2[40] * tFact_26 + nrgamma2[41] * tFact_40 +
																		nrgamma2[42] * tFact_58))))))))))))))));
	return ret_val;
} /* gamma2rpi_ */

/* ------------------------------------------------------------------------- */
/* derivative of gamma with respect to pi region 2 */
double H2O::gamma2Pi(double pNorm, double tNorm)
{
	double ret_val;
	double gamma0;

	gamma0 = 1. / pNorm;
	ret_val = gamma0 + gamma2rPi(pNorm, tNorm);
	return ret_val;
} /* gamma2pi_ */

/* --------------------------------------------------------------------------- */
/* derivative of gamma with respect to pi and tau region 2 */
double H2O::gamma2rPiTau(double pNorm, double tNorm)
{
	double ret_val;
	static double tFact, pFact_2, pFact_4, tFact_2, tFact_3,
		tFact_5, tFact_6, tFact_7, tFact_9, tFact_10,
		tFact_20, tFact_12, tFact_13, tFact_15, tFact_24,
		tFact_25, tFact_34, tFact_19, tFact_28, tFact_35,
		tFact_38, tFact_39, tFact_47, tFact_49, tFact_52,
		tFact_56, tFact_57;

	tFact = tNorm - .5;
	tFact_2 = tFact * tFact;
	tFact_3 = tFact * tFact_2;
	tFact_5 = tFact_2 * tFact_3;
	tFact_6 = tFact_5 * tFact;
	tFact_7 = tFact_6 * tFact;
	tFact_9 = tFact_7 * tFact_2;
	tFact_10 = tFact_9 * tFact;
	tFact_12 = tFact_10 * tFact_2;
	tFact_13 = tFact_12 * tFact;
	tFact_15 = tFact_13 * tFact_2;
	tFact_19 = tFact_13 * tFact_6;
	tFact_20 = tFact_10 * tFact_10;
	tFact_24 = tFact_19 * tFact_5;
	tFact_25 = tFact_20 * tFact_5;
	tFact_28 = tFact_25 * tFact_3;
	tFact_34 = tFact_28 * tFact_6;
	tFact_35 = tFact_34 * tFact;
	tFact_38 = tFact_35 * tFact_3;
	tFact_39 = tFact_38 * tFact;
	tFact_47 = tFact_38 * tFact_9;
	tFact_49 = tFact_47 * tFact_2;
	tFact_52 = tFact_49 * tFact_3;
	tFact_56 = tFact_49 * tFact_7;
	tFact_57 = tFact_56 * tFact;
	pFact_2 = pNorm * pNorm;
	pFact_4 = pFact_2 * pFact_2;
	ret_val = nrgamma2[1] + nrgamma2[2] * 2. * tFact + nrgamma2[3] * 3. *
		tFact_2 + nrgamma2[4] * 6. * tFact_5 + pNorm * (nrgamma2[5] *
			2. + nrgamma2[6] * 4. * tFact + nrgamma2[7] * 8. * tFact_3 +
			nrgamma2[8] * 14. * tFact_6 + nrgamma2[9] * 72. * tFact_35 +
			pNorm * (nrgamma2[11] * 3. + nrgamma2[12] * 9. * tFact_2 +
				nrgamma2[13] * 18. * tFact_5 + nrgamma2[14] * 105. * tFact_34
				+ pNorm * (nrgamma2[15] * 4. + nrgamma2[16] * 8. * tFact +
					nrgamma2[17] * 12. * tFact_2 + pNorm * (nrgamma2[18] * 35. *
						tFact_6 + pNorm * (nrgamma2[19] * 18. * tFact_2 + nrgamma2[
							20] * 96. * tFact_15 + nrgamma2[21] * 210. * tFact_34 +
								pNorm * (nrgamma2[23] * 77. * tFact_10 + nrgamma2[24] * 175. *
									tFact_24 + pNorm * (nrgamma2[25] * 64. * tFact_7 + nrgamma2[
										26] * 288. * tFact_35 + pNorm * (nrgamma2[27] * 117. *
											tFact_12 + pNorm * (nrgamma2[28] * 40. * tFact_3 + nrgamma2[
												29] * 100. * tFact_9 + nrgamma2[30] * 140. * tFact_13 +
													pFact_2 * pFact_4 * (nrgamma2[31] * 464. * tFact_28 +
														nrgamma2[32] * 800. * tFact_49 + pFact_2 * (nrgamma2[33] *
															1026. * tFact_56 + pFact_2 * (nrgamma2[34] * 400. *
																tFact_19 + nrgamma2[35] * 700. * tFact_34 + nrgamma2[36] *
																960. * tFact_47 + pNorm * (nrgamma2[37] * 441. * tFact_20 +
																	pNorm * (nrgamma2[38] * 1166. * tFact_52 + pNorm * (nrgamma2[
																		39] * 897. * tFact_38 + pNorm * (nrgamma2[40] * 624. *
																			tFact_25 + nrgamma2[41] * 960. * tFact_39 + nrgamma2[42] *
																			1392. * tFact_57))))))))))))))));
	return ret_val;
} /* gamma2rpitau_ */

/* ----------------------------------------------------------------------- */
/* derivative of gamma with respect to pi and tau region 2 */
double H2O::gamma2PiTau(double pNorm, double tNorm)
{
	double ret_val;

	ret_val = gamma2rPiTau(pNorm, tNorm);
	return ret_val;
} /* gamma2pitau_ */

/* -------------------------------------------------------------------------- */
/*second derivative of gamma with respect to pi region 2 */
double H2O::gamma2rPiPi(double pNorm, double tNorm)
{
	double ret_val;
	static double tFact, pFact_2, pFact_4, tFact_2, tFact_3,
		tFact_4, tFact_6, tFact_7, tFact_8, tFact_10,
		tFact_11, tFact_20, tFact_13, tFact_14, tFact_21,
		tFact_16, tFact_25, tFact_26, tFact_35, tFact_29,
		tFact_36, tFact_39, tFact_40, tFact_48, tFact_50,
		tFact_53, tFact_57, tFact_58;

	tFact = tNorm - .5;
	tFact_2 = tFact * tFact;
	tFact_3 = tFact * tFact_2;
	tFact_4 = tFact_2 * tFact_2;
	tFact_6 = tFact_4 * tFact_2;
	tFact_7 = tFact_6 * tFact;
	tFact_8 = tFact_7 * tFact;
	tFact_10 = tFact_8 * tFact_2;
	tFact_11 = tFact_10 * tFact;
	tFact_13 = tFact_11 * tFact_2;
	tFact_14 = tFact_13 * tFact;
	tFact_16 = tFact_10 * tFact_6;
	tFact_20 = tFact_10 * tFact_10;
	tFact_21 = tFact_20 * tFact;
	tFact_25 = tFact_21 * tFact_4;
	tFact_26 = tFact_25 * tFact;
	tFact_29 = tFact_26 * tFact_3;
	tFact_35 = tFact_29 * tFact_6;
	tFact_36 = tFact_35 * tFact;
	tFact_39 = tFact_36 * tFact_3;
	tFact_40 = tFact_39 * tFact;
	tFact_48 = tFact_40 * tFact_8;
	tFact_50 = tFact_48 * tFact_2;
	tFact_53 = tFact_50 * tFact_3;
	tFact_57 = tFact_53 * tFact_4;
	tFact_58 = tFact_57 * tFact;
	pFact_2 = pNorm * pNorm;
	pFact_4 = pFact_2 * pFact_2;
	ret_val = (nrgamma2[5] * tFact + nrgamma2[6] * tFact_2 + nrgamma2[7] *
		tFact_4 + nrgamma2[8] * tFact_7 + nrgamma2[9] * tFact_36) *
		2. + pNorm * ((nrgamma2[10] + nrgamma2[11] * tFact + nrgamma2[12]
			* tFact_3 + nrgamma2[13] * tFact_6 + nrgamma2[14] *
			tFact_35) * 6. + pNorm * ((nrgamma2[15] * tFact + nrgamma2[16]
				* tFact_2 + nrgamma2[17] * tFact_3) * 12. + pNorm * (
					nrgamma2[18] * 20. * tFact_7 + pNorm * ((nrgamma2[19] *
						tFact_3 + nrgamma2[20] * tFact_16 + nrgamma2[21] * tFact_35)
						* 30. + pNorm * ((nrgamma2[22] + nrgamma2[23] * tFact_11 +
							nrgamma2[24] * tFact_25) * 42. + pNorm * ((nrgamma2[25] *
								tFact_8 + nrgamma2[26] * tFact_36) * 56. + pNorm * (nrgamma2[
									27] * 72. * tFact_13 + pNorm * ((nrgamma2[28] * tFact_4 +
										nrgamma2[29] * tFact_10 + nrgamma2[30] * tFact_14) * 90. +
										pFact_2 * pFact_4 * ((nrgamma2[31] * tFact_29 + nrgamma2[32]
											* tFact_50) * 240. + pFact_2 * (nrgamma2[33] * 306. *
												tFact_57 + pFact_2 * ((nrgamma2[34] * tFact_20 + nrgamma2[35]
													* tFact_35 + nrgamma2[36] * tFact_48) * 380. + pNorm * (
														nrgamma2[37] * 420. * tFact_21 + pNorm * (nrgamma2[38] * 462. *
															tFact_53 + pNorm * (nrgamma2[39] * 506. * tFact_39 + pNorm
																* 552. * (nrgamma2[40] * tFact_26 + nrgamma2[41] * tFact_40 +
																	nrgamma2[42] * tFact_58)))))))))))))));
	return ret_val;
} /* gamma2rpipi_ */

/* --------------------------------------------------------------------------- */
/*second derivative of gamma with respect to pi region 2 */
double H2O::gamma2PiPi(double pNorm, double tNorm)
{
	double ret_val;
	double gamma0;

	gamma0 = -1. / pNorm / pNorm;
	ret_val = gamma0 + gamma2rPiPi(pNorm, tNorm);
	return ret_val;
} /* gamma2pipi_ */

/* ---------------------------------------------------------------------------- */
/* derivative of gamma with respect to tau region 2 */
double H2O::gamma2rTau(double pNorm, double tNorm)
{
	double ret_val;
	static double tFact, pFact_2, pFact_4, tFact_2, tFact_3,
		tFact_4, tFact_5, tFact_6, tFact_7, tFact_9, tFact_10,
		tFact_20, tFact_12, tFact_13, tFact_15, tFact_24,
		tFact_25, tFact_34, tFact_19, tFact_28, tFact_35,
		tFact_38, tFact_39, tFact_47, tFact_49, tFact_52,
		tFact_56, tFact_57;

	tFact = tNorm - .5;
	tFact_2 = tFact * tFact;
	tFact_3 = tFact * tFact_2;
	tFact_4 = tFact_2 * tFact_2;
	tFact_5 = tFact_4 * tFact;
	tFact_6 = tFact_5 * tFact;
	tFact_7 = tFact_6 * tFact;
	tFact_9 = tFact_7 * tFact_2;
	tFact_10 = tFact_9 * tFact;
	tFact_12 = tFact_10 * tFact_2;
	tFact_13 = tFact_12 * tFact;
	tFact_15 = tFact_13 * tFact_2;
	tFact_19 = tFact_10 * tFact_9;
	tFact_20 = tFact_19 * tFact;
	tFact_24 = tFact_20 * tFact_4;
	tFact_25 = tFact_24 * tFact;
	tFact_28 = tFact_25 * tFact_3;
	tFact_34 = tFact_28 * tFact_6;
	tFact_35 = tFact_34 * tFact;
	tFact_38 = tFact_35 * tFact_3;
	tFact_39 = tFact_38 * tFact;
	tFact_47 = tFact_38 * tFact_9;
	tFact_49 = tFact_47 * tFact_2;
	tFact_52 = tFact_49 * tFact_3;
	tFact_56 = tFact_52 * tFact_4;
	tFact_57 = tFact_56 * tFact;
	pFact_2 = pNorm * pNorm;
	pFact_4 = pFact_2 * pFact_2;
	ret_val = pNorm * (nrgamma2[1] + nrgamma2[2] * 2. * tFact + nrgamma2[3] *
		3. * tFact_2 + nrgamma2[4] * 6. * tFact_5 + pNorm * (
			nrgamma2[5] + nrgamma2[6] * 2. * tFact + nrgamma2[7] * 4. *
			tFact_3 + nrgamma2[8] * 7. * tFact_6 + nrgamma2[9] * 36. *
			tFact_35 + pNorm * (nrgamma2[11] + nrgamma2[12] * 3. *
				tFact_2 + nrgamma2[13] * 6. * tFact_5 + nrgamma2[14] * 35. *
				tFact_34 + pNorm * (nrgamma2[15] + nrgamma2[16] * 2. * tFact +
					nrgamma2[17] * 3. * tFact_2 + pNorm * (nrgamma2[18] * 7. *
						tFact_6 + pNorm * (nrgamma2[19] * 3. * tFact_2 + nrgamma2[20]
							* 16. * tFact_15 + nrgamma2[21] * 35. * tFact_34 + pNorm * (
								nrgamma2[23] * 11. * tFact_10 + nrgamma2[24] * 25. * tFact_24
								+ pNorm * (nrgamma2[25] * 8. * tFact_7 + nrgamma2[26] * 36. *
									tFact_35 + pNorm * (nrgamma2[27] * 13. * tFact_12 + pNorm *
										(nrgamma2[28] * 4. * tFact_3 + nrgamma2[29] * 10. * tFact_9 +
											nrgamma2[30] * 14. * tFact_13 + pFact_4 * pFact_2 * (
												nrgamma2[31] * 29. * tFact_28 + nrgamma2[32] * 50. * tFact_49
												+ pFact_2 * (nrgamma2[33] * 57. * tFact_56 + pFact_2 * (
													nrgamma2[34] * 20. * tFact_19 + nrgamma2[35] * 35. * tFact_34
													+ nrgamma2[36] * 48. * tFact_47 + pNorm * (nrgamma2[37] * 21. *
														tFact_20 + pNorm * (nrgamma2[38] * 53. * tFact_52 + pNorm *
															(nrgamma2[39] * 39. * tFact_38 + pNorm * (nrgamma2[40] * 26. *
																tFact_25 + nrgamma2[41] * 40. * tFact_39 + nrgamma2[42] *
																58. * tFact_57)))))))))))))))));
	return ret_val;
} /* gamma2rtau_ */

/* ------------------------------------------------------------------------- */
/* derivative of gamma with respect to tau region 2 */
double H2O::gamma2Tau(double pNorm, double tNorm)
{
	double ret_val;
	static double gamma0, ktNorm;

	ktNorm = 1. / tNorm;
	gamma0 = tNorm * tNorm * (n0gamma2[8] * 3. + ktNorm * (n0gamma2[7] * 2.
		+ ktNorm * (n0gamma2[1] + ktNorm * ktNorm * (n0gamma2[6] * -1. +
			ktNorm * (n0gamma2[5] * -2. + ktNorm * (n0gamma2[4] * -3. +
				ktNorm * (n0gamma2[3] * -4. + ktNorm * (n0gamma2[2] * -5.))))))));
	ret_val = gamma0 + gamma2rTau(pNorm, tNorm);
	return ret_val;
} /* gamma2tau_ */

/* --------------------------------------------------------------------------- */
/*second  derivative of gamma with respect to tau region 2 */
double H2O::gamma2rTauTau(double pNorm, double tNorm)
{
	double ret_val;
	static double tFact, pFact_2, pFact_4, tFact_2, tFact_4,
		tFact_5, tFact_6, tFact_8, tFact_9, tFact_11,
		tFact_12, tFact_14, tFact_23, tFact_24, tFact_33,
		tFact_18, tFact_19, tFact_27, tFact_34, tFact_37,
		tFact_38, tFact_46, tFact_48, tFact_51, tFact_55,
		tFact_56;

	tFact = tNorm - .5;
	tFact_2 = tFact * tFact;
	tFact_4 = tFact_2 * tFact_2;
	tFact_5 = tFact_4 * tFact;
	tFact_6 = tFact_5 * tFact;
	tFact_8 = tFact_6 * tFact_2;
	tFact_9 = tFact_8 * tFact;
	tFact_11 = tFact_9 * tFact_2;
	tFact_12 = tFact_11 * tFact;
	tFact_14 = tFact_12 * tFact_2;
	tFact_18 = tFact_14 * tFact_4;
	tFact_19 = tFact_18 * tFact;
	tFact_23 = tFact_19 * tFact_4;
	tFact_24 = tFact_23 * tFact;
	tFact_27 = tFact_23 * tFact_4;
	tFact_33 = tFact_27 * tFact_6;
	tFact_34 = tFact_33 * tFact;
	tFact_37 = tFact_33 * tFact_4;
	tFact_38 = tFact_37 * tFact;
	tFact_46 = tFact_38 * tFact_8;
	tFact_48 = tFact_46 * tFact_2;
	tFact_51 = tFact_46 * tFact_5;
	tFact_55 = tFact_51 * tFact_4;
	tFact_56 = tFact_55 * tFact;
	pFact_2 = pNorm * pNorm;
	pFact_4 = pFact_2 * pFact_2;
	ret_val = pNorm * (nrgamma2[2] * 2. + nrgamma2[3] * 6. * tFact +
		nrgamma2[4] * 30. * tFact_4 + pNorm * (nrgamma2[6] * 2. +
			nrgamma2[7] * 12. * tFact_2 + nrgamma2[8] * 42. * tFact_5 +
			nrgamma2[9] * 1260. * tFact_34 + pNorm * (nrgamma2[12] * 6. *
				tFact + nrgamma2[13] * 30. * tFact_4 + nrgamma2[14] * 1190. *
				tFact_33 + pNorm * (nrgamma2[16] * 2. + nrgamma2[17] * 6. *
					tFact + pNorm * (nrgamma2[18] * 42. * tFact_5 + pNorm * (
						nrgamma2[19] * 6. * tFact + nrgamma2[20] * 240. * tFact_14 +
						nrgamma2[21] * 1190. * tFact_33 + pNorm * (nrgamma2[23] * 110.
							* tFact_9 + nrgamma2[24] * 600. * tFact_23 + pNorm * (
								nrgamma2[25] * 56. * tFact_6 + nrgamma2[26] * 1260. *
								tFact_34 + pNorm * (nrgamma2[27] * 156. * tFact_11 + pNorm *
									(nrgamma2[28] * 12. * tFact_2 + nrgamma2[29] * 90. * tFact_8
										+ nrgamma2[30] * 182. * tFact_12 + pFact_4 * pFact_2 * (
											nrgamma2[31] * 812. * tFact_27 + nrgamma2[32] * 2450. *
											tFact_48 + pFact_2 * (nrgamma2[33] * 3192. * tFact_55 +
												pFact_2 * (nrgamma2[34] * 380. * tFact_18 + nrgamma2[35] *
													1190. * tFact_33 + nrgamma2[36] * 2256. * tFact_46 + pNorm *
													(nrgamma2[37] * 420. * tFact_19 + pNorm * (nrgamma2[38] *
														2756. * tFact_51 + pNorm * (nrgamma2[39] * 1482. * tFact_37
															+ pNorm * (nrgamma2[40] * 650. * tFact_24 + nrgamma2[41] *
																1560. * tFact_38 + nrgamma2[42] * 3306. * tFact_56))))))))))))
						)))));
	return ret_val;
} /* gamma2rtautau_ */

/* ------------------------------------------------------------------------------- */
/*second  derivative of gamma with respect to tau region 2 */
double H2O::gamma2TauTau(double pNorm, double tNorm)
{
	double ret_val;
	double gamma0, ktNorm;

	ktNorm = 1. / tNorm;
	gamma0 = n0gamma2[7] * 2. + n0gamma2[8] * 6. * tNorm + ktNorm * ktNorm *
		ktNorm * (n0gamma2[6] * 2. + ktNorm * (n0gamma2[5] * 6. + ktNorm *
			(n0gamma2[4] * 12. + ktNorm * (n0gamma2[3] * 20. + ktNorm * 30. *
				n0gamma2[2]))));
	ret_val = gamma0 + gamma2rTauTau(pNorm, tNorm);
	return ret_val;
} /* gamma2tautau_ */

/* ----------------------------------------------------------------------------------- */
double H2O::phi3(double rhoNorm, double tNorm)
{
	double ret_val;
	static double tNorm_2, tNorm_3, tNorm_4, tNorm_6, tNorm_7,
		tNorm_10, tNorm_12, tNorm_22, tNorm_23, tNorm_15,
		tNorm_16, tNorm_17, tNorm_26;

	tNorm_2 = tNorm * tNorm;
	tNorm_3 = tNorm_2 * tNorm;
	tNorm_4 = tNorm_2 * tNorm_2;
	tNorm_6 = tNorm_3 * tNorm_3;
	tNorm_7 = tNorm_6 * tNorm;
	tNorm_10 = tNorm_6 * tNorm_4;
	tNorm_12 = tNorm_10 * tNorm_2;
	tNorm_15 = tNorm_12 * tNorm_3;
	tNorm_16 = tNorm_15 * tNorm;
	tNorm_17 = tNorm_16 * tNorm;
	tNorm_22 = tNorm_10 * tNorm_12;
	tNorm_23 = tNorm_22 * tNorm;
	tNorm_26 = tNorm_23 * tNorm_3;
	ret_val = n3[0] * log(rhoNorm) + n3[1] + n3[2] * tNorm + n3[3] *
		tNorm_2 + n3[4] * tNorm_7 + n3[5] * tNorm_10 + n3[6] *
		tNorm_12 + n3[7] * tNorm_23 + rhoNorm * (n3[8] * tNorm_2 +
			n3[9] * tNorm_6 + n3[10] * tNorm_15 + n3[11] * tNorm_17 +
			rhoNorm * (n3[12] + n3[13] * tNorm_2 + n3[14] * tNorm_6 + n3[15]
				* tNorm_7 + n3[16] * tNorm_22 + n3[17] * tNorm_26 +
				rhoNorm * (n3[18] + n3[19] * tNorm_2 + n3[20] * tNorm_4 + n3[21]
					* tNorm_16 + n3[22] * tNorm_26 + rhoNorm * (n3[23] + n3[24]
						* tNorm_2 + n3[25] * tNorm_4 + n3[26] * tNorm_26 +
						rhoNorm * (n3[27] * tNorm + n3[28] * tNorm_3 + n3[29] *
							tNorm_26 + rhoNorm * (n3[30] + n3[31] * tNorm_2 + n3[32] *
								tNorm_26 + rhoNorm * (n3[33] * tNorm_2 + rhoNorm * (n3[34] *
									tNorm_26 + rhoNorm * (n3[35] * tNorm_2 + n3[36] *
										tNorm_26 + rhoNorm * (n3[37] + n3[38] * tNorm + rhoNorm * n3[39]
											* tNorm_26))))))))));
	return ret_val;
} /* phi3_ */

/* ----------------------------------------------------------------------- */
double H2O::phi3Delta(double rhoNorm, double tNorm)
{
	double ret_val;
	static double tNorm_2, tNorm_3, tNorm_4, tNorm_6, tNorm_7,
		tNorm_10, tNorm_12, tNorm_22, tNorm_23, tNorm_15,
		tNorm_16, tNorm_17, tNorm_26;

	tNorm_2 = tNorm * tNorm;
	tNorm_3 = tNorm_2 * tNorm;
	tNorm_4 = tNorm_2 * tNorm_2;
	tNorm_6 = tNorm_3 * tNorm_3;
	tNorm_7 = tNorm_6 * tNorm;
	tNorm_10 = tNorm_6 * tNorm_4;
	tNorm_12 = tNorm_10 * tNorm_2;
	tNorm_15 = tNorm_12 * tNorm_3;
	tNorm_16 = tNorm_15 * tNorm;
	tNorm_17 = tNorm_16 * tNorm;
	tNorm_22 = tNorm_10 * tNorm_12;
	tNorm_23 = tNorm_22 * tNorm;
	tNorm_26 = tNorm_23 * tNorm_3;
	ret_val = n3[0] / rhoNorm + n3[8] * tNorm_2 + n3[9] * tNorm_6 +
		n3[10] * tNorm_15 + n3[11] * tNorm_17 + rhoNorm * ((n3[12] +
			n3[13] * tNorm_2 + n3[14] * tNorm_6 + n3[15] * tNorm_7 + n3[16]
			* tNorm_22 + n3[17] * tNorm_26) * 2. + rhoNorm * ((n3[18] +
				n3[19] * tNorm_2 + n3[20] * tNorm_4 + n3[21] * tNorm_16 +
				n3[22] * tNorm_26) * 3. + rhoNorm * ((n3[23] + n3[24] *
					tNorm_2 + n3[25] * tNorm_4 + n3[26] * tNorm_26) * 4. +
					rhoNorm * ((n3[27] * tNorm + n3[28] * tNorm_3 + n3[29] *
						tNorm_26) * 5. + rhoNorm * ((n3[30] + n3[31] * tNorm_2 +
							n3[32] * tNorm_26) * 6. + rhoNorm * (n3[33] * 7. * tNorm_2 +
								rhoNorm * (n3[34] * 8. * tNorm_26 + rhoNorm * ((n3[35] *
									tNorm_2 + n3[36] * tNorm_26) * 9. + rhoNorm * ((n3[37] +
										n3[38] * tNorm) * 10. + rhoNorm * 11. * n3[39] * tNorm_26))))))))
			);
	return ret_val;
} /* phi3delta_ */

/* ---------------------------------------------------------------------- */
double H2O::phi3Tau(double rhoNorm, double tNorm)
{
	double ret_val;
	static double tNorm_1, tNorm_2, tNorm_3, tNorm_5, tNorm_6,
		tNorm_9, tNorm_11, tNorm_21, tNorm_22, tNorm_14,
		tNorm_15, tNorm_16, tNorm_25;

	tNorm_2 = tNorm * tNorm;
	tNorm_3 = tNorm_2 * tNorm;
	tNorm_5 = tNorm_3 * tNorm_2;
	tNorm_6 = tNorm_5 * tNorm;
	tNorm_9 = tNorm_6 * tNorm_3;
	tNorm_11 = tNorm_9 * tNorm_2;
	tNorm_14 = tNorm_11 * tNorm_3;
	tNorm_15 = tNorm_14 * tNorm;
	tNorm_16 = tNorm_15 * tNorm;
	tNorm_21 = tNorm_15 * tNorm_6;
	tNorm_22 = tNorm_21 * tNorm;
	tNorm_25 = tNorm_22 * tNorm_3;
	tNorm_1 = tNorm * 2.;
	tNorm_2 *= 3.;
	tNorm_3 *= 4.;
	tNorm_5 *= 6.;
	tNorm_6 *= 7.;
	tNorm_9 *= 10.;
	tNorm_11 *= 12.;
	tNorm_14 *= 15.;
	tNorm_15 *= 16.;
	tNorm_16 *= 17.;
	tNorm_21 *= 22.;
	tNorm_22 *= 23.;
	tNorm_25 *= 26.;
	ret_val = n3[2] + n3[3] * tNorm_1 + n3[4] * tNorm_6 + n3[5] *
		tNorm_9 + n3[6] * tNorm_11 + n3[7] * tNorm_22 + rhoNorm * (
			n3[8] * tNorm_1 + n3[9] * tNorm_5 + n3[10] * tNorm_14 + n3[11]
			* tNorm_16 + rhoNorm * (n3[13] * tNorm_1 + n3[14] *
				tNorm_5 + n3[15] * tNorm_6 + n3[16] * tNorm_21 + n3[17] *
				tNorm_25 + rhoNorm * (n3[19] * tNorm_1 + n3[20] * tNorm_3
					+ n3[21] * tNorm_15 + n3[22] * tNorm_25 + rhoNorm * (n3[24] *
						tNorm_1 + n3[25] * tNorm_3 + n3[26] * tNorm_25 + rhoNorm *
						(n3[27] + n3[28] * tNorm_2 + n3[29] * tNorm_25 + rhoNorm * (
							n3[31] * tNorm_1 + n3[32] * tNorm_25 + rhoNorm * (n3[33] *
								tNorm_1 + rhoNorm * (n3[34] * tNorm_25 + rhoNorm * (n3[35] *
									tNorm_1 + n3[36] * tNorm_25 + rhoNorm * (n3[38] + rhoNorm *
										n3[39] * tNorm_25))))))))));
	return ret_val;
} /* phi3tau_ */

/* -------------------------------------------------------------------- */
double H2O::phi3TauTau(double rhoNorm, double tNorm)
{
	double ret_val;
	static double tNorm_1, tNorm_2, tNorm_4, tNorm_5, tNorm_8,
		tNorm_10, tNorm_20, tNorm_21, tNorm_13, tNorm_14,
		tNorm_15, tNorm_24;

	tNorm_2 = tNorm * tNorm;      /* **2 */
	tNorm_4 = tNorm_2 * tNorm_2;  /* **4 */
	tNorm_5 = tNorm_4 * tNorm;    /* **5 */
	tNorm_8 = tNorm_4 * tNorm_4;  /* **8 */
	tNorm_10 = tNorm_8 * tNorm_2; /* **10 */
	tNorm_13 = tNorm_8 * tNorm_5; /* **13 */
	tNorm_14 = tNorm_13 * tNorm;  /* **14 */
	tNorm_15 = tNorm_14 * tNorm;  /* **15 */
	tNorm_20 = tNorm_15 * tNorm_5;/* **20 */
	tNorm_21 = tNorm_20 * tNorm;  /* **21 */
	tNorm_24 = tNorm_20 * tNorm_4;/* **24 */
	tNorm_1 = tNorm * 6.;
	tNorm_2 *= 12.;
	tNorm_4 *= 30.;
	tNorm_5 *= 42.;
	tNorm_8 *= 90.;
	tNorm_10 *= 132.;
	tNorm_13 *= 210.;
	tNorm_14 *= 240.;
	tNorm_15 *= 272.;
	tNorm_20 *= 462.;
	tNorm_21 *= 506.;
	tNorm_24 *= 650.;
	ret_val = n3[3] * 2. + n3[4] * tNorm_5 + n3[5] * tNorm_8 + n3[6] *
		tNorm_10 + n3[7] * tNorm_21 + rhoNorm * (n3[8] * 2. + n3[9] *
			tNorm_4 + n3[10] * tNorm_13 + n3[11] * tNorm_15 + rhoNorm
			* (n3[13] * 2. + n3[14] * tNorm_4 + n3[15] * tNorm_5 + n3[16]
				* tNorm_20 + n3[17] * tNorm_24 + rhoNorm * (n3[19] * 2. +
					n3[20] * tNorm_2 + n3[21] * tNorm_14 + n3[22] * tNorm_24 +
					rhoNorm * (n3[24] * 2. + n3[25] * tNorm_2 + n3[26] * tNorm_24
						+ rhoNorm * (n3[28] * tNorm_1 + n3[29] * tNorm_24 + rhoNorm
							* (n3[31] * 2. + n3[32] * tNorm_24 + rhoNorm * (n3[33] * 2. +
								rhoNorm * (n3[34] * tNorm_24 + rhoNorm * (n3[35] * 2. + n3[36]
									* tNorm_24 + rhoNorm * rhoNorm * n3[39] * tNorm_24)))))))));
	return ret_val;
} /* phi3tautau_ */

/* --------------------------------------------------------------------- */
double H2O::phi3DeltaTau(double rhoNorm, double tNorm)
{
	double ret_val;
	static double tNorm_2, tNorm_3, tNorm_4, tNorm_6, tNorm_7,
		tNorm_10, tNorm_12, tNorm_22, tNorm_23, tNorm_15,
		tNorm_16, tNorm_17, tNorm_26;

	tNorm_3 = tNorm * tNorm;
	tNorm_4 = tNorm_3 * tNorm;
	tNorm_6 = tNorm_4 * tNorm * tNorm;
	tNorm_7 = tNorm_6 * tNorm;
	tNorm_10 = tNorm_7 * tNorm_4;
	tNorm_12 = tNorm_10 * tNorm_3;
	tNorm_15 = tNorm_12 * tNorm_4;
	tNorm_16 = tNorm_15 * tNorm;
	tNorm_17 = tNorm_16 * tNorm;
	tNorm_22 = tNorm_17 * tNorm_6;
	tNorm_23 = tNorm_22 * tNorm;
	tNorm_26 = tNorm_23 * tNorm_4;
	tNorm_2 = tNorm * 2.;
	tNorm_3 *= 3.;
	tNorm_4 *= 4.;
	tNorm_6 *= 6.;
	tNorm_7 *= 7.;
	tNorm_10 *= 10.;
	tNorm_12 *= 12.;
	tNorm_15 *= 15.;
	tNorm_16 *= 16.;
	tNorm_17 *= 17.;
	tNorm_22 *= 22.;
	tNorm_23 *= 23.;
	tNorm_26 *= 26.;
	ret_val = n3[8] * tNorm_2 + n3[9] * tNorm_6 + n3[10] * tNorm_15 +
		n3[11] * tNorm_17 + rhoNorm * ((n3[13] * tNorm_2 + n3[14] *
			tNorm_6 + n3[15] * tNorm_7 + n3[16] * tNorm_22 + n3[17] *
			tNorm_26) * 2. + rhoNorm * ((n3[19] * tNorm_2 + n3[20] *
				tNorm_4 + n3[21] * tNorm_16 + n3[22] * tNorm_26) * 3. +
				rhoNorm * ((n3[24] * tNorm_2 + n3[25] * tNorm_4 + n3[26] *
					tNorm_26) * 4. + rhoNorm * ((n3[27] + n3[28] * tNorm_3 + n3[29]
						* tNorm_26) * 5. + rhoNorm * ((n3[31] * tNorm_2 + n3[32]
							* tNorm_26) * 6. + rhoNorm * (n3[33] * 7. * tNorm_2 +
								rhoNorm * (n3[34] * 8. * tNorm_26 + rhoNorm * ((n3[35] *
									tNorm_2 + n3[36] * tNorm_26) * 9. + rhoNorm * (n3[38] * 10.
										+ rhoNorm * 11. * n3[39] * tNorm_26)))))))));
	return ret_val;
} /* phi3deltatau_ */

/* ------------------------------------------------------------------ */
double H2O::phi3DeltaDelta(double rhoNorm, double tNorm)
{
	double ret_val;
	static double tNorm_2, tNorm_3, tNorm_4, tNorm_6, tNorm_7,
		tNorm_10, tNorm_12, tNorm_22, tNorm_23, tNorm_15,
		tNorm_16, tNorm_26;

	tNorm_2 = tNorm * tNorm;
	tNorm_3 = tNorm_2 * tNorm;
	tNorm_4 = tNorm_2 * tNorm_2;
	tNorm_6 = tNorm_3 * tNorm_3;
	tNorm_7 = tNorm_6 * tNorm;
	tNorm_10 = tNorm_6 * tNorm_4;
	tNorm_12 = tNorm_10 * tNorm_2;
	tNorm_15 = tNorm_12 * tNorm_3;
	tNorm_16 = tNorm_15 * tNorm;
	tNorm_22 = tNorm_10 * tNorm_12;
	tNorm_23 = tNorm_22 * tNorm;
	tNorm_26 = tNorm_23 * tNorm_3;
	ret_val = -n3[0] / rhoNorm / rhoNorm + (n3[12] + n3[13] * tNorm_2 +
		n3[14] * tNorm_6 + n3[15] * tNorm_7 + n3[16] * tNorm_22 +
		n3[17] * tNorm_26) * 2. + rhoNorm * ((n3[18] + n3[19] *
			tNorm_2 + n3[20] * tNorm_4 + n3[21] * tNorm_16 + n3[22] *
			tNorm_26) * 6. + rhoNorm * ((n3[23] + n3[24] * tNorm_2 + n3[25]
				* tNorm_4 + n3[26] * tNorm_26) * 12. + rhoNorm * ((n3[27]
					* tNorm + n3[28] * tNorm_3 + n3[29] * tNorm_26) * 20. +
					rhoNorm * ((n3[30] + n3[31] * tNorm_2 + n3[32] * tNorm_26) *
						30. + rhoNorm * (n3[33] * 42. * tNorm_2 + rhoNorm * (n3[34] *
							56. * tNorm_26 + rhoNorm * ((n3[35] * tNorm_2 + n3[36] *
								tNorm_26) * 72. + rhoNorm * ((n3[37] + n3[38] * tNorm) * 90.
									+ rhoNorm * 110. * n3[39] * tNorm_26))))))));
	return ret_val;
} /* phi3deltadelta_ */

/**   ******************************************************
*  Calculation of pressure from density and temperature for region  3
*  Parameter:
*  P3_H2O: Pressure in MPa
*  tK    : Temperature in K
*  rho   : Density in kg/m3
 ****************************************************** */
double H2O::p3_H2O(double tK, double rho)
{
	double ret_val;
	double tNorm, rhoNorm;

	rhoNorm = rho / 322.;
	tNorm = 647.096 / tK;
	ret_val = rhoNorm * rho * .461526 * tK *
		phi3Delta(rhoNorm, tNorm) * .001;
	return ret_val;
} /* p3_H2O__ */

/**
* \brief  Approximation of density in region 3 (water sub-region)
*
* Approximation of density in region 3 to give start value for iteration\n
* Boundaries: above 165.292 bar, below critical pressure,
*             above 350 degC, below tsatC
*
* \param pBar pressure [bar]
* \param TempK temperature [K]
* \return rho density [kg/m3]                          */
/* ------------------------------------------------------- */
double H2O::approx3_WATER(double pBar, double TempK)
{
	double x1, z1, z2, z3, z4, rho, TempC;

	TempC = TempK - 273.15;
	x1 = log(pBar);
	z1 = x1 * (x1 * .4900453016111941 - 5.496985340097434) +
		744.4969839686646;
	z2 = TempC * (TempC * (TempC * -1.440177149457524e-5 + .01600869442473704) -
		5.922278143963672);
	z3 = x1 * (x1 * 9.412755797690677e-4 - .0104541834490465) + 1.;
	z4 = TempC * (TempC * (TempC * -1.860355519379978e-8 + 2.090553662234177e-5) -
		.007812763220552083);
	rho = (z1 + z2) / (z3 + z4);
	return rho;
} /* approx3_WATER__ */

/**
* \brief  Approximation of density in region 3 (steam sub-region)
*
* Approximation of density in region 3 to give start value for iteration\n
* Boundaries: above 165.292 bar, below critical pressure,
*             above tsatC, below Boundary to region 2
* \param pBar pressure [bar]
* \param TempK temperature [K]
* \return rho density [kg/m3]                          */
/* ------------------------------------------------------- */
double H2O::approx3_STEAM(double pBar, double TempK)
{
	double x1, y1, z1, z2, z3, z4, rho;

	x1 = log(pBar);
	y1 = log(TempK - 273.15);
	z1 = x1 * 3.684149332679729 - 398.3392644916269;
	z2 = y1 * (y1 * -14.10375588932333 + 147.2371696470623);
	z3 = x1 * .06386412164467998 + 1.;
	z4 = y1 * -.2278839823653747;
	rho = (z1 + z2) / (z3 + z4);
	return rho;
} /* approx3_STEAM__ */

/* -------------------------------------------------------------------- */
double H2O::enth(double tK, double pMPa, int medium)
{
	double ret_val;
	double rho;
	double pNorm, tNorm, tBoundary;
	double rhoNorm;

	if (pMPa <= 16.5292) {
		if (medium == STEAM) {
			/*  calculation of enthalpy for superheated steam */
			/*  Region  2 */
			pNorm = pMPa;
			tNorm = 540. / tK;
			ret_val = tNorm * gamma2Tau(pNorm, tNorm) * GasConstH2O * tK;
			return ret_val;
		}
		else if (medium == WATER) {
			/* calculation of enthalpy for water */
			/* region 1 */
			if (tK < 273.15) {
				cout << "\nError in H2O.enthalpy t < 0 degC" << endl;
				prot << "\nError in H2O.enthalpy t < 0 degC" << endl;
				exit(1);
				//            return 0.;
			}
			pNorm = pMPa / 16.53;
			tNorm = 1386. / tK;
			ret_val = tNorm * gamma1Tau(pNorm, tNorm) * GasConstH2O * tK;
			return ret_val;
		}
		else {
			cout << "\nH2O.enth:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else if (pMPa < 22.064) {
		if (medium == STEAM) {
			/*  calculation of enthalpy for superheated steam */
			/*  checking boundary -> region 2 or 3 */
			tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
				572.54459862746;
			if (tK > tBoundary) {
				/*     region 2 */
				pNorm = pMPa;
				tNorm = 540. / tK;
				ret_val = tNorm * gamma2Tau(pNorm, tNorm) * GasConstH2O * tK;
				return ret_val;
			}
			else {
				/*   region 3 only valid above 350 degC */
				rho = 1. / specVol(tK, pMPa, STEAM);
				rhoNorm = rho / 322.;
				tNorm = 647.096 / tK;
				ret_val = GasConstH2O * tK * (tNorm * phi3Tau(rhoNorm, tNorm)
					+ rhoNorm * phi3Delta(rhoNorm, tNorm));
				return ret_val;
			}
		}
		else if (medium == WATER) {
			/* calculation of enthalpy for water */
			if (tK < 273.15) {
				cout << "\nError in H2O.enthalpy t < 0 degC" << endl;
				prot << "\nError in H2O.enthalpy t < 0 degC" << endl;
				exit(1);
				//            return 0.;
			}
			else if (tK < 623.15) {  // 350+273.15
		 /* region 1 */
				pNorm = pMPa / 16.53;
				tNorm = 1386. / tK;
				ret_val = tNorm * gamma1Tau(pNorm, tNorm) * GasConstH2O * tK;
				return ret_val;
			}
			else {
				/*  region 3 */
				rho = 1. / specVol(tK, pMPa, WATER);
				rhoNorm = rho / 322.;
				tNorm = 647.096 / tK;
				ret_val = GasConstH2O * tK * (tNorm * phi3Tau(rhoNorm, tNorm)
					+ rhoNorm * phi3Delta(rhoNorm, tNorm));
				return ret_val;
			}
		}
		else {
			cout << "\nH2O.enth:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else {
		/* Fluid */
		tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
			572.54459862746;
		if (tK > tBoundary) {
			/*   Even if supercritical still region 2 */
			pNorm = pMPa;
			tNorm = 540. / tK;
			ret_val = tNorm * gamma2Tau(pNorm, tNorm) * GasConstH2O * tK;
			return ret_val;
		}
		else if (tK >= 623.15) {// 350+273.15
		/* region 3 */
			rho = 1. / specVol(tK, pMPa, FLUID);
			rhoNorm = rho / 322.;
			tNorm = 647.096 / tK;
			ret_val = GasConstH2O * tK * (tNorm * phi3Tau(rhoNorm, tNorm)
				+ rhoNorm * phi3Delta(rhoNorm, tNorm));
			return ret_val;
		}
		else if (tK >= 273.15) {
			/*  below 350 degC calculated as water (region 1) */
			pNorm = pMPa / 16.53;
			tNorm = 1386. / tK;
			ret_val = tNorm * gamma1Tau(pNorm, tNorm) * GasConstH2O * tK;
			return ret_val;
		}
		else {
			cout << "\nError in H2O.enthalpy t < 0 degC" << endl;
			exit(1);
		}
	}
	//   return ret_val;
} /* enth_H2O */

/* ---------------------------------------------------------------- */

double H2O::entropy(double tK, double pMPa, int medium)
{
	double ret_val;
	double rho;
	double pNorm, tNorm;
	double tBoundary;
	double rhoNorm;

	if (pMPa <= 16.5292) {
		if (medium == STEAM) {
			/*  calculation of entropy for superheated steam */
			/*  Region  2 */
			pNorm = pMPa;
			tNorm = 540. / tK;
			ret_val = GasConstH2O * (tNorm * gamma2Tau(pNorm, tNorm) -
				gamma2(pNorm, tNorm));
			return ret_val;
		}
		else if (medium == WATER) {
			/* calculation of entropy for water */
		  /* region 1 */
			if (tK < 273.15) {
				cout << "\nError in H2O.entropy t < 0 degC" << endl;
				return 0.;
			}
			pNorm = pMPa / 16.53;
			tNorm = 1386. / tK;
			ret_val = GasConstH2O * (tNorm * gamma1Tau(pNorm, tNorm) -
				gamma1(pNorm, tNorm));
			return ret_val;
		}
		else {
			cout << "\nH2O.entropy:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else if (pMPa < 22.064) {
		if (medium == STEAM) {
			/*  calculation of entropy for superheated steam */
			/*  checking boundary -> region 2 or 3 */
			tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
				572.54459862746;
			if (tK > tBoundary) {
				/* region 2 */
				pNorm = pMPa;
				tNorm = 540. / tK;
				ret_val = GasConstH2O * (tNorm * gamma2Tau(pNorm, tNorm) -
					gamma2(pNorm, tNorm));
				return ret_val;
			}
			else {
				/*   region 3 only valid above 350 deg C */
				rho = 1. / specVol(tK, pMPa, STEAM);
				rhoNorm = rho / 322.;
				tNorm = 647.096 / tK;
				ret_val = GasConstH2O * (tNorm * phi3Tau(rhoNorm, tNorm) -
					phi3(rhoNorm, tNorm));
				return ret_val;
			}
		}
		else if (medium == WATER) {
			/* entropy of water */
			if (tK < 273.15) {
				cout << "\nError in H2O.entropy t < 0 degC" << endl;
				return 0.;
			}
			else if (tK < 623.15) {
				/*  region 1 */
				pNorm = pMPa / 16.53;
				tNorm = 1386. / tK;
				ret_val = GasConstH2O * (tNorm * gamma1Tau(pNorm, tNorm) -
					gamma1(pNorm, tNorm));
				return ret_val;
			}
			else {
				/* region 3 */
				rho = 1. / specVol(tK, pMPa, WATER);
				rhoNorm = rho / 322.;
				tNorm = 647.096 / tK;
				ret_val = GasConstH2O * (tNorm * phi3Tau(rhoNorm, tNorm) -
					phi3(rhoNorm, tNorm));
				return ret_val;
			}
		}
		else {
			cout << "\nH2O.entropy:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else {
		/* Fluid */
		tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
			572.54459862746;
		if (tK > tBoundary) {
			/*   even if supercritical still region 2 */
			pNorm = pMPa;
			tNorm = 540. / tK;
			ret_val = GasConstH2O * (tNorm * gamma2Tau(pNorm, tNorm) -
				gamma2(pNorm, tNorm));
			return ret_val;
		}
		else if (tK >= 623.15) {
			/* region 3 */
			rho = 1. / specVol(tK, pMPa, FLUID);
			rhoNorm = rho / 322.;
			tNorm = 647.096 / tK;
			ret_val = GasConstH2O * (tNorm * phi3Tau(rhoNorm, tNorm) -
				phi3(rhoNorm, tNorm));
			return ret_val;
		}
		else if (tK >= 273.15) {
			/*   below 350 deg C calculated as water (region 1) */
			pNorm = pMPa / 16.53;
			tNorm = 1386. / tK;
			ret_val = GasConstH2O * (tNorm * gamma1Tau(pNorm, tNorm) -
				gamma1(pNorm, tNorm));
			return ret_val;
		}
		else {
			cout << "\nError in H2O.entropy t < 0 degC" << endl;
			return 0.;
		}
	}
} /* entropy_H2O__ */

/* ------------------------------------------------------------------- */

double H2O::cp(double tK, double pMPa, int medium)
{
	double ret_val, d__1;
	double rho;
	double pNorm, tNorm;
	double tBoundary, rhoNorm;

	if (pMPa <= 16.5292) {
		if (medium == STEAM) {
			/* calculation of spec. heat capacity at constant pressure (cp) for */
			/* steam region 2 */
			pNorm = pMPa;
			tNorm = 540. / tK;
			ret_val = GasConstH2O * (-tNorm) * tNorm * gamma2TauTau(pNorm, tNorm);
			return ret_val;
		}
		else if (medium == WATER) {
			/* calculation of spec. heat capacity at constant pressure (cp) for */
			/* water region 1 */
			if (tK < 273.15) {
				cout << "\nError in H2O.cp t < 0 degC" << endl;
				ret_val = 4.1868;
				return ret_val;
			}
			pNorm = pMPa / 16.53;
			tNorm = 1386. / tK;
			ret_val = GasConstH2O * (-tNorm) * tNorm * gamma1TauTau(pNorm, tNorm);
			return ret_val;
		}
		else {
			cout << "\nH2O.cp:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else if (pMPa < 22.064) {
		if (medium == STEAM) {
			/* calculation of spec. heat capacity at constant pressure (cp) for */
			/* steam. check boundary -> region 2 or 3 */
			tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
				572.54459862746;
			if (tK > tBoundary) {
				/* region 2 */
				pNorm = pMPa;
				tNorm = 540. / tK;
				ret_val = GasConstH2O * (-tNorm) * tNorm * gamma2TauTau(
					pNorm, tNorm);
				return ret_val;
			}
			else {
				/* region 3 only valid above 350 deg C */
				rho = 1. / specVol(tK, pMPa, STEAM);
				rhoNorm = rho / 322.;
				tNorm = 647.096 / tK;
				/* Computing 2nd power */
				d__1 = rhoNorm * phi3Delta(rhoNorm, tNorm) - rhoNorm *
					tNorm * phi3DeltaTau(rhoNorm, tNorm);
				ret_val = GasConstH2O * (-tNorm * tNorm * phi3TauTau(
					rhoNorm, tNorm) + d__1 * d__1 / (rhoNorm * 2. *
						phi3Delta(rhoNorm, tNorm) + rhoNorm * rhoNorm *
						phi3DeltaDelta(rhoNorm, tNorm)));
				return ret_val;
			}
		}
		else if (medium == WATER) {
			/* calculation of spec. heat capacity at constant pressure (cp) for */
			/* water */
			if (tK < 273.15) {
				cout << "\nError in H2O.cp t < 0 degC" << endl;
				ret_val = 4.1868;
				return ret_val;
			}
			else if (tK < 625.15) {
				/*  region 1 */
				pNorm = pMPa / 16.53;
				tNorm = 1386. / tK;
				ret_val = GasConstH2O * (-tNorm) * tNorm * gamma1TauTau(
					pNorm, tNorm);
				return ret_val;
			}
			else {
				/*  region 3 */
				rho = 1. / specVol(tK, pMPa, WATER);
				rhoNorm = rho / 322.;
				tNorm = 647.096 / tK;
				/* Computing 2nd power */
				d__1 = rhoNorm * phi3Delta(rhoNorm, tNorm) - rhoNorm *
					tNorm * phi3DeltaTau(rhoNorm, tNorm);
				ret_val = GasConstH2O * (-tNorm * tNorm * phi3TauTau(
					rhoNorm, tNorm) + d__1 * d__1 / (rhoNorm * 2. *
						phi3Delta(rhoNorm, tNorm) + rhoNorm * rhoNorm *
						phi3DeltaDelta(rhoNorm, tNorm)));
				return ret_val;
			}
		}
		else {
			cout << "\nH2O.cp:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else {
		/* Fluid */
		tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
			572.54459862746;
		if (tK > tBoundary) {
			/* even if supercritical still region 2 */
			pNorm = pMPa;
			tNorm = 540. / tK;
			ret_val = GasConstH2O * (-tNorm) * tNorm * gamma2TauTau(pNorm, tNorm);
			return ret_val;
		}
		else if (tK >= 623.15) {
			/*   region 3 */
			rho = 1. / specVol(tK, pMPa, FLUID);
			rhoNorm = rho / 322.;
			tNorm = 647.096 / tK;
			/* Computing 2nd power */
			d__1 = rhoNorm * phi3Delta(rhoNorm, tNorm) - rhoNorm * tNorm *
				phi3DeltaTau(rhoNorm, tNorm);
			ret_val = GasConstH2O * (-tNorm * tNorm * phi3TauTau(rhoNorm,
				tNorm) + d__1 * d__1 / (rhoNorm * 2. * phi3Delta(
					rhoNorm, tNorm) + rhoNorm * rhoNorm * phi3DeltaDelta(
						rhoNorm, tNorm)));
			return ret_val;
		}
		else if (tK >= 273.15) {
			/*   below 350 deg C calculated as water (region 1) */
			pNorm = pMPa / 16.53;
			tNorm = 1386. / tK;
			ret_val = GasConstH2O * (-tNorm) * tNorm * gamma1TauTau(pNorm, tNorm);
			return ret_val;
		}
		else {
			cout << "\nError in H2O.cp t < 0 degC" << endl;
			ret_val = 4.1868;
			return ret_val;
		}
	}
	//   return ret_val;
} /* cp_H2O__ */

/* ------------------------------------------------------------------- */

double H2O::cv(double tK, double pMPa, int medium)
{
	double ret_val, d__1;
	double rho;
	double pNorm, tNorm;
	double tBoundary, rhoNorm;

	if (pMPa <= 16.5292) {
		if (medium == STEAM) {
			/* calculation of spec. heat capacity at constant volume (cv) for */
			/* steam region 2 */
			pNorm = pMPa;
			tNorm = 540. / tK;
			d__1 = 1. + pNorm * gamma2rPi(pNorm, tNorm) - tNorm * pNorm * gamma2rPiTau(pNorm, tNorm);
			ret_val = GasConstH2O * ((-tNorm) * tNorm * gamma2TauTau(pNorm, tNorm) - d__1 * d__1 /
				(1. - pNorm * pNorm * gamma2rPiPi(pNorm, tNorm)));
			return ret_val;
		}
		else if (medium == WATER) {
			/* calculation of spec. heat capacity at constant volume (cv) for */
			/* water region 1 */
			if (tK < 273.15) {
				cout << "\nError in H2O.cv t < 0 degC" << endl;
				ret_val = 4.1868;
				return ret_val;
			}
			pNorm = pMPa / 16.53;
			tNorm = 1386. / tK;
			d__1 = gamma1Pi(pNorm, tNorm) - tNorm * gamma1PiTau(pNorm, tNorm);
			ret_val = GasConstH2O * ((-tNorm) * tNorm * gamma1TauTau(pNorm, tNorm) +
				d__1 * d__1 / gamma1PiPi(pNorm, tNorm));
			return ret_val;
		}
		else {
			cout << "\nH2O.cv:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else if (pMPa < 22.064) {
		if (medium == STEAM) {
			/* calculation of spec. heat capacity at constant volume (cv) for */
			/* steam. check boundary -> region 2 or 3 */
			tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
				572.54459862746;
			if (tK > tBoundary) {
				/* region 2 */
				pNorm = pMPa;
				tNorm = 540. / tK;
				d__1 = 1. + pNorm * gamma2rPi(pNorm, tNorm) - tNorm * pNorm * gamma2rPiTau(pNorm, tNorm);
				ret_val = GasConstH2O * ((-tNorm) * tNorm * gamma2TauTau(pNorm, tNorm) - d__1 * d__1 /
					(1. - pNorm * pNorm * gamma2rPiPi(pNorm, tNorm)));
				return ret_val;
			}
			else {
				/* region 3 only valid above 350 deg C */
				rho = 1. / specVol(tK, pMPa, STEAM);
				rhoNorm = rho / 322.;
				tNorm = 647.096 / tK;
				ret_val = GasConstH2O * (-tNorm) * tNorm * phi3TauTau(rhoNorm, tNorm);
				return ret_val;
			}
		}
		else if (medium == WATER) {
			/* calculation of spec. heat capacity at constant volume (cv) for */
			/* water */
			if (tK < 273.15) {
				cout << "\nError in H2O.cv t < 0 degC" << endl;
				ret_val = 4.1868;
				return ret_val;
			}
			else if (tK < 623.15) {
				/*  region 1 */
				pNorm = pMPa / 16.53;
				tNorm = 1386. / tK;
				d__1 = gamma1Pi(pNorm, tNorm) - tNorm * gamma1PiTau(pNorm, tNorm);
				ret_val = GasConstH2O * ((-tNorm) * tNorm * gamma1TauTau(pNorm, tNorm) +
					d__1 * d__1 / gamma1PiPi(pNorm, tNorm));
				return ret_val;
			}
			else {
				/*  region 3 */
				rho = 1. / specVol(tK, pMPa, WATER);
				rhoNorm = rho / 322.;
				tNorm = 647.096 / tK;
				ret_val = GasConstH2O * (-tNorm) * tNorm * phi3TauTau(rhoNorm, tNorm);
				return ret_val;
			}
		}
		else {
			cout << "\nH2O.cv:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else {
		/* Fluid */
		tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
			572.54459862746;
		if (tK > tBoundary) {
			/* even if supercritical still region 2 */
			pNorm = pMPa;
			tNorm = 540. / tK;
			d__1 = 1. + pNorm * gamma2rPi(pNorm, tNorm) - tNorm * pNorm * gamma2rPiTau(pNorm, tNorm);
			ret_val = GasConstH2O * ((-tNorm) * tNorm * gamma2TauTau(pNorm, tNorm) - d__1 * d__1 /
				(1. - pNorm * pNorm * gamma2rPiPi(pNorm, tNorm)));
			return ret_val;
		}
		else if (tK >= 623.15) {
			/*   region 3 */
			rho = 1. / specVol(tK, pMPa, FLUID);
			rhoNorm = rho / 322.;
			tNorm = 647.096 / tK;
			ret_val = GasConstH2O * (-tNorm) * tNorm * phi3TauTau(rhoNorm, tNorm);
			return ret_val;
		}
		else if (tK >= 273.15) {
			/*   below 350 deg C calculated as water (region 1) */
			pNorm = pMPa / 16.53;
			tNorm = 1386. / tK;
			d__1 = gamma1Pi(pNorm, tNorm) - tNorm * gamma1PiTau(pNorm, tNorm);
			ret_val = GasConstH2O * ((-tNorm) * tNorm * gamma1TauTau(pNorm, tNorm) +
				d__1 * d__1 / gamma1PiPi(pNorm, tNorm));
			return ret_val;
		}
		else {
			cout << "\nError in H2O.cv t < 0 degC" << endl;
			ret_val = 4.1868;
			return ret_val;
		}
	}
	//   return ret_val;
} /* cv_H2O__ */

/* ------------------------------------------------------------------ */

double H2O::specVol(double tK, double pMPa, int medium)
{
	double ret_val;
	int i;
	double rho = 0., upperRho, lowerRho, bfakt;
	double Press, pNorm, tNorm, upperInterval;
	double tempSat, upperPress, lowerPress, tBoundary, lowerInterval, PressDiff;

	if (pMPa <= 16.5292) {
		if (medium == STEAM) {
			/* spec. volume for steam */
			/* region  2 */
//         pNorm = pMPa;
			tNorm = 540. / tK;
			//         ret_val = pNorm * gamma2Pi(pNorm, tNorm) * GasConstH2O * tK / 
			//                   (pMPa * 1e3);
			ret_val = gamma2Pi(pMPa, tNorm) * GasConstH2O * tK / 1e3;
			return ret_val;
		}
		else if (medium == WATER) {
			/* spec. volume for water region 1*/
			if (tK < 273.15) {
				cout << "\nError in H2O.specVolume t < 0 degC" << endl;
				prot << "\nError in H2O.specVolume t < 0 degC" << endl;
				ret_val = 0.;
				return ret_val;
			}
			else {
				/* region 1 */
				pNorm = pMPa / 16.53;
				tNorm = 1386. / tK;
				ret_val = pNorm * gamma1Pi(pNorm, tNorm) * GasConstH2O * tK
					/ (pMPa * 1e3);
				return ret_val;
			}
		}
		else {
			cout << "\nH2O.spec.Vol:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else if (pMPa <= 22.064) {
		if (medium == STEAM) {
			/*   check boundary -> region 2 or 3 */
			tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
				572.54459862746;
			if (tK > tBoundary) {
				/* region 2 */
//            pNorm = pMPa;
				tNorm = 540. / tK;
				//            ret_val = pNorm * gamma2Pi(pNorm, tNorm) * GasConstH2O * tK 
				//                      / (pMPa * 1e3);
				ret_val = gamma2Pi(pMPa, tNorm) * GasConstH2O * tK / 1e3;
				return ret_val;
			}
			else {
				/*  region 3 */
				tempSat = satTemp(pMPa);
				if (pMPa < 19.) {
					upperInterval = .326;
					lowerInterval = .1615;
				}
				else if (tK > tempSat + 2.) {
					upperInterval = .29;
					lowerInterval = 2.12;
				}
				else {
					upperInterval = 90.76;
					lowerInterval = .534;
				}
				rho = approx3_STEAM(pMPa * 10., tK);
				Press = p3_H2O(tK, rho);
				PressDiff = Press - pMPa;
				if (FABS(PressDiff) < 1e-6) {
					return 1. / rho;
				}
				if (PressDiff > 0.) {
					upperRho = rho;
					upperPress = PressDiff;
					lowerRho = rho - lowerInterval;
					Press = p3_H2O(tK, lowerRho);
					lowerPress = Press - pMPa;
					if (FABS(lowerPress) < 1e-6) {
						return 1. / lowerRho;
					}
				}
				else {
					lowerRho = rho;
					lowerPress = PressDiff;
					upperRho = rho + upperInterval;
					Press = p3_H2O(tK, upperRho);
					PressDiff = Press - pMPa;
					upperPress = PressDiff;
					if (FABS(PressDiff) < 1e-6) {
						return 1. / upperRho;
					}
				}
				for (i = 1; i <= 1000; ++i) {
					rho = upperRho - upperPress * (upperRho - lowerRho) /
						(upperPress - lowerPress);
					Press = p3_H2O(tK, rho);
					PressDiff = Press - pMPa;
					if (FABS(PressDiff) < 1e-6) {
						break;
					}
					if (PressDiff * upperPress < 0.) {
						lowerRho = upperRho;
						lowerPress = upperPress;
					}
					else {
						bfakt = 1. - PressDiff / upperPress;
						if (bfakt <= 0.) {
							bfakt = .5;
						}
						lowerPress = bfakt * lowerPress;
					}
					upperRho = rho;
					upperPress = PressDiff;
				}
				ret_val = 1. / rho;
				return ret_val;
			}
		}
		else if (medium == WATER) {
			/* spec. volume for water */
			if (tK < 273.15) {
				cout << "\nError in H2O.specVolume t < 0 degC" << endl;
				ret_val = 0.;
				return ret_val;
			}
			else if (tK < 623.15) {
				/* region 1 */
				pNorm = pMPa / 16.53;
				tNorm = 1386. / tK;
				ret_val = pNorm * gamma1Pi(pNorm, tNorm) * GasConstH2O * tK
					/ (pMPa * 1e3);
				return ret_val;
			}
			else {
				/* region 3, WATER */
				upperInterval = 24.6;
				lowerInterval = 17.11;
				rho = approx3_WATER(pMPa * 10., tK);
				Press = p3_H2O(tK, rho);
				PressDiff = Press - pMPa;
				if (FABS(PressDiff) < 1e-6) {
					return 1. / rho;
				}
				if (PressDiff > 0.) {
					upperRho = rho;
					upperPress = PressDiff;
					lowerRho = rho - lowerInterval;
					Press = p3_H2O(tK, lowerRho);
					lowerPress = Press - pMPa;
					if (FABS(lowerPress) < 1e-6) {
						return 1. / lowerRho;
					}
				}
				else {
					lowerRho = rho;
					lowerPress = PressDiff;
					upperRho = rho + upperInterval;
					Press = p3_H2O(tK, upperRho);
					PressDiff = Press - pMPa;
					upperPress = PressDiff;
					if (FABS(PressDiff) < 1e-6) {
						return 1. / upperRho;
					}
				}
				for (i = 1; i <= 1000; ++i) {
					rho = upperRho - upperPress * (upperRho - lowerRho) /
						(upperPress - lowerPress);
					Press = p3_H2O(tK, rho);
					PressDiff = Press - pMPa;
					if (FABS(PressDiff) < 1e-6) {
						break;
					}
					if (PressDiff * upperPress < 0.) {
						lowerRho = upperRho;
						lowerPress = upperPress;
					}
					else {
						bfakt = 1. - PressDiff / upperPress;
						if (bfakt <= 0.) {
							bfakt = .5;
						}
						lowerPress = bfakt * lowerPress;
					}
					upperRho = rho;
					upperPress = PressDiff;
				}
				ret_val = 1. / rho;
				return ret_val;
			}
		}
		else {
			cout << "\nH2O.specVol:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else {
		/* super-critical check for region 1, 2 or 3 */
		tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
			572.54459862746;
		if (tK < 273.15) {
			cout << "\nError in H2O.specVolume t < 0 degC" << endl;
			ret_val = 0.;
			return ret_val;
		}
		else if (tK < 623.15) {
			/* region 1 */
			pNorm = pMPa / 16.53;
			tNorm = 1386. / tK;
			ret_val = pNorm * gamma1Pi(pNorm, tNorm) * GasConstH2O * tK / (
				pMPa * 1e3);
			return ret_val;
		}
		else if (tK <= tBoundary) {
			/* region 3 */
//          pNorm = pMPa;
			tNorm = 540. / tBoundary;
			lowerRho = 1e3 / (gamma2Pi(pMPa, tNorm) * GasConstH2O * tBoundary);
			lowerPress = p3_H2O(tK, lowerRho) - pMPa;
			pNorm = pMPa / 16.53;
			tNorm = 2.2242192765670636;
			upperRho = pMPa * 1e3 / (pNorm * gamma1Pi(pNorm, tNorm) *
				GasConstH2O * 623.15);
			upperPress = p3_H2O(tK, upperRho) - pMPa;
			for (i = 1; i <= 1000; ++i) {
				rho = upperRho - upperPress * (upperRho - lowerRho) /
					(upperPress - lowerPress);
				Press = p3_H2O(tK, rho);
				PressDiff = Press - pMPa;
				if (FABS(PressDiff) < 1e-6) {
					break;
				}
				if (PressDiff * upperPress < 0.) {
					lowerRho = upperRho;
					lowerPress = upperPress;
				}
				else {
					bfakt = 1. - PressDiff / upperPress;
					if (bfakt <= 0.) {
						bfakt = .5;
					}
					lowerPress = bfakt * lowerPress;
				}
				upperRho = rho;
				upperPress = PressDiff;
			}
			ret_val = 1. / rho;
			return ret_val;
		}
		else {
			/* region 2 */
//          pNorm = pMPa;
			tNorm = 540. / tK;
			ret_val = gamma2Pi(pMPa, tNorm) * GasConstH2O * tK / 1e3;
			//         ret_val = pNorm * gamma2Pi(pNorm, tNorm) * GasConstH2O * tK / (
			//                    pMPa * 1e3);
			return ret_val;
		}
	}
	//    return ret_val;
} /* specVol_H2O__ */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     Function for calculation saturation pressure of water   */
/*     at given saturation temperature                         */
/*     Function of 1997 IFC Formulation for industrial use     */
/*     Input: saturation temperature in K                      */
/*     Return: saturation pressure in MPa                      */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
double H2O::satPress(double tK)
{
	/* Initialized data */
	static double nSat[10] = { 1167.0521452767,-724213.16703206,
		-17.073846940092,12020.82470247,-3232555.0322333,14.91510861353,
		-4823.2657361591,405113.40542057,-.23855557567849,650.17534844798 };

	double ret_val, d__1;
	double a, b, c, ps, theta;

	//   if(ts > 373.946) return(-1.);
	theta = tK + nSat[8] / (tK - nSat[9]);
	a = (theta + nSat[0]) * theta + nSat[1];
	b = (nSat[2] * theta + nSat[3]) * theta + nSat[4];
	c = (nSat[5] * theta + nSat[6]) * theta + nSat[7];
	ps = c * 2. / (-b + sqrt(b * b - a * 4. * c));
	/* Computing 4th power */
	d__1 = ps, d__1 *= d__1;
	ret_val = d__1 * d__1;
	return ret_val;
} /* satPress_ */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     Function for calculation saturation temperature of water  */
/*     at given saturation  pressure                             */
/*     Function of 1997 IFC Formulation for industrial use       */
/*     Input: saturation pressure in MPa                         */
/*     Return: saturation temperature in K                       */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
double H2O::satTemp(double pMPa)
{
	/* Initialized data */
	static double nSat[10] = { 1167.0521452767,-724213.16703206,
		-17.073846940092,12020.82470247,-3232555.0322333,14.91510861353,
		-4823.2657361591,405113.40542057,-.23855557567849,650.17534844798 };
	double ret_val;
	double d, e, f, g, nd, ts, beta;

	if (pMPa > 22.064) return(-1.);

	beta = sqrt(sqrt(pMPa));
	e = nSat[5] + beta * (nSat[2] + beta);
	f = nSat[6] + beta * (nSat[3] + nSat[0] * beta);
	g = nSat[7] + beta * (nSat[4] + nSat[1] * beta);
	d = g * 2. / (-f - sqrt(f * f - e * 4. * g));
	nd = nSat[9] + d;
	ts = (nd - sqrt(nd * nd - (nSat[8] + nSat[9] * d) * 4.)) / 2.;
	ret_val = ts;
	return ret_val;
} /* satTemp_ */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     calculation of dynamic viscosity of water or steam     */
/*     Function of 1997 IFC Formulation for industrial use    */
/*     Parameters :                                           */
/*     temp:  Temperature in K                                */
/*     specVol : Specific Volume in m3/kg                     */
/*     Return :                                               */
/*     dynamic Viscosity in Pa s                              */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

double H2O::dynVisc(double tK, double specVol)
{
	/* Initialized data */
	static double n0[4] = { 1.,.978197,.579829,-.202354 };
	static double n1[19] = { .5132047,.3205656,-.7782567,.1885447,
		.2151778,.7317883,1.241044,1.476783,-.2818107,-1.070786,-1.263184,
		.1778064,.460504,.2340379,-.4924179,-.0417661,.1600435,-.01578386,
		-.003629481 };

	double ret_val;
	double tau, ksi0, ksi1, tau1, delta, tau1_2, tau1_3,
		tau1_4, tau1_5, delta1;

	tau = 647.226 / tK;
	delta = 1. / (specVol * 317.763);
	ksi0 = 1. / (sqrt(tau) * (n0[0] + tau * (n0[1] + tau * (n0[2] + tau * n0[
		3]))));
	tau1 = tau - 1.;
	delta1 = delta - 1.;
	tau1_2 = tau1 * tau1;
	tau1_3 = tau1 * tau1_2;
	tau1_4 = tau1 * tau1_3;
	tau1_5 = tau1 * tau1_4;
	ksi1 = n1[0] + tau1 * n1[1] + tau1_4 * n1[2] + tau1_5 * n1[3] +
		delta1 * (n1[4] + tau1 * n1[5] + tau1_2 * n1[6] + tau1_3 * n1[
			7] + delta1 * (n1[8] + tau1 * n1[9] + tau1_2 * n1[10] + delta1 *
				(n1[11] + tau1 * n1[12] + tau1_2 * n1[13] + tau1_3 * n1[14]
					+ delta1 * (n1[15] + tau1_3 * n1[16] + delta1 * (tau1 * n1[17]
						+ delta1 * tau1_3 * n1[18])))));
	ret_val = ksi0 * exp(delta * ksi1) * 5.5071e-5;
	return ret_val;
} /* dynVisc_H2O__ */


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     calculation of thermal conductivity of water and steam   */
/*     Function of 1997 IFC Formulation for industrial use      */
/*     Parameter :                                              */
/*     tC:  Temperature in K                                    */
/*     pMPa : Pressure in MPa                                   */
/*     medium : WATER, STEAM or FLUID                           */
/*     Return :                                                 */
/*     thermal conductivity in W/mK                             */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
double H2O::conduct(double tK, double pMPa, int medium)
{
	/* Initialized data */
	static double n0[4] = { 1.,6.978267,2.599096,-.998254 };
	static double n1[30]	/* was [5][6] */ = { 1.3293046,1.7018363,
		5.2246158,8.7127675,-1.8525999,-.40452437,-2.2156845,-10.124111,
		-9.5000611,.9340469,.2440949,1.6511057,4.9874687,4.3786606,0.0,
		.018660751,-.76736002,-.27297694,-.91783782,0.0,-.12961068,
		.37283344,-.43083393,0.0,0.0,.044809953,-.1120316,.13333849 };

	double ret_val;
	double rho, tau;
	double tau1, pib = 0., taub = 0., deltab = 0.;
	double tau11, delta, delta1, pitau1 = 0.;
	double delta12, delta14, lambda0, lambda1, lambda2,
		tBoundary;
	double gPi = 0., gPiPi = 0., gPiTau = 0., pDelta = 0.;
	int region = 0;
	double deltapi = 0.;

	tau = 647.226 / tK;
	rho = 1. / specVol(tK, pMPa, medium);
	delta = rho / 317.763;
	//pi = pMPa / 22.115;
	lambda0 = 1. / (sqrt(tau) * (n0[0] + tau * (n0[1] + tau * (n0[2] + tau *
		n0[3]))));

	tau1 = tau - 1.;
	delta1 = delta - 1.;
	lambda1 = n1[0] + delta1 * (n1[5] + delta1 * (n1[10] + delta1 * (n1[15] +
		delta1 * (n1[20] + delta1 * n1[25])))) + tau1 * (n1[1] + delta1 *
			(n1[6] + delta1 * (n1[11] + delta1 * (n1[16] + delta1 * (n1[21] +
				delta1 * n1[26])))) + tau1 * (n1[2] + delta1 * (n1[7] + delta1 * (
					n1[12] + delta1 * (n1[17] + delta1 * (n1[22] + delta1 * n1[27]))))
					+ tau1 * (n1[3] + delta1 * (n1[8] + delta1 * (n1[13] + delta1 *
						n1[18])) + tau1 * (n1[4] + delta1 * n1[9]))));

	lambda1 = exp(delta * lambda1);

	if (pMPa <= 16.5292) {
		if (medium == STEAM) {
			/* region 2 */
			region = 2;
		}
		else if (medium == WATER) {
			/* region 1 */
			region = 1;
		}
		else {
			cout << "\nH2O.conduct:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else if (pMPa < 22.064) {
		if (medium == STEAM) {
			/*     check boundary -> region 2 or 3 */
			tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
				572.54459862746;
			if (tK > tBoundary) {
				/*   region 2 */
				region = 2;
			}
			else {
				/*   region 3 */
				region = 3;
			}
		}
		else if (medium == WATER) {
			if (tK < 623.15) {
				/*     region 1 */
				region = 1;
			}
			else {
				/*     region 3 */
				region = 3;
			}
		}
		else {
			cout << "\nH2O.conduct:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else {
		/* Fluid */
		tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
			572.54459862746;
		if (tK > tBoundary) {
			/*   even if supercritical still region 2 */
			region = 2;
		}
		else if (tK > 623.15) {
			/*   region 3 */
			region = 3;
		}
		else if (tK >= 273.15) {
			/*   below 350 degC calculated as water (region 1) */
			region = 1;
		}
		else {
			/*      error */
			return 0.;
		}
	}
	if (region == 1) {
		pib = pMPa / 16.53;
		taub = 1386. / tK;
		gPiPi = gamma1PiPi(pib, taub);
		gPi = gamma1Pi(pib, taub);
		pitau1 = (gPi - 1386. / tK * gamma1PiTau(pib, taub)) *
			483.77326610897586 / (gPiPi * tK);
		deltapi = -22.115 * gPiPi / (317.763 * GasConstH2O * tK * gPi * gPi);
	}
	else if (region == 2) {
		pib = pMPa;
		taub = 540. / tK;
		gPiPi = gamma2PiPi(pib, taub);
		gPi = gamma2Pi(pib, taub);
		gPiTau = gamma2PiTau(pib, taub);
		pitau1 = 647.226 * (gPiTau * 540. - gPi * tK) / (22.115 * tK * tK * gPiPi);
		deltapi = -22115. * gPiPi / (317.763 * GasConstH2O * tK * gPi * gPi);
	}
	else if (region == 3) {
		deltab = rho / 322.;
		taub = 647.096 / tK;
		pDelta = phi3Delta(deltab, taub);
		pitau1 = 647.226 * GasConstH2O * rho * rho * (pDelta - 647.096 / tK *
			phi3DeltaTau(deltab, taub)) / (22115. * 322.);
		deltapi = 22115. * 322. / (317.763 * rho * GasConstH2O * tK *
			(pDelta * 2. + rho / 322. * phi3DeltaDelta(deltab, taub)));
	}
	else {
		/*      error */
	}
	tau11 = (1. / tau - 1.);
	delta12 = delta1 * delta1;
	delta14 = delta12 * delta12;
	lambda2 = 7.6262320800000008e-8 / dynVisc(tK, 1. / rho) / (tau *
		tau * delta * delta) * pitau1 * pitau1 * pow(delta * deltapi, 0.4678) *
		sqrt(delta) * exp(tau11 * -18.66 * tau11 - delta14);
	ret_val = (lambda0 * lambda1 + lambda2) * .4945;
	return ret_val;
} /* conduct_H2O__ */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*    calculation of surface tension                     */
/*    at boundary of water and steam                     */
/*    Source: 1997 IFC Formulations for industrial use   */
/*    Parameter :                                        */
/*    tSatK: saturation temperature in K                 */
/*    Return :                                           */
/*    surface tension in mN/m                            */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++ */
double H2O::sigma(double tSatK)
{
	double ret_val, d__1;
	double theta;

	theta = tSatK / 647.096;
	if (theta > 1.) return (-1.);
	d__1 = 1. - theta;
	ret_val = pow(d__1, 1.256) * 235.8 * (1. - d__1 * .625);
	return ret_val;
} /* sigma_H2O__ */

/* ------------------------------------------------------------------- */

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*    calculation of speed of sound                      */
/*    Source: B. Spang, www.cheresources.com             */
/*    Parameter :                                        */
/*    tK:   temperature in K                             */
/*    pMPa: saturation temperature in degC               */
/*    medium Water, Steam or Fluid                       */
/*    Return :                                           */
/*    speed of sound in m/s                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++ */

double H2O::SoS(double tK, double pMPa, int medium)
{
	double ret_val, d__1, d__2;
	int region;
	double rho;
	double pNorm, tNorm;
	double tBoundary, rhoNorm;

	region = 0;
	if (pMPa <= 16.5292) {
		if (medium == STEAM) {
			/* calculation of speed of sound for */
			/* steam region 2 */
			region = 2;
		}
		else if (medium == WATER) {
			/* calculation of speed of sound for */
			/* water region 1 */
			if (tK > 273.15) {
				region = 1;
			}
		}
		else {
			cout << "\nH2O.SoS:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else if (pMPa < 22.064) {
		if (medium == STEAM) {
			/* calculation of speed of sound for */
			/* steam. check boundary -> region 2 or 3 */
			tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
				572.54459862746;
			if (tK > tBoundary) {
				/* region 2 */
				region = 2;
			}
			else {
				/* region 3 only valid above 350 deg C */
				region = 3;
			}
		}
		else if (medium == WATER) {
			/* calculation of speed of sound) for */
			/* water */
			if (tK > 273.15 && tK < 623.15) {
				/*  region 1 */
				region = 1;
			}
			else {
				/*  region 3 */
				region = 3;
			}
		}
		else {
			cout << "\nH2O.SoS:wrong medium" << medium << " P= " << pMPa << " t= " << tK << endl;
			return(-1.);
		}
	}
	else {
		/* Fluid */
		tBoundary = sqrt((pMPa - 13.91883977887) / .0010192970039326) +
			572.54459862746;
		if (tK > tBoundary) {
			/* even if supercritical still region 2 */
			region = 2;
		}
		else if (tK >= 623.15) {
			/*   region 3 */
			region = 3;
		}
		else if (tK >= 273.15) {
			/*   below 350 deg C calculated as water (region 1) */
			region = 1;
		}
	}
	switch (region) {
	case 1:
		pNorm = pMPa / 16.53;
		tNorm = 540. / tK;
		d__1 = gamma1Pi(pNorm, tNorm) - tNorm * gamma1PiTau(pNorm, tNorm);
		d__2 = gamma1Pi(pNorm, tNorm);
		ret_val = sqrt(1000. * GasConstH2O * tK * (d__2 * d__2 / (d__1 * d__1 /
			(tNorm * tNorm * gamma1TauTau(pNorm, tNorm)) - gamma1PiPi(pNorm, tNorm))));
		break;
	case 2:
		pNorm = pMPa;
		tNorm = 540. / tK;
		d__1 = pNorm * gamma2rPi(pNorm, tNorm);
		d__2 = (1. + pNorm * gamma2rPi(pNorm, tNorm) - tNorm * pNorm * gamma2rPiTau(pNorm, tNorm));
		ret_val = sqrt(1000. * GasConstH2O * tK * (1. + 2. * d__1 + d__1 * d__1) /
			((1. - pNorm * pNorm * gamma2rPiPi(pNorm, tNorm)) + d__2 * d__2 / (
				tNorm * tNorm * gamma2TauTau(pNorm, tNorm))));
		break;
	case 3:
		rho = 1. / specVol(tK, pMPa, STEAM);
		rhoNorm = rho / 322.;
		tNorm = 647.096 / tK;
		d__1 = (rhoNorm * phi3Delta(rhoNorm, tNorm) - rhoNorm * tNorm * phi3DeltaTau(rhoNorm, tNorm));
		ret_val = sqrt(1000. * GasConstH2O * tK * (2. * rhoNorm * phi3Delta(rhoNorm, tNorm)
			+ rhoNorm * rhoNorm * phi3DeltaDelta(rhoNorm, tNorm) - d__1 * d__1 /
			(tNorm * tNorm * phi3TauTau(rhoNorm, tNorm))));
		break;
	default:
		cout << "\nError in H2O.SoS t < 0 degC" << endl;
		ret_val = 0.;
		break;
	}
	return ret_val;
} /* SoS_H2O__ */


/* +++++++++++++++++++++++++++++++++++++++++++++ */
/*   function for calculation of temperatures    */
/*   at given pressure and enthalpy              */
/*   return value: tK in K                       */
/*   Source :                                    */
/*   Properties of water and steam  (IAPWS-IF97) */
/* +++++++++++++++++++++++++++++++++++++++++++++ */
double H2O::temp(double ha, double pMPa)
{
	/* Initialized data */
	static double nda[34] = { 1089.8952318288,849.51654495535,
		-107.81748091826,33.153654801263,-7.4232016790248,11.765048724356,
		1.844574935579,-4.1792700549624,6.2478196935812,-17.344563108114,
		-200.58176862096,271.96065473796,-455.11318285818,3091.9688604755,
		252266.40357872,-.0061707422868339,-.31078046629583,
		11.670873077107,128127984.04046,-985549096.23276,2822454697.3002,
		-3594897141.0703,1722734991.3197,-13551.334240775,12848734.66465,
		1.3865724283226,235988.32556514,-13105236.545054,7399.9835474766,
		-551966.9703006,3715408.5996233,19127.72923966,-415351.64835634,
		-62.459855192507 };
	static double ndb[38] = { 1489.5041079516,743.07798314034,
		-97.708318797837,2.4742464705674,-.63281320016026,1.1385952129658,
		-.47811863648625,.0085208123431544,.93747147377932,
		3.3593118604916,3.3809355601454,.16844539671904,.73875745236695,
		-.47128737436186,.15020273139707,-.002176411421975,
		-.021810755324761,-.10829784403677,-.046333324635812,
		7.1280351959551e-5,1.1032831789999e-4,1.8955248387902e-4,
		.0030891541160537,.0013555504554949,2.8640237477456e-7,
		-1.0779857357512e-5,-7.6462712454814e-5,1.4052392818316e-5,
		-3.1083814331434e-5,-1.0302738212103e-6,2.821728163504e-7,
		1.2704902271945e-6,7.3803353468292e-8,-1.1030139238909e-8,
		-8.1456365207833e-14,-2.5180545682962e-11,-1.7565233969407e-18,
		8.6934156344163e-15 };
	static double ndc[23] = { -3236839855524.2,7326335090218.1,
		358250899454.47,-583401318515.9,-10783068217.47,20825544563.171,
		610747.83564516,859777.2253558,-25745.72360417,31081.088422714,
		1208.2315865936,482.19755109255,3.7966001272486,-10.842984880077,
		-.04536417267666,1.4559115658698e-13,1.126159740723e-12,
		-1.7804982240686e-11,1.2324579690832e-7,-1.1606921130984e-6,
		2.7846367088554e-5,-5.9270038474176e-4,.0012918582991878 };
	static double nw[20] = { -238.72489924521,404.21188637945,
		113.49746881718,-5.8457616048039,-1.528548241314e-4,
		-1.0866707695377e-6,-13.391744872602,43.211039183559,
		-54.010067170506,30.535892203916,-6.5964749423638,
		.0093965400878363,1.157364750534e-7,-2.5858641282073e-5,
		-4.0644363084799e-9,6.6456186191635e-8,8.0670734103027e-11,
		-9.3477771213947e-13,5.8265442020601e-15,-1.5020185953503e-17 };

	double tBoundary2_3;
	int i, medium;
	static double upperh, tK, pi, lowerh, pi1, pi3, pi4, pi7, enth_Boundary,
		iterLimit, eta, htK, eta1, eta2, eta3, eta4, eta6, eta7, eta8,
		eta9, eta10, eta11, eta12, eta20, eta32, eta24, eta34, eta40,
		eta18, eta28, eta36, eta38, eta42, eta44, eta16, eta22, hdiff,
		lowerTemp, upperTemp, hdiffu, tSatK, hSatSTEAM, hSatWATER;
	double pBoundary;

	iterLimit = 1e-4;
	/* first check for saturation below critical pressure */
	if (pMPa < 22.064) {
		tSatK = satTemp(pMPa);
		hSatSTEAM = enth(tSatK, pMPa, STEAM);
		hSatWATER = enth(tSatK, pMPa, WATER);
		if (ha >= hSatWATER && ha <= hSatSTEAM) {
			/* Saturation */
			return tSatK;
		}
	}
	if (pMPa < 16.529) {
		if (ha > hSatSTEAM) {
			/* STEAM */
			if (pMPa <= 4.) {
				/* sub-region 2a */
				eta = ha / 2e3;
				pi = pMPa;
				eta1 = eta - 2.1;
				eta2 = eta1 * eta1;
				eta3 = eta2 * eta1;
				eta7 = eta3 * eta3 * eta1;
				eta9 = eta7 * eta2;
				eta11 = eta9 * eta2;
				eta12 = eta11 * eta1;
				eta18 = eta11 * eta7;
				eta20 = eta18 * eta2;
				eta24 = eta20 * eta3 * eta1;
				eta28 = eta24 * eta3 * eta1;
				eta32 = eta28 * eta3 * eta1;
				eta34 = eta32 * eta2;
				eta36 = eta34 * eta2;
				eta38 = eta36 * eta2;
				eta40 = eta38 * eta2;
				eta42 = eta40 * eta2;
				eta44 = eta42 * eta2;
				tK = nda[0] + nda[1] * eta1 + nda[2] * eta2 + nda[3] *
					eta3 + nda[4] * eta7 + nda[5] * eta20 + pi * (nda[6]
						+ nda[7] * eta1 + nda[8] * eta2 + nda[9] * eta3 +
						nda[10] * eta7 + nda[11] * eta9 + nda[12] * eta11 +
						nda[13] * eta18 + nda[14] * eta44 + pi * (nda[15] +
							nda[16] * eta2 + nda[17] * eta7 + nda[18] * eta36 +
							nda[19] * eta38 + nda[20] * eta40 + nda[21] * eta42 +
							nda[22] * eta44 + pi * (nda[23] * eta24 + nda[24] * eta44
								+ pi * (nda[25] * eta12 + nda[26] * eta32 + nda[27] *
									eta44 + pi * (nda[28] * eta32 + nda[29] * eta36 +
										nda[30] * eta42 + pi * (nda[31] * eta34 + nda[32] * eta44
											+ pi * nda[33] * eta28))))));
				lowerTemp = tK - .00911;
				upperTemp = tK + .00923;
				medium = STEAM;
			}
			else {
				/*  check for sub-region 2b or 2c */
				pBoundary = ha * (ha * 1.2809002730136e-4 - .67955786399241) +
					905.84278514723;
				eta = ha / 2e3;
				pi = pMPa;
				if (pi <= pBoundary) {
					/* sub-region 2b */
					pi1 = pi - 2.;
					eta1 = eta - 2.6;
					eta2 = eta1 * eta1;
					eta6 = eta2 * eta2 * eta2;
					eta8 = eta6 * eta2;
					eta12 = eta6 * eta6;
					eta18 = eta12 * eta6;
					eta24 = eta18 * eta6;
					eta28 = eta12 * eta8 * eta8;
					eta40 = eta28 * eta12;
					tK = ndb[0] + ndb[1] * eta1 + ndb[2] * eta2 + ndb[3] *
						eta12 + ndb[4] * eta18 + ndb[5] * eta24 + ndb[6] *
						eta28 + ndb[7] * eta40 + pi1 * (ndb[8] + ndb[9] *
							eta2 + ndb[10] * eta6 + ndb[11] * eta12 + ndb[12]
							* eta18 + ndb[13] * eta24 + ndb[14] * eta28 +
							ndb[15] * eta40 + pi1 * (ndb[16] * eta2 + ndb[17]
								* eta8 + ndb[18] * eta18 + ndb[19] * eta40 + pi1 *
								(ndb[20] * eta1 + ndb[21] * eta2 + ndb[22] *
									eta12 + ndb[23] * eta24 + pi1 * (ndb[24] * eta2 +
										ndb[25] * eta12 + ndb[26] * eta18 + ndb[27] *
										eta24 + ndb[28] * eta28 + ndb[29] * eta40 + pi1 *
										(ndb[30] * eta18 + ndb[31] * eta24 + ndb[32] *
											eta40 + pi1 * (ndb[33] * eta28 + pi1 * (ndb[34] *
												eta2 + ndb[35] * eta28 + pi1 * pi1 * (ndb[36] *
													eta1 + ndb[37] * eta40))))))));
				}
				else {
					/*  sub-region 2c */
					pi1 = pi + 25.;
					eta1 = eta - 1.8;
					pi3 = pi1 * pi1 * pi1;
					pi4 = pi3 * pi1;
					pi7 = pi3 * pi4;
					eta2 = eta1 * eta1;
					eta4 = eta2 * eta2;
					eta8 = eta4 * eta4;
					eta10 = eta8 * eta2;
					eta12 = eta10 * eta2;
					eta16 = eta12 * eta4;
					eta20 = eta16 * eta4;
					eta22 = eta20 * eta2;
					tK = (ndc[0] + ndc[1] * eta4 + pi1 * (ndc[2] + ndc[3] *
						eta2 + pi1 * (ndc[4] + ndc[5] * eta2 + pi3 * (
							ndc[6] + ndc[7] * eta1 + pi1 * (ndc[8] + ndc[9] *
								eta2 + pi1 * (ndc[10] + ndc[11] * eta1 + pi1 * (
									ndc[12] * eta4 + ndc[13] * eta8 + pi1 * (ndc[14] *
										eta4 + pi4 * (ndc[15] + ndc[16] * eta1 + ndc[17]
											* eta4 + ndc[18] * eta10 + ndc[19] * eta12 +
											ndc[20] * eta16 + ndc[21] * eta20 + ndc[22] * eta22)))
									)))))) / pi7;
				}
				lowerTemp = tK - .019;
				upperTemp = tK + .022;
				medium = STEAM;
			}
		}
		else {
			/* below 16.5 MPa only remaining is water */
				/* region 1 */
			eta = ha / 2500.;
			pi = pMPa;
			eta1 = eta + 1.;
			eta2 = eta1 * eta1;
			eta6 = eta2 * eta2 * eta2;
			eta10 = eta6 * eta2 * eta2;
			eta32 = eta10 * eta10 * eta10 * eta2;
			tK = nw[0] + nw[1] * eta1 + nw[2] * eta2 + nw[3] * eta6 +
				nw[4] * eta32 / eta10 + nw[5] * eta32 + pi * (nw[6] +
					nw[7] * eta1 + nw[8] * eta2 + nw[9] * eta2 * eta1 +
					nw[10] * eta2 * eta2 + nw[11] * eta10 + nw[12] *
					eta32 + pi * (nw[13] * eta10 + nw[14] * eta32 + pi * (
						nw[15] * eta10 + eta32 * (nw[16] + pi * (nw[17] + pi *
							(nw[18] + pi * nw[19]))))));
			lowerTemp = tK - .0235;
			upperTemp = tK + .0235;
			medium = WATER;
		}
	}
	else {
		/*     pressure is above 165.29 bar -> check for region */
		enth_Boundary = enth(623.15, pMPa, WATER);
		if (FABS(ha - enth_Boundary) < iterLimit) {
			tK = 623.15;
			return tK;
		}
		if (ha < enth_Boundary) {
			/*  region 1 */
			eta = ha / 2500.;
			pi = pMPa;
			eta1 = eta + 1.;
			eta2 = eta1 * eta1;
			eta6 = eta2 * eta2 * eta2;
			eta10 = eta6 * eta2 * eta2;
			eta32 = eta10 * eta10 * eta10 * eta2;
			tK = nw[0] + nw[1] * eta1 + nw[2] * eta2 + nw[3] * eta6 + nw[4]
				* eta32 / eta10 + nw[5] * eta32 + pi * (nw[6] + nw[7] *
					eta1 + nw[8] * eta2 + nw[9] * eta2 * eta1 + nw[10] * eta2
					* eta2 + nw[11] * eta10 + nw[12] * eta32 + pi * (nw[13] *
						eta10 + nw[14] * eta32 + pi * (nw[15] * eta10 + eta32 * (
							nw[16] + pi * (nw[17] + pi * (nw[18] + pi * nw[19]))))));
			lowerTemp = tK - .0235;
			upperTemp = tK + .0235;
			medium = WATER;
		}
		else {
			/* check for region 2 or 3 */
			tBoundary2_3 = sqrt((pMPa - 13.91883977887) / .0010192970039326)
				+ 572.54459862746;
			enth_Boundary = enth(tBoundary2_3, pMPa, STEAM);
			if (FABS(ha - enth_Boundary) < iterLimit) {
				return tBoundary2_3;
			}
			if (ha > enth_Boundary) {
				/*     region 2, check for sub-region 2b or 2c */
				pBoundary = ha * (ha * 1.2809002730136e-4 - .67955786399241) +
					905.84278514723;
				eta = ha / 2e3;
				pi = pMPa;
				if (pi <= pBoundary) {
					/* sub-region 2b  */
					pi1 = pi - 2.;
					eta1 = eta - 2.6;
					eta2 = eta1 * eta1;
					eta6 = eta2 * eta2 * eta2;
					eta8 = eta6 * eta2;
					eta12 = eta6 * eta6;
					eta18 = eta12 * eta6;
					eta24 = eta18 * eta6;
					eta28 = eta12 * eta8 * eta8;
					eta40 = eta28 * eta12;
					tK = ndb[0] + ndb[1] * eta1 + ndb[2] * eta2 + ndb[3] *
						eta12 + ndb[4] * eta18 + ndb[5] * eta24 + ndb[6] *
						eta28 + ndb[7] * eta40 + pi1 * (ndb[8] + ndb[9] *
							eta2 + ndb[10] * eta6 + ndb[11] * eta12 + ndb[12]
							* eta18 + ndb[13] * eta24 + ndb[14] * eta28 +
							ndb[15] * eta40 + pi1 * (ndb[16] * eta2 + ndb[17]
								* eta8 + ndb[18] * eta18 + ndb[19] * eta40 + pi1 *
								(ndb[20] * eta1 + ndb[21] * eta2 + ndb[22] *
									eta12 + ndb[23] * eta24 + pi1 * (ndb[24] * eta2 +
										ndb[25] * eta12 + ndb[26] * eta18 + ndb[27] *
										eta24 + ndb[28] * eta28 + ndb[29] * eta40 + pi1 *
										(ndb[30] * eta18 + ndb[31] * eta24 + ndb[32] *
											eta40 + pi1 * (ndb[33] * eta28 + pi1 * (ndb[34] *
												eta2 + ndb[35] * eta28 + pi1 * pi1 * (ndb[36] *
													eta1 + ndb[37] * eta40))))))));
				}
				else {
					/*  sub-region 2c */
					pi1 = pi + 25.;
					eta1 = eta - 1.8;
					pi3 = pi1 * pi1 * pi1;
					pi4 = pi3 * pi1;
					pi7 = pi3 * pi4;
					eta2 = eta1 * eta1;
					eta4 = eta2 * eta2;
					eta8 = eta4 * eta4;
					eta10 = eta8 * eta2;
					eta12 = eta10 * eta2;
					eta16 = eta12 * eta4;
					eta20 = eta16 * eta4;
					eta22 = eta20 * eta2;
					tK = (ndc[0] + ndc[1] * eta4 + pi1 * (ndc[2] + ndc[3] *
						eta2 + pi1 * (ndc[4] + ndc[5] * eta2 + pi3 * (
							ndc[6] + ndc[7] * eta1 + pi1 * (ndc[8] + ndc[9] *
								eta2 + pi1 * (ndc[10] + ndc[11] * eta1 + pi1 * (
									ndc[12] * eta4 + ndc[13] * eta8 + pi1 * (ndc[14] *
										eta4 + pi4 * (ndc[15] + ndc[16] * eta1 + ndc[17]
											* eta4 + ndc[18] * eta10 + ndc[19] * eta12 +
											ndc[20] * eta16 + ndc[21] * eta20 + ndc[22] * eta22)))
									)))))) / pi7;
				}
				lowerTemp = tK - .019;
				upperTemp = tK + .022;
				medium = STEAM;
			}
			else {
				/*     region 3 */
				if (pMPa < 22.064) {
					if (ha < hSatWATER) {
						lowerTemp = 623.15;
						upperTemp = tSatK;
						medium = WATER;
					}
					else {
						lowerTemp = tSatK;
						upperTemp = tBoundary2_3;
						medium = STEAM;
					}
					while (upperTemp - lowerTemp > 1e-6) {
						tK = (lowerTemp + upperTemp) / 2.;
						htK = enth(tK, pMPa, medium);
						hdiff = htK - ha;
						if (FABS(hdiff) < iterLimit) {
							return tK;
						}
						if (hdiff > 0.) {
							upperTemp = tK;
						}
						else {
							lowerTemp = tK;
						}
					}
				}
				else {
					/* above critical pressure, region 3 */
					lowerTemp = 623.15;
					upperTemp = tBoundary2_3;
					while (upperTemp - lowerTemp > 1e-6) {
						tK = (lowerTemp + upperTemp) / 2.;
						htK = enth(tK, pMPa, FLUID);
						hdiff = htK - ha;
						if (FABS(hdiff) < iterLimit) {
							return tK;
						}
						if (hdiff > 0.) {
							upperTemp = tK;
						}
						else {
							lowerTemp = tK;
						}
					}
				}
				return tK;
			}
		}
	}

	htK = enth(tK, pMPa, medium);
	hdiff = htK - ha;
	if (FABS(hdiff) < iterLimit) {
		return tK;
	}
	if (hdiff > 0.) {
		upperTemp = tK;
		upperh = htK;
		lowerh = enth(lowerTemp, pMPa, medium);
		hdiffu = lowerh - ha;
		if (FABS(hdiffu) < iterLimit) {
			return lowerTemp;
		}
	}
	else {
		lowerTemp = tK;
		lowerh = htK;
		upperh = enth(upperTemp, pMPa, medium);
		hdiff = upperh - ha;
		if (FABS(hdiff) < iterLimit) {
			return upperTemp;
		}
	}
	for (i = 1; i <= 100; ++i) {
		tK = upperTemp - hdiff * (upperTemp - lowerTemp) / (upperh - lowerh);
		htK = enth(tK, pMPa, medium);
		hdiff = htK - ha;
		if (FABS(hdiff) < iterLimit) {
			return tK;
		}
		if (hdiff > 0.) {
			upperTemp = tK;
			upperh = htK;
		}
		else {
			lowerTemp = tK;
			lowerh = htK;
		}
	}
	return tK;
} /* temp_H2O__ */

