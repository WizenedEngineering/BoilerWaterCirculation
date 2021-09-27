//#include "stdafx.h"
#include "CommonHeader.h"

double _tube::FrictFact(double Reynolds, double RelRough) {
	double zeta;
	/*   ----------------------------------------------------- */
	/*   friction factor according Moody's diagram             */
	/*   ----------------------------------------------------- */

	if (Reynolds < 2100.) {
		zeta = 64. / Reynolds;
	}
	else {
		double A = pow((log(1. / (pow((7. / Reynolds), 0.9) + 0.27 * RelRough)) * 2.457), 16.);
		double B = pow((37530. / Reynolds), 16);
		zeta = 8. * pow((1. / pow((A + B), 1.5) + pow((8. / Reynolds), 12.)), (1. / 12.));
	}
	//    prot<<" Frict: reyn "<<Reynolds<< " rough/Diam " << RelRough<<" zeta "<< zeta<< endl;

	return zeta;
}
