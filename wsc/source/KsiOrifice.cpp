//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

// Brandt, F., Dampferzeuger: Kesselsysteme, Energiebilanz, Stroemungstechnik. FDBR Fachbuchreihe Band 3, Vulkan-Verlag Essen, 1992  

double _tube::ksiOrifice(double LengthOrifice, double visc) {
	double ReOrifice, m1, m2, tau, eps, ksi0, ksiphi;
	double DiaOrifice, zeta, ret_val;

	/// If the length of orifice > 0 we have inlet orifice (typical smaller hole in header,\n
	/// length of orifice is header wall thickness)
	if (LengthOrifice > 1e-3) {
		DiaOrifice = DiaOrificeIn;
		m1 = 0.;
		m2 = DiaOrifice / Dia;
		tau = 0.5 * (1. - tanh((LengthOrifice / DiaOrifice - 0.657) * 2.264));
	}
	else {
		DiaOrifice = DiaOrificeOut;
		m1 = DiaOrifice / Dia;
		m2 = m1;
		tau = 1.35;
	}

	ksi0 = 0.5 * (1. - m1) + (1. - m2) * (1. - m2) + tau * sqrt(1. - m1) * (1. - m2);

	ReOrifice = MassVel / (m2 * m2) * DiaOrifice / visc;
	zeta = FrictFact(ReOrifice, Base.Rough / DiaOrifice);

	if (ReOrifice < 1e5) {
		ksiphi = 0.39 * exp(-18.34 * pow(m1, 3.5)) * exp(-1.078 * (log(ReOrifice) - 3.));
		eps = 0.26 * tanh(1.04 * (log(ReOrifice) - 4.)) + 0.74;
		ksi0 = ksiphi + eps * ksi0;
	}
	ret_val = (ksi0 + zeta * fabs(LengthOrifice) / DiaOrifice) * (1. / (m2 * m2));
	return ret_val;
}

