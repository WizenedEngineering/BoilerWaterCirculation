//#include "stdafx.h"
#undef MAINFUNCTION

#include "CommonHeader.h"

double TotalDissipation() {
	double dissipation = 0.;
	for (auto& iTube : Tubes) {
		if (iTube.PointOut == DRUM) {
			dissipation += Drum.dpDyn * iTube.Flow / iTube.rhoMean;
		}
		else if (iTube.PointIn != DRUM &&
			(iTube.dpdyn + iTube.dpIn + iTube.dpOut) > 1e-3) { // dissipation energy is 0 for no dyn. pressure drop
			dissipation += 0.5 * (iTube.ksiIn / iTube.rhoIn / iTube.rhoIn +
				iTube.ksiTube / iTube.rhoMean / iTube.rhoMean +
				iTube.ksiOut / iTube.rhoOut / iTube.rhoOut) *
				iTube.MassVel * iTube.MassVel * iTube.Flow * iTube.NoParallel;
		}
		//      prot << "\n tube " << iTube.Number << " dissi " << dissipation<< " dpdyn "<< iTube.dpdyn + iTube.dpIn + iTube.dpOut;
	}
	prot << "\n iterg " << Base.iterg << " dissipation energy " << dissipation;

	//   cout << "\n dissi " << dissipation;
	return dissipation;
}
