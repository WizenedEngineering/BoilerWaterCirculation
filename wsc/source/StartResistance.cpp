/*****************************************************************//*
 * \file   StartResistance.cpp
 * \brief calculates a resistance of branches to be used as "weight" for path finding
 *
 * for the determination of shortest (or split) path the branches have to have a "weight"
 *
 * in this case it is a pseudo resistance factor that is only needed in initFlow function
 *
 * it correspond to the dynamic pressure drop at g = 1 kg/s
 *
 * the branch class member dPConstant is used for this purpose
 *
 * \return int error code
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/
//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

int StartResistance() {

	size_t PtTbIn,
		PtTbOut,
		iTbInBr,
		iTb,
		iTbPrev = 0;
	_tube* iTube;
	double ksiIn,
		ksiOut;
	/// go through all branches
	for (auto& iBranch : Branches) {
		iBranch.dPdyn = 0.;
		iBranch.g = 500. * iBranch.minArea;
		size_t NodeIn = iBranch.NbNdIn;
		/// go through all tubes in branch
		for (iTbInBr = 0; iTbInBr <= iBranch.mTbInBr; ++iTbInBr) {//iTbInBr needed for easy check for first or last tube in branch
			iTb = iBranch.NbTbInBr[iTbInBr];
			iTube = &Tubes[iTb];
			PtTbIn = iTube->PointIn;
			PtTbOut = iTube->PointOut;
			//			betaIn = 0.;
			ksiIn = 0.;
			iTube->dpIn = 0.;
			ksiOut = 0.;
			iTube->dpOut = 0.;
			//    -----------------
			//    Handle drum first
			//    -----------------
///if tube is between drum center and pipe connection to drum (basically inside drum) no dynamic pressure loss 
			if (PtTbIn == DRUM) { // Drum is Point 0 and Node 0
				iTube->dpdyn = 0.;
				iTube->dpstat = (Points[PtTbOut].zCoord / 1000. - (100. + Drum.LevelW / 1000.)) *
					9.80665 * Drum.rhoW;
				iTube->pPaIn = Drum.pMPa * 1e6;
				iTube->pPaOut = iTube->pPaIn - iTube->dpstat;
				iTube->EnthIn = Drum.enthW;
				iTube->EnthOut = Drum.enthW;
			}
			else if (PtTbOut == DRUM) { //drum riser
				/// if tube is between pipe connection to drum and drum center (basically inside drum) then the given drum pressure drop is used
				iBranch.dPdyn += Drum.dpDyn;
			}
			else {
				/**
				* all other tubes:\n
				* 1) Pressure drop at inlet of branch\n
				* handled by friction factor ksiIn for sharp edged deflection
				 */
				if (iTbInBr == 0) { //first tube in branch
					if (Nodes[NodeIn].IsT &&
						(iBranch.Number == Nodes[NodeIn].NbBrTStraight[0] ||
							iBranch.Number == Nodes[NodeIn].NbBrTStraight[1])) {
						ksiIn = 0.2;
					}
					else {
						ksiIn = 1.3; //standard value for sharp edged deflection
					}
					iTube->pPaIn = Drum.pMPa * 1e6;
					iTube->EnthIn = Drum.enthW;
				}
				else {
					/**
					 * 2) if tube is downcomer connected to drum ksiIn is value for outflow from vessel
					 */
					if (NodeIn == DRUM && iTbInBr == 1) { //drum downcomers
						ksiIn = 1.5; //first tube is drum, second downcomer
					}
					else {
						/**
						 * 3) all other tubes use ksiIn like a 90 deg smooth deflection
						 */
						if (iTube->RadiusBend < 1e-3 && Tubes[iTbPrev].RadiusBend < 1e-3) {
							ksiIn = 0.5;
						}
					}
					iTube->pPaIn = Tubes[iTbPrev].pPaOut;
					iTube->EnthIn = Tubes[iTbPrev].EnthOut;
				}
				/*     ---------------- */
				/*     Outlet loss      */
				/*     ---------------- */
				if (iTbInBr == iBranch.mTbInBr) {
					/**
					 * 4) outlet losses are taken as ksiOut = 0.3
					 */
					ksiOut = 0.3;
				}
				iTube->NoSections = 1;
				iTube->LengthSection = iTube->Length;
				iTube->HeightSection = 0.;
				iTube->HeatSection = iTube->q;
				iTube->dpTube(ksiIn, ksiOut, iBranch.g);
				//               iTube->pPaOut = iTube->pPaIn - iTube->dpdyn - iTube->dpstat - iTube->dpIn - iTube->dpOut;
			}
			iBranch.dPdyn += iTube->dpdyn + iTube->dpIn + iTube->dpOut;
			iTbPrev = iTb;
		}
		iBranch.dPConstant = iBranch.dPdyn / iBranch.g / iBranch.g; // should be dpDyn/ g^2 
		if (Base.showMeshDetail) {
			prot << "\n branch dpdyn " << iBranch.dPdyn;
		}
		// reset g
		iBranch.g = 0.;
	}
	return 0;
}
