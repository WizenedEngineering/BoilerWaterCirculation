/*****************************************************************//**
 * \file   dpBranch.cpp
 * \brief  calculation of pressure difference in branch
 *
 *     Calculation of pressure difference in branch
 *     and coefficients of approximate characteristic curve
 *     dp = dPLinear * g + dPConstant
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

int _branch::dpBranch() {
	/* Local variables */
	size_t PtTbIn,
		PtTbOut,
		iTbPrev,
		iBr0,
		iBr1,
		iBrOff;
	_node* NodeIn = &Nodes[NbNdIn];
	_node* NodeOut = &Nodes[NbNdOut];
	_tube* iTubez;
	_tube* iTubeOff;
	double ksiIn = 0.,
		ksiOut = 0.,
		gz = 0.,
		angle = 0.,
		gfact = 0.,
		gf1 = 0.,
		gf2 = 0.,
		gf3 = 0.,
		gOriginal = 0.,
		betaIn = 0.,
		gcalc = 0.,
		g1 = 0.,
		g2 = 0.,
		dPdyn2 = 0.,
		dPdyn1 = 0.,
		dPstat2 = 0.,
		dPstat1 = 0.;
	//		dpIn1,		dpIn2, dpOut1, dpOut2,
	bool isZero = false;

	if (Base.iterg == 1 && g / minArea < 0.1) {
		/**
		 * --------
		 *  In the first flow iteration step it can happen that the flow is ~0 but we need some flow for the approximate characteristic curve.
		 */
		gOriginal = g;
		if (fabs(deltaH) > 0.1) {
			g = 300. * minArea; /// in vertical branches set maximum mass velocity = 300
		}
		else {
			/**
			 * In horizontal branches (mostly in headers in the middle between 2 supply pipes or connections to drum)\n
			 * the flow is close to 0. A low starting value is used in this case: mass velocity = 50 \n
			 */
			g = 50. * minArea;
		}
		isZero = true;
		/**
		 * This flow is only used for the first approximate characteristic curve and reset later.\n
		 *
		 *  -----------
		 */
	}
	// prot << "\nmassvel " << g/minArea <<    endl;
	if (g / minArea > 0.1) { /// No calculation if max. mass velocity less than 0.1 the data from last flow iteration step are kept
		gf1 = gPrev1 / g;
		gf2 = gPrev2 / g;
		gf3 = gPrev3 / g;
		/**
		 * Looping 3 times through the calculation for 3 flows to determine factors of approximate characteristic curve.
		 */
		for (int iz = 1; iz <= 3; ++iz) {
			switch (iz) {
			case 1:
				/**
				 * The first calculation is done with 90 % of actual flow.\n
				 * If in the last 3 flow iteration steps a flow was closer to the actual flow, this flow is used.\n
				 * This reduces the interval if we are close to final flow solution.
				 */
				gfact = 0.9;
				if (gf1 > 0.9 && gf1 < 0.99) gfact = gf1;
				if (gf2 > gfact && gf2 < 0.99) gfact = gf2;
				if (gf3 > gfact && gf3 < 0.99) gfact = gf3;
				break;
			case 2:
				/**
				* The second calculation is done with 110 % of actual flow.\n
				* If in the last 3 flow iteration steps a flow was closer to the actual flow, this flow is used.\n
				* This reduces the interval if we are close to final flow solution
				*/
				gfact = 1.1;
				if (gf1 < 1.1 && gf1 > 1.01) gfact = gf1;
				if (gf2 < gfact && gf2 > 1.01) gfact = gf2;
				if (gf3 < gfact && gf3 > 1.01) gfact = gf3;
				break;
			case 3:
				/**
				* The last calculation is done with the actual flow.\n
				* All data produced for this branch and even tubes during this calculation are kept
				*/
				gfact = 1.;
			}
			gcalc = gfact * g;
			dPIn = 0.;
			dPOut = 0.;
			dPdyn = 0.;
			dPstat = 0.;
			iTbPrev = MINUS1;
			/**
			 * Iterating through all tubes in this branch to calculate the pressure difference
			 */
			for (size_t iTbInBr = 0; iTbInBr <= mTbInBr; ++iTbInBr) { // iTbInBr needed for easy check for first or last tube in branch
				size_t iTb = NbTbInBr[iTbInBr];
				_tube* iTube = &Tubes[iTb];
				if (Base.showDPBranch || Base.showDPTube || Base.showDPTubeDetail) {
					prot << "\n#Branch " << Number << " iTbInBr " << iTbInBr << " iTb "
						<< iTb << " gcalc " << gcalc << "massvel "
						<< gcalc / iTube->area / iTube->NoParallel << endl;
				}
				PtTbIn = iTube->PointIn;
				PtTbOut = iTube->PointOut;
				betaIn = 0.;
				ksiIn = 0.;
				iTube->dpIn = 0.;
				ksiOut = 0.;
				iTube->dpOut = 0.;
				//    -----------------
				///   1) if tube is drum handle first
				//    -----------------
				if (PtTbIn == DRUM) { // Drum is Point 0 and Node 0
					iTube->dpdyn = 0.;
					iTube->dpIn = 0.;
					iTube->dpOut = 0.;
					dPIn = 0.;
					iTube->dpstat = (Points[PtTbOut].zCoord / 1000. - (100. + Drum.LevelW / 1000.)) *
						9.80665 * Drum.rhoW;
					iTube->pPaIn = Drum.pMPa * 1e6;
					iTube->pPaOut = iTube->pPaIn - iTube->dpstat;
					iTube->EnthIn = Drum.enthW;
					iTube->EnthOut = Drum.enthW;
					iTube->rhoIn = Drum.rhoW;
					iTube->rhoOut = Drum.rhoW;
					iTube->rhoMean = Drum.rhoW;
					if (Base.showDPTube) {
						prot << "\n tube number " << iTube->Number << " drum  dpStat " << iTube->dpstat << "[Pa] dpdyn " << iTube->dpdyn << " [Pa]";
					}
				}
				else if (PtTbOut == DRUM) { // drum riser
				//           prot << "\n\niTb" << iTb << " drum riser";
					iTube->dpdyn = Drum.dpDyn;
					iTube->dpIn = 0.;
					iTube->dpOut = 0.;
					dPOut = 0.;
					iTube->dpstat = (Points[PtTbIn].zCoord / 1000. - (100. + Drum.LevelW / 1000.)) *
						9.80665 * Drum.rhoW;
					iTube->pPaOut = iTube->pPaIn - iTube->dpdyn;
					iTube->EnthIn = Tubes[iTbPrev].EnthOut;
					iTube->EnthOut = Tubes[iTbPrev].EnthOut;
					iTube->rhoIn = Tubes[iTbPrev].rhoOut;
					iTube->rhoOut = Tubes[iTbPrev].rhoOut;
					iTube->rhoMean = Tubes[iTbPrev].rhoOut;
					if (Base.showDPTube) {
						prot << "\n tube number " << iTube->Number << " drum  dpStat " << iTube->dpstat << "[Pa] dpdyn " << iTube->dpdyn << " [Pa]";
					}
				}
				else if (iTube->Dia > MinDrumDiameter) { // tube seems to be mud drum 
					iTube->dpdyn = 0.;
					iTube->dpIn = 0.;
					iTube->dpOut = 0.;
					if (iTbInBr == 0) { // first tube in branch
						dPIn = 0.;
						iTube->pPaIn = NodeIn->pNode + Drum.pMPa * 1e6;
						iTube->EnthIn = enthIn;
						iTube->EnthOut = enthIn;
						double pMPaIn = iTube->pPaIn / 1e6;
						double Dens = 1. / H2O::specVol(H2O::temp(iTube->EnthIn, pMPaIn), pMPaIn, WATER);
						iTube->rhoIn = Dens;
						iTube->rhoMean = Dens;
						iTube->rhoOut = Dens;
						iTube->dpstat = (Points[PtTbOut].zCoord / 1000. - Points[PtTbIn].zCoord / 1000.) *
							9.80665 * Dens;
						iTube->pPaOut = iTube->pPaIn - iTube->dpstat;
						if (Base.showDPTube) {
							prot << "\n pMpaIn " << pMPaIn << " dens " << Dens << " height " << (Points[PtTbOut].zCoord / 1000. - Points[PtTbIn].zCoord / 1000.);
							prot << "\n tube number " << iTube->Number << " mud drum ? dpStat " << iTube->dpstat << "[Pa]";
						}
					}
					else if (iTbInBr == mTbInBr) {
						iTube->pPaIn = Tubes[iTbPrev].pPaOut;
						iTube->EnthIn = Tubes[iTbPrev].EnthOut;
						iTube->EnthOut = iTube->EnthIn;
						dPOut = 0.;
						double pMPaMudDrum = NodeOut->pNode / 1e6 + Drum.pMPa;
						double Dens = 1. / H2O::specVol(H2O::temp(NodeOut->enth, pMPaMudDrum),
							pMPaMudDrum, WATER);
						iTube->dpstat = (Points[PtTbOut].zCoord / 1000. - Points[PtTbIn].zCoord / 1000.) *
							9.80665 * Dens;
						iTube->pPaOut = iTube->pPaIn - iTube->dpstat;
						if (Base.showDPTube) {
							prot << "\n pMpaMuddrum " << pMPaMudDrum << " dens " << Dens << " height " << (Points[PtTbOut].zCoord / 1000. - Points[PtTbIn].zCoord / 1000.);
							prot << "\n tube number " << iTube->Number << " mud drum ? dpStat " << iTube->dpstat << "[Pa]";
						}
					}
					//if (Base.showDPTube) {
					//   prot<<"\n pMpaMuddrum "<<
					//   prot << "\n tube number " << iTube->Number << " mud drum ? dpStat " << iTube->dpstat << "[Pa]";
					//}
				}
				else {
					/**
					 * 2) First tube in branch\n
					 * The pressure difference at inlet of branch is handled in the first tube.\n
					 * This pressure difference is mostly caused by Tee junctions or flow deflections at nodes.\n
					 * A inlet resistance factor ksiIn is used. The standard value is 1.3 as for sharp edged deflection. If there is a Tee this value will be overwritten.
					 */
					if (iTbInBr == 0) { // first tube in branch
						iTube->pPaIn = NodeIn->pNode + Drum.pMPa * 1e6;
						/**
						 * If the steam flow into this branch/tube is greater than 0 the inlet enthalpy is calculated because it depends on the flow.
						 */
						if (gSteamIn < 1e-6) {
							iTube->EnthIn = enthIn;
						}
						else {
							double pMPaIn = iTube->pPaIn * 1e-6;
							double tSatIn = H2O::satTemp(pMPaIn);
							double enthWSatIn = H2O::enth(tSatIn, pMPaIn, WATER);
							double enthSSatIn = H2O::enth(tSatIn, pMPaIn, STEAM);
							iTube->EnthIn = gSteamIn / gcalc * (enthSSatIn - enthWSatIn) + enthWSatIn;
						}
						if (NodeIn->mBrInNd == 1) { // this can only happen if one branch is set to zero
															 // mBrInNd is index, the actual number of branches in this node is 2   
							ksiIn = 0.; // node with 2 branches set to 0
						}
						else {
							ksiIn = 1.3; // standard value for 90 deg sharp-edged deflection will be overwritten

							if (iTube->RadiusBend < 1e-3)
								iTube->beta = 0.; // the first tube can't have an angle to the preceding tube but can be a bend
							//   bool IsT; //Whether Node is a Tee piece
							//   int NbTStraight[2]; // Branch numbers that are straight (inline) at a Tee piece
							//   int NbTOff; // Branch number that branches off at a Tee piece (only one branch)
/**
 * If the node at branch inlet is a Tee more precise inlet resistance factors can be calculated.\n
 * The function checks for the different cases and calls the appropriate case in KsiTee function
 */
							if (NodeIn->IsT) {
								iBr0 = NodeIn->NbBrTStraight[0];
								iBr1 = NodeIn->NbBrTStraight[1];
								iBrOff = NodeIn->NbBrTOff;
								if (Base.showDPBranch) {
									prot << "\nNumber " << Number << " ibr0 " << iBr0 << " ibr1 "
										<< iBr1 << " ibrOff " << iBrOff;
								}
								if (Number == iBrOff) {
									if (Branches[iBr0].NbNdOut == NbNdIn &&
										Branches[iBr1].NbNdOut == NbNdIn) {
										// union of flow 1 and 0
										//	OffUnionInletOff
										angle = fmax(Tubes[NodeIn->NbTbTStraight[0]].Angle3d(*iTube),
											Tubes[NodeIn->NbTbTStraight[1]].Angle3d(*iTube));
										ksiIn = KsiTee(TFlow::OffUnionInletOff, 1.,
											1., 1., 1., 1., 1.0, 1.);
										if (Base.showDPBranch) {
											prot << "inlet flow union from both straight tubes ";
											prot << "ksiIn " << ksiIn << endl;
										}
									}
									else if (Branches[iBr0].NbNdOut == NbNdIn &&
										Branches[iBr1].NbNdIn == NbNdIn) {
										// main flow from iBr0 to iBr1, branching to iBrOff
										// StraightSeparationInletOff, //was 1
										iTubez = &Tubes[NodeIn->NbTbTStraight[0]];
										gz = (Branches[iBr1].g * gfact + gcalc) / iTubez->NoParallel;
										// make sure that total flow is bigger than flow in off branch
										/**
										 * Angle between "straight" branch and "off" branch has to be calculated each time, angle other than 90deg\n
										 * change depending on flow direction in "straight" branch
										 */
										angle = iTubez->Angle3d(*iTube);
										ksiIn = KsiTee(TFlow::StraightSeparationInletOff, gz,
											gcalc / iTube->NoParallel, iTubez->area,
											iTube->area, 0, angle, iTube->Dia);
										if (Base.showDPBranch) {
											prot << " StraightSeparationInletOff:branch-off inlet 1 iTBz " << iTubez->Number;
											prot << "ksiIn " << ksiIn << endl;
										}
									}
									else if (Branches[iBr0].NbNdIn == NbNdIn &&
										Branches[iBr1].NbNdOut == NbNdIn) {
										// main flow from iBr1 to iBr0, branching to iBrOff
										// StraightSeparationInletOff, //was 1
										iTubez = &Tubes[NodeIn->NbTbTStraight[1]];
										gz = (Branches[iBr0].g * gfact + gcalc) / iTubez->NoParallel;
										//make sure that total flow is bigger than flow in off branch
										angle = iTubez->Angle3d(*iTube);
										ksiIn = KsiTee(TFlow::StraightSeparationInletOff, gz,
											gcalc / iTube->NoParallel, iTubez->area,
											iTube->area, 0, angle, iTube->Dia);
										if (Base.showDPBranch) {
											prot << " StraightSeparationInletOff: branch-off inlet 0 iTBz " << iTubez->Number;
											prot << "ksiIn " << ksiIn << endl;
										}
									}
								}
								else { // inlet in straight tube
									if (Branches[iBrOff].NbNdIn == NbNdIn) { /// inlet in straight tube at flow separation
										iTubeOff = &Tubes[NodeIn->NbTbTOff];
										angle = iTubeOff->Angle3d(*iTube);
										//                             if (Number == iBr0) {
										iTubez = &Tubes[NodeIn->NbTbTStraight[1]];
										/*   } else */ if (Number == iBr1) {
											iTubez = &Tubes[NodeIn->NbTbTStraight[0]];
										}
										// make sure that total flow is bigger than flow in off branch
										gz = (Branches[iBrOff].g + g) / iTubez->NoParallel;
										ksiIn = KsiTee(TFlow::StraightSeparationInletStraight, gz,
											Branches[iBrOff].g / iTubeOff->NoParallel,
											iTubez->area, iTubeOff->area, 0, angle,
											iTubeOff->Dia);
										/**
										 * In the case of flow separation the inlet resistance factor of the main flow branch can be negative.\n
										 * To get a dPLinear that is positive, i.e. more flow -> higher dPdyn we have to fix dpIn in this step\n
										 * otherwise we get problems solving the equation system. \n
										 * if we get less pressure drop (or even win) with different flow it has to be handled in next flow iteration step
										 */
										ksiIn /= gfact * gfact;
										if (Base.showDPBranch) {
											prot << " straight iTbOff " << iTubeOff->Number
												<< " StraightSeparationInletStraight ";
											prot << " ksiIn " << ksiIn << endl;
										}
									}
									else if (Branches[iBrOff].NbNdOut == NbNdIn) {
										if (Branches[iBr0].NbNdIn == NbNdIn &&
											Branches[iBr1].NbNdIn == NbNdIn) { // outlet of iBrOff separates to iBr0 and iBr1
									 //  OffSeparationInletStraight
											iTubeOff = &Tubes[NodeIn->NbTbTOff];
											angle = iTubeOff->Angle3d(*iTube);
											ksiIn = KsiTee(TFlow::OffSeparationInletStraight, 1., 1., 1., 1.,
												1., angle, 1.);
											if (Base.showDPBranch) {
												prot << " OffSeparationInletStraight: outlet of iBrOff "
													"separates to iBr0 and iBr1 angle "
													<< ksiIn << endl;
											}
										}
									}
								}
							}
						}
					}
					else {
						/**
						 * 3) Other tubes except last.
						 */
						iTube->pPaIn = Tubes[iTbPrev].pPaOut;
						iTube->EnthIn = Tubes[iTbPrev].EnthOut;
						if (NbNdIn == DRUM && iTbInBr == 1) { // drum downcomers
							ksiIn = 1.5; // first tube is drum, second downcomer
							if (iTube->EnthInGiven < 0.) {
								iTube->EnthIn = Drum.enthW + iTube->EnthInGiven;
							}
							betaIn = iTube->beta;
							iTube->beta = 0.;
						}
						else if (iTbInBr == 1 && Tubes[iTbPrev].Dia > MinDrumDiameter) { // connection to mud drum
							ksiIn = 1.5; // first tube is mud drum, second is leaving
							betaIn = iTube->beta;
							iTube->beta = 0.;
						}
						else if (iTube->RadiusBend < 1e-3 && // no bend
							iTube->Dia < 0.9 * Tubes[iTbPrev].Dia && 
							iTube->beta > 10. && 
							Tubes[iTbPrev].Dia < MinDrumDiameter) { // iTbPrev is header but in one Branch
					 // is handled like Tee with full flow in off tube
							ksiIn = KsiTee(TFlow::StraightSeparationInletOff,
								gcalc / Tubes[iTbPrev].NoParallel,
								gcalc / iTube->NoParallel, Tubes[iTbPrev].area,
								iTube->area, 0, iTube->beta * M_PI / 180., iTube->Dia);
							if (Base.showDPBranch) {
								prot << "\n Header -> tube connection (no Tee) angle "
									<< iTube->beta << " ksiIn " << ksiIn << endl;
							}
							// to avoid any double calculation the angle beta has to be set to
							// zero and later restored
							betaIn = iTube->beta;
							iTube->beta = 0.;
						}
						else if (iTube->RadiusBend < 1e-3 && // no bend
							iTube->Dia > 1.1 * Tubes[iTbPrev].Dia &&
							iTube->beta > 10. &&
							iTube->Dia < MinDrumDiameter) { // iTb is header but in one Branch
					 // is handled like Tee with full flow in incoming tube
					 // a bit tricky: the pressure drop is taken as outlet of previous
					 // but calculated as inlet loss of this tube
							ksiIn = KsiTee(TFlow::StraightUnionOutletOff,
								gcalc / iTube->NoParallel,
								gcalc / Tubes[iTbPrev].NoParallel, iTube->area,
								Tubes[iTbPrev].area, 0, iTube->beta * M_PI / 180.,
								iTube->Dia);
							// now ksi has to be corrected back for main pipe
							double AreaRatio = iTube->area / Tubes[iTbPrev].area;
							ksiIn *= (AreaRatio * AreaRatio);

							if (Base.showDPBranch) {
								prot << "\n  tube connection -> Header(no Tee) angle "
									<< iTube->beta << " ksiIn " << ksiIn << endl;
							}
							// to avoid any double calculation the angle beta has to be set to
							// zero and later restored
							betaIn = iTube->beta;
							iTube->beta = 0.;
						}
					}
					/**
					 * 4) Last tube in branch\n
					 * determination of branch outlet pressure difference
					 */
					 /*     ---------------- */
					 /*     Outlet loss      */
					 /*     ---------------- */
					if (PtTbOut == NbPtOut) {
						if (NodeOut->mBrInNd == 1) { // node with 2 branches set to 0
							ksiOut = 0.;
						}
						else {
							ksiOut = 1.;
							/**
							* If the node at branch outlet is a Tee more precise outlet resistance factors can be calculated.\n
							* The function checks for the different cases and calls the appropriate case in KsiTee function
							*/
							//  bool IsT; Whether Node is a Tee piece
							//  int NbTStraight[2]; Branch numbers that are straight(inline) at a Tee piece
							//  int NbTOff;  Branch number that branches off at a Tee piece (only one branch)
							if (NodeOut->IsT) {
								iBr0 = NodeOut->NbBrTStraight[0];
								iBr1 = NodeOut->NbBrTStraight[1];
								iBrOff = NodeOut->NbBrTOff;
								if (Base.showDPBranch) {
									prot << "\nNumber " << Number << " ibr0 " << iBr0 << " ibr1 "
										<< iBr1 << " ibrOff " << iBrOff;
								}
								if (Number == iBrOff) {
									if (Branches[iBr0].NbNdIn == NbNdOut &&
										Branches[iBr1].NbNdIn == NbNdOut) {
										// separation of flow to 1 and 0, handled by inlet of straight tubes
										//     OffSeparationOutletOff
										ksiOut = 0.;
										if (Base.showDPBranch) {
											prot << " 	OffSeparationOutletOff exit with flow split ";
											prot << "ksiOut " << ksiOut << endl;
										}
									}
									else if (Branches[iBr0].NbNdOut == NbNdOut &&
										Branches[iBr1].NbNdIn == NbNdOut) {
										// main flow from iBr0 to iBr1, union from iBrOff
										iTubez = &Tubes[NodeOut->NbTbTStraight[1]];
										//make sure that total flow is bigger than flow in off branch
										gz = (Branches[iBr0].g + g) / iTubez->NoParallel;
										//										prot << "\n1 is z g.ibr0 " << Branches[iBr0].g << " g.iBr1 "
										//											<< Branches[iBr1].g << " g " << g << " gz " << gz;
										// angle has to be calculated each time, angle other than 90 deg
										// change depending on flow in main branch
										angle = iTubez->Angle3d(*iTube);
										ksiOut = KsiTee(TFlow::StraightUnionOutletOff, gz,
											gcalc / iTube->NoParallel, iTubez->area,
											iTube->area, 0, angle, iTube->Dia);
										if (Base.showDPBranch) {
											prot << " union inlet 1 iTBz " << iTubez->Number
												<< " StraightUnionOutletOff ";
											prot << "ksiOut " << ksiOut << endl;
										}
									}
									else if (Branches[iBr0].NbNdIn == NbNdOut &&
										Branches[iBr1].NbNdOut == NbNdOut) {
										//  main flow from iBr1 to iBr0, union from iBrOff
										iTubez = &Tubes[NodeOut->NbTbTStraight[0]];
										//make sure that total flow is bigger than flow in off branch
										gz = (Branches[iBr1].g + g) / iTubez->NoParallel;
										//										prot << "\n0 is z g.ibr0 " << Branches[iBr0].g << " g.iBr1 "
										//											<< Branches[iBr1].g << " g " << g << " gz " << gz;
										angle = iTubez->Angle3d(*iTube);
										ksiOut = KsiTee(TFlow::StraightUnionOutletOff, gz,
											gcalc / iTube->NoParallel, iTubez->area,
											iTube->area, 0, angle, iTube->Dia);
										if (Base.showDPBranch) {
											prot << " union inlet 0 iTBz " << iTubez->Number
												<< " StraightUnionOutletOff ";
											prot << " ksiOut " << ksiOut << endl;
										}
									}
								}
								else { // outlet from straight tube
									if (Number == iBr0 && NbNdOut == Branches[iBr1].NbNdIn &&
										NbNdOut == Branches[iBrOff].NbNdIn) {
										ksiOut = KsiTee(TFlow::StraightSeparationOutletStraight, 1.,
											1., 1., 1., 0, 1., 1.);
										if (Base.showDPBranch) {
											prot << " StraightSeparationOutletStraight  ksiout "
												<< ksiOut << endl;
										}
									}
									else if (Number == iBr1 && NbNdOut == Branches[iBr0].NbNdIn &&
										NbNdOut == Branches[iBrOff].NbNdIn) {
										ksiOut = KsiTee(TFlow::StraightSeparationOutletStraight, 1.,
											1., 1., 1., 0, 1., 1.);
										if (Base.showDPBranch) {
											prot << " StraightSeparationOutletStraight  ksiout "
												<< ksiOut << endl;
										}
									}
									else if ((Number == iBr0 &&
										NbNdOut == Branches[iBr1].NbNdOut &&
										NbNdOut == Branches[iBrOff].NbNdIn) ||
										(Number == iBr1 &&
											NbNdOut == Branches[iBr0].NbNdOut &&
											NbNdOut == Branches[iBrOff].NbNdIn)) {
										iTubez = &Tubes[NodeOut->NbTbTStraight[0]];
										angle = iTubez->Angle3d(*iTube);
										ksiOut = KsiTee(TFlow::OffUnionOutletStraight, 1., 1., 1.,
											1., 0., angle, 1.);
										if (Base.showDPBranch) {
											prot << " OffUnionOutletStraight  ksiout " << ksiOut
												<< endl;
										}
									}
									else {
										iTubeOff = &Tubes[NodeOut->NbTbTOff];
										if (Branches[iBrOff].NbNdOut == NbNdOut) {
											angle = iTubeOff->Angle3d(*iTube);
											//if (Number == iBr0) {
											iTubez = &Tubes[NodeOut->NbTbTStraight[1]];
											/*} else */if (Number == iBr1) {
												iTubez = &Tubes[NodeOut->NbTbTStraight[0]];
											}
											//make sure that total flow is bigger than flow in off branch
											gz = (Branches[iBrOff].g * gfact + gcalc) / iTubez->NoParallel;
											ksiOut = KsiTee(TFlow::StraightUnionOutletStraight, gz,
												Branches[iBrOff].g * gfact / iTubeOff->NoParallel,
												iTubez->area, iTubeOff->area, 0, angle,
												iTubeOff->Dia);
											if (Base.showDPBranch) {
												prot << " union straight  iTbOff " << iTubeOff->Number
													<< " StraightUnionOutletStraight ";
												prot << " ksiOut " << ksiOut << endl;
											}
										}
									}
								}
							}
						}
					}
					/**
				  * 5) After the inlet and outlet coefficient are set the pressure difference in tube is calculated.\n
				  * The number of sections are determined for the lowest flow (highest enthalpy difference)
				  */

					if (iz == 1) {
						iTube->DetermineNoSections(gcalc);
					}
					iTube->dpTube(ksiIn, ksiOut, gcalc);
					//               iTube->pPaOut = iTube->pPaIn - iTube->dpdyn -
					//               iTube->dpstat - iTube->dpIn - iTube->dpOut;
					if (iTbInBr == 0) dPIn = iTube->dpIn;
					if (iTbInBr == mTbInBr) dPOut = iTube->dpOut;
					if (betaIn != 0.) {
						iTube->beta = betaIn;
						betaIn = 0.;
					}
				}
				dPdyn += iTube->dpdyn + iTube->dpIn + iTube->dpOut;
				dPstat += iTube->dpstat;
				iTbPrev = iTb;
			}
			if (Base.showDPBranch) {
				prot << "\n gcalc " << gcalc << " branch dpdyn " << dPdyn
					<< " branch dpstat " << dPstat << " branch dpIn " << dPIn
					<< " branch dpOut " << dPOut;
			}
			switch (iz) {
			case 1:
				dPdyn1 = dPdyn;
				dPstat1 = dPstat;
				g1 = gcalc;
				break;
			case 2:
				dPdyn2 = dPdyn;
				dPstat2 = dPstat;
				g2 = gcalc;
				break;
			case 3:
				break;
			}
		} // end of loop for the 3 flows
 /**
  * 6) All 3 flows are handled and dPLinear and dPConstant are determined.
  */
		dPLinear = (dPdyn2 + dPstat2 - dPdyn1 - dPstat1) / (g2 - g1);
		dPConstant = (dPdyn + dPstat) - g * dPLinear;

	} //else { // max. mass velocity < 0.1
		// the factors, pressure drop, outlet pressure etc. is kept from previous step
//   }
	if (Base.showDPBranch) {
		double dPtotal = dPdyn + dPstat;
		prot << "\nBranch : " << Number << " g =  " << g << " dPdyn = " << dPdyn
			<< " dPstat = " << dPstat;
		prot << "\ndPtotal " << dPtotal << " dynLin " << dPLinear << " dPConstant "
			<< dPConstant;
	}
	if (isZero) {
		// Flow is ~0, but we needed some flow for the characteristic curve ->has to be set back
		g = gOriginal;
	}

	return 0;
} /* dpBranch */
