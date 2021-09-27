/*****************************************************************//**
 * \file   ENTH.cpp
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/

//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

const int MaxIterH = 400; ///< max number of iterations

#if 0
As equal distribution to all leaving branches is standard case
no extra determination of the case needed
/**
* @brief Chawla gives a limit of finely dispersed bubble flow
*
* with this flow pattern no separation of water and steam is likely
*
* Source: VDI Heat Atlas, H3.2
*
* @param x steam quality [-]
* @param rhoWSat Density of water at saturation [kg/m3]
* @param rhoSSat Density of steam at saturation [kg/m3]
* @param Dia Diameter [m]
* @param MassVel mass velocity [kg/m2 s]
* @return bool
*/
bool isChawlaDispersed(double x, double rhoWSat, double rhoSSat, double Dia, double MassVel) {
	double beta_1 = (x * rhoWSat) / ((1. - x) * rhoSSat); // to handle x = 0
	double d__1 = (MassVel * x) / rhoWSat;
	double Fr = d__1 * d__1 * rhoWSat / (rhoSSat * 9.80665 * Dia);
	double srFr = sqrt(Fr);
	return (beta_1 <= 12. * srFr / (1. + srFr / 7.));
}
#endif
/**
* \fn int CalcEnthalpyNodes(void)
* @brief calculation of enthalpy in nodes and inlet steam flow in branches
*
* iterative calculation of enthalpies in nodes
* distribution of steam to different branches according some criteria based on flow pattern
* @return int error code
*/
int CalcEnthalpyNodes(void) {
	/* Local variables */
	size_t iTbArrive, iBrArrive;
	_branch* iBranchLeave;
	double  enthNdPrev, enthWSatNd, enthEvapNd, tSatNd,
		pNd;
	double SeparationFactor = 1.; ///< it is the part of steam flow in node that goes in the branch with less steam. A factor of 1 means equal inlet enthalpy in all branches
	int iterh;

	//prot << "\n start enth drum enthalpy " << Drum.enthW << "\n";
/**
 *-# setting branch inlet enthalpy to drum water enthalpy as start value\n
 *    and branch outlet enthalpy by heat and flow\n
 *    just to assure a minimum value in the node\n
 *    will be overwritten during iteration
 */
	for (auto& iBranch : Branches) {
		iBranch.enthIn = Drum.enthW;
		if (iBranch.g > 1e-5) {
			iBranch.enthOut = Drum.enthW + iBranch.qSum / iBranch.g;
		}
	}
	/**
	 *-# first, all branches leaving drum have no steam at inlet but can be sub-cooled
	 */
	for (const auto& Leave : Nodes[DRUM].NbBrLeave) {
		Branches[Leave].gSteamIn = 0.;
		size_t iTb = Branches[Leave].NbTbInBr[1];
		if (Tubes[iTb].EnthInGiven < 0.) {
			Branches[Leave].enthIn = Drum.enthW + Tubes[iTb].EnthInGiven;
		}
		if (Branches[Leave].g > 1e-5) {
			Branches[Leave].enthOut = Branches[Leave].enthIn + Branches[Leave].qSum / Branches[Leave].g;
		}
		else {
			Branches[Leave].enthOut = Branches[Leave].enthIn;
		}
	}
	/**
	 *-# flows are constant in this step \n
	 * gSumArriving and gSumLeave in each node are also constant and should be the same
	 */
	for (auto& iNode : Nodes) {
		iNode.gSumArrive = 0.;
		iNode.gSumLeave = 0.;
		for (const auto& Arrive : iNode.NbBrArrive) {
			iNode.gSumArrive += Branches[Arrive].g;
		}
		for (const auto& Leave : iNode.NbBrLeave) {
			iNode.gSumLeave += Branches[Leave].g;
		}
	}
	/**
	 *-# the enthalpy in the node depends on the enthalpies of the nodes connected to this node\n
	 *    certainly only the incoming (arriving) branches and their start nodes\n
	 *    as each node depends on most other nodes there is an iteration loop until the differences are too small
	 */
	for (iterh = 1; iterh <= MaxIterH; ++iterh) {
		double diffSum = 0.;
		/**
		*-# Loop over all nodes except drum
		*/
		for (_node& iNode : Nodes) {
			if (iNode.Number > 0) {
				/**
				*	-# determine enthalpy and steam in node by summing up outlet enthalpy * flow for all arriving branches\n
				*		and divide by sum of incoming flows
				*/
				enthNdPrev = iNode.enth;
				iNode.enth = 0.;
				for (const auto& Arrive : iNode.NbBrArrive) {
					iNode.enth += Branches[Arrive].enthOut * Branches[Arrive].g;
					if (Base.showEnth) {
						prot << "\niterh " << iterh << " iNd " << iNode.Number << " arrive " << Arrive
							<< " Br_enthOut " << Branches[Arrive].enthOut << " Br_g " << Branches[Arrive].g << " enthNode " << iNode.enth;
					}
				}
				if (iNode.enth < 1e-6 || iNode.gSumArrive < 1e-6) {
					iNode.enth = Branches[iNode.NbBrArrive[0]].enthOut;
				}
				else {
					iNode.enth /= iNode.gSumArrive; //mixture enthalpy
				}
				if (Base.showEnth) {
					prot << "\n iNd " << iNode.Number << " enthNd " << iNode.enth << " flowdiff " << iNode.gSumArrive - iNode.gSumLeave;
				}
				/**
				 * -# check if enthalpy is above saturation water enthalpy\n
				 * if yes calculate steam flow in node
				 */
				pNd = iNode.pNode / 1e6 + Drum.pMPa;
				tSatNd = H2O::satTemp(pNd);
				enthWSatNd = H2O::enth(tSatNd, pNd, WATER);
				enthEvapNd = H2O::enth(tSatNd, pNd, STEAM) - enthWSatNd;
				if (iNode.enth > enthWSatNd) {
					iNode.gSteam = (iNode.enth - enthWSatNd) / enthEvapNd * iNode.gSumArrive;
				}
				else {
					iNode.gSteam = 0.;
				}
				if (Base.showEnth) {
					prot << "\niterh " << iterh << " iNd " << iNode.Number << " gsumArrive " << iNode.gSumArrive
						<< " gSteam " << iNode.gSteam << " enthPrev " << enthNdPrev << " enthNode " << iNode.enth;
				}
				diffSum += fabs(iNode.enth - enthNdPrev);
				//         prot << "\n diffsum " << diffSum;
/**
 * -# the mass balance in node is given by the path search for the first iterg-step\n
 * in following iterg-steps the mass balance should also be close\n
 * we only have to deal with the distribution of steam (if there is any)
 *
 * -# setting same inlet enthalpy and corresponding steam flow in all leaving branches\n
 * as standard case
 *
 */
				for (auto& Leave : iNode.NbBrLeave) {
					iBranchLeave = &Branches[Leave];
					iBranchLeave->enthIn = iNode.enth;
					iBranchLeave->gSteamIn = iNode.gSteam * iBranchLeave->g / iNode.gSumLeave;
					if (iBranchLeave->g < 1e-5) {
						iBranchLeave->enthOut = iBranchLeave->enthIn;
					}
					else {
						iBranchLeave->enthOut = iBranchLeave->enthIn + iBranchLeave->qSum / iBranchLeave->g;
					}
				}
				/**
				 * -# easiest case first: only one branch leaving\n
				 * no further checking needed
				 */
				if (iNode.mBrLeave == 0) { // only one leaving branch
					if (Base.showEnth) {
						prot << "\n single leave ";
					}
					continue;
				}
				/**
				 * -# check for special case only if steam flow in node is greater 0
				 */
				if (iNode.gSteam > 0.) {
					for (auto& Leave : iNode.NbBrLeave) {
						iBranchLeave = &Branches[Leave];
						/**
						 * -# if there are more leaving branches but only one arriving and we have steam flow in node the steam is distributed according following criteria
						 *		- equal distribution to all leaving branches (all tubes same inlet enthalpy, already covered):
						 *			- flow pattern is bubble or mist flow for horizontal tubes (Steiner) or finely dispersed bubble (Chawla)
						 *			- distance to the last node on arriving branch is short -> turbulence from mixing at last node prevents separation
						 *			- for Tees: all tubes are horizontal, the arriving branch is TOff and orientation is OffHor, OffVerUp, OffVerDown or Straight is horizontal ad Off is vetical up or down
						 *		- if we have a Tee flow distribution according to orientation and flow pattern
						 *			- orientation straight horizontal Off Vertical up
						 *				- flow pattern stratified or wavy (steam is at top of tube): all steam to Off
						 *				- flow pattern slug/plug or annular: 20 % of steam to Off
						 *			- orientation straight horizontal Off Vertical down
						 *				- flow pattern stratified or wavy (steam is at top of tube): no steam to Off
						 *				- flow pattern slug/plug or annular: 20 % of steam to Off
						 *			- orientation straight vertical off horizontal
						 *				- 30 % of steam to Off
						 *			- straight any orientation off horizontal
						 *				- equal distribution (all tubes same inlet enthalpy)
						 *			- straight any orientation off vertical up
						 *				- 60% of steam goes to off branch
						 *			- straight any orientation off vertical down
						 *				- flow pattern stratified or wavy (steam is at top of tube): all steam to Off
						 *				- flow pattern slug/plug or annular: 20 % of steam to Off
						 *			- all other orientation
						 *				- equal distribution (all tubes same inlet enthalpy, already covered)
						 * -# all other cases: equal distribution (all tubes same inlet enthalpy, already covered)
						 */
						 // distribution of steam to different branches according some criteria (flow pattern, steiner, Chawla)
						if (iNode.mBrArrive == 0) { //we have only one branch arriving but several leaving -> looking for the flow pattern
							iBrArrive = iNode.NbBrArrive[0];
							iTbArrive = Branches[iBrArrive].NbTbInBr[Branches[iBrArrive].mTbInBr];
							Tubes[iTbArrive].Steiner(pNd * 1e6, Tubes[iTbArrive].xOut);
							/* As equal distribution to all leaving branches is standard case
							no extra determination of the case needed
							 this part not need because it's the same as standard case
														if (isChawlaDispersed(Tubes[iTbArrive].xOut, 1. / H2O::specVol(tSatNd, pNd, WATER), 1. / H2O::specVol(tSatNd, pNd, STEAM),
															Tubes[iTbArrive].Dia, Tubes[iTbArrive].MassVel) || // bubbles finely dispersed
															Tubes[iTbArrive].steiner == FlowPattern::bubble || // steiner bubble flow
															Tubes[iTbArrive].steiner == FlowPattern::mist || // steiner mist flow
															(Branches[iBrArrive].mTbInBr == 0 && Tubes[iTbArrive].Length < Tubes[iTbArrive].Dia) ||// last arriving branch to iBrArrive caused mixing
															iNode.TOrientation == TeeOrientation::Undetermined ||
															iNode.TOrientation == TeeOrientation::StraightHorOffHor || // all tubes in Tee are horizontal
															((iNode.TOrientation == TeeOrientation::StraightHorOffVerUp ||
																iNode.TOrientation == TeeOrientation::StraightHorOffVerDown ||
																iNode.TOrientation == TeeOrientation::OffHor ||
																iNode.TOrientation == TeeOrientation::OffVerUp ||
																iNode.TOrientation == TeeOrientation::OffVerDown) && iBrArrive == iNode.NbBrTOff)
															) {//inlet enthalpy in branches is same
															if (Base.showEnth) {

																	prot << "\n same inlet enthalpy";
																	prot << "\nBranch " << Leave << " enthIn " << iBranchLeave->enthIn << " enthOut " << iBranchLeave->enthOut <<
																			" gSteamIn " << iBranchLeave->gSteamIn << " g " << iBranchLeave->g << " gSumLeave " << iNode.gSumLeave;
															}
															continue;
															}
														else */
							if (iNode.IsT) {
								switch (iNode.TOrientation) {
								case TeeOrientation::StraightHorOffVerUp:
									if (Base.showEnth) {
										prot << "\n number " << iNode.Number << " straightHor Off Ver up";
									}
									if (Tubes[iTbArrive].steiner == FlowPattern::stratified || // steiner stratified flow
										Tubes[iTbArrive].steiner == FlowPattern::wavy) { // steiner wavy flow)
										if (Leave == iNode.NbBrTOff) {
											if (iNode.gSteam <= iBranchLeave->g) {
												iBranchLeave->gSteamIn = iNode.gSteam;
											}
											else {
												iBranchLeave->gSteamIn = iBranchLeave->g;
											}
										}
										else {
											if (iNode.gSteam <= Branches[iNode.NbBrTOff].g) {
												iBranchLeave->gSteamIn = 0.;
											}
											else {
												iBranchLeave->gSteamIn = iNode.gSteam - Branches[iNode.NbBrTOff].g;
											}
										}
									}
									else if (Tubes[iTbArrive].steiner == FlowPattern::plugSlug || // steiner slug/plug flow
										Tubes[iTbArrive].steiner == FlowPattern::annular) { // steiner annular flow)
										SeparationFactor = 0.2;
										if (Leave == iNode.NbBrTOff) {
											iBranchLeave->gSteamIn = iNode.gSteam * ((1. - SeparationFactor) + SeparationFactor * iBranchLeave->g / iNode.gSumLeave);
										}
										else {
											iBranchLeave->gSteamIn = iNode.gSteam * SeparationFactor * iBranchLeave->g / iNode.gSumLeave;
										}
									}
									break;
								case TeeOrientation::StraightHorOffVerDown:
									if (Base.showEnth) {
										prot << "\n number " << iNode.Number << " straightHor Off Ver down";
									}
									if (Tubes[iTbArrive].steiner == FlowPattern::stratified || // steiner stratified flow
										Tubes[iTbArrive].steiner == FlowPattern::wavy) { // steiner wavy flow)
										if (Leave == iNode.NbBrTOff) {
											iBranchLeave->gSteamIn = 0.;
										}
										else {
											iBranchLeave->gSteamIn = iNode.gSteam;
										}
									}
									else if (Tubes[iTbArrive].steiner == FlowPattern::plugSlug || // steiner slug/plug flow
										Tubes[iTbArrive].steiner == FlowPattern::annular) { // steiner annular flow)
										SeparationFactor = 0.2;
										if (Leave == iNode.NbBrTOff) {
											iBranchLeave->gSteamIn = iNode.gSteam * SeparationFactor * iBranchLeave->g / iNode.gSumLeave;
										}
										else {
											iBranchLeave->gSteamIn = iNode.gSteam * ((1. - SeparationFactor) + SeparationFactor * iBranchLeave->g / iNode.gSumLeave);
										}
									}
									break;
								case TeeOrientation::StraightVerOffHor:
									if (Base.showEnth) {
										prot << "\n number " << iNode.Number << " straight Ver Off Hor";
									}
									SeparationFactor = 0.3;
									if (Leave == iNode.NbBrTOff) {
										iBranchLeave->gSteamIn = iNode.gSteam * SeparationFactor * iBranchLeave->g / iNode.gSumLeave;
									}
									else {
										iBranchLeave->gSteamIn = iNode.gSteam * ((1. - SeparationFactor) + SeparationFactor * iBranchLeave->g / iNode.gSumLeave);
									}
									break;
								case TeeOrientation::OffHor:
									iBranchLeave->gSteamIn = iNode.gSteam * iBranchLeave->g / iNode.gSumLeave;
									break;
								case TeeOrientation::OffVerUp:
									SeparationFactor = 0.6;
									if (Leave == iNode.NbBrTOff) {
										iBranchLeave->gSteamIn = iNode.gSteam * ((1. - SeparationFactor) + SeparationFactor * iBranchLeave->g / iNode.gSumLeave);
									}
									else {
										iBranchLeave->gSteamIn = iNode.gSteam * SeparationFactor * iBranchLeave->g / iNode.gSumLeave;
									}
									break;
								case TeeOrientation::OffVerDown:
									SeparationFactor = 0.1;
									if (Leave == iNode.NbBrTOff) {
										iBranchLeave->gSteamIn = iNode.gSteam * SeparationFactor * iBranchLeave->g / iNode.gSumLeave;
									}
									else {
										iBranchLeave->gSteamIn = iNode.gSteam * ((1. - SeparationFactor) + SeparationFactor * iBranchLeave->g / iNode.gSumLeave);
									}
									break;
									//									default:  already handled
									//										iBranchLeave->gSteamIn = iNode.gSteam * iBranchLeave->g / iNode.gSumLeave;
									//										break;
								}
							}
							//								else { already handled 
																//TODO: node with  2 leaving branches but not a Tee, maybe same criteria as Tees
							//									iBranchLeave->gSteamIn = iNode.gSteam * iBranchLeave->g / iNode.gSumLeave;
							//								}
							if (iBranchLeave->g < 1e-5) {
								iBranchLeave->enthIn = iNode.enth;
								iBranchLeave->enthOut = iBranchLeave->enthIn;
							}
							else {
								iBranchLeave->enthIn = iBranchLeave->gSteamIn / iBranchLeave->g * enthEvapNd + enthWSatNd;
								iBranchLeave->enthOut = iBranchLeave->enthIn + iBranchLeave->qSum / iBranchLeave->g;
							}
						}
						/*							else { // more than one branch arriving
														if (iBranchLeave->g < 1e-5) {
															iBranchLeave->enthIn = iNode.enth;
														}
														else {
															iBranchLeave->gSteamIn = iNode.gSteam * iBranchLeave->g / iNode.gSumLeave;
															iBranchLeave->enthIn = iBranchLeave->gSteamIn / iBranchLeave->g * enthEvapNd + enthWSatNd;
														}
													} */
						if (Base.showEnth) {
							prot << "\nBranch " << Leave << " enthIn " << iBranchLeave->enthIn << " enthOut " << iBranchLeave->enthOut <<
								" gSteamIn " << iBranchLeave->gSteamIn << " g " << iBranchLeave->g << " gSumLeave " << iNode.gSumLeave;
						}
					}
				}
				if (Base.showEnth) {
					prot << "\n\n";
				}
			}
		} // end nodes-loop
/**
 * -# end of loop over all nodes
 * -# check for condition to break main iteration loop\n
 * sum of all enthalpy differences between this step and previous step should be less than 1 kJ/kg
 * -# end of main iteration loop
 */
		if (Base.showEnth) {
			prot << "\n iterh = " << iterh << " diffSum " << diffSum;
		}
		if (diffSum < 1.) break;
		if (iterh >= MaxIterH - 1) {
			Base.showEnth = true;
			if (iterh == MaxIterH) {
				cout << "\n maximum number of Iterations in enth exceeded" << endl;
				prot << "\n maximum number of Iterations in enth exceeded" << endl;
				return 1;
			}
		}
	} // end iterh-loop

	if (Base.showEnth) {
		prot << "\n    *********************";
		prot << "\n    * Enthalpy function *";
		prot << "\n    *********************\n Node    gSteam enthalpy\n";
		for (size_t iNd = 1; iNd <= mNd; ++iNd) { //without drum
			prot << setw(5) << iNd << " " << setprecision(4) << setw(8) << Nodes[iNd].gSteam << setprecision(4) << setw(8) << Nodes[iNd].enth << endl;
		}
		prot << "\nBranch  gSteamIn  EnthIn   EnthOut " << endl;
		for (const auto& iBranch : Branches) {
			prot << setw(5) << iBranch.Number << "   " << setprecision(4) << setw(8) << iBranch.gSteamIn
				<< "   " << setprecision(4) << setw(8) << iBranch.enthIn
				<< "   " << setprecision(4) << setw(8) << iBranch.enthOut << endl;
		}
		prot << endl;
	}
	return 0;
} /* enth_ */

