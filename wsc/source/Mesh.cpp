/*****************************************************************//**
 * \file   Mesh.cpp
 * \brief  subroutine to determine the mesh (nodes and branches)
 *
 * \author Rainer_Jordan@<very, very warm>mail.com
 * ## Licence
 * Licensed under the EUPL, Version 1.2 or later
 * \date   September 2021
 *********************************************************************/
//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

int findMaxNodesConnected();

/*
 * \fn _node::isTee(_tube& firstTube, _tube& secondTube, _tube& thirdTube)
 * \brief checks if this node is a Tee and sets the orientation of the Tee
 *
 * All 3 tubes are straight (no bends)\n
 * 2 tubes should be in line (straight) and have the same diameter\n 
 * the remaining tube should have equal or smaller diameter, the angle to the 2 other tubes should be between 45 deg and 135 deg  
 * \param [in] firstTube reference to first tube
 * \param [in] secondTube reference to second tube
 * \param [in] thirdTube reference to third tube
 * \return
 */
bool _node::isTee(_tube& firstTube, _tube& secondTube, _tube& thirdTube) {
	double firstAngle, secondAngle, thirdAngle;

	firstAngle = firstTube.Angle3d(secondTube);
	secondAngle = firstTube.Angle3d(thirdTube);
	thirdAngle = secondTube.Angle3d(thirdTube);
	if (Base.showMeshDetail) {
		prot << "\n iNode " << Number;
		prot << "\n firstTb " << firstTube.Number << " secondTb " << secondTube.Number << " thirdTb " << thirdTube.Number;
		prot << "\n angle1 " << firstAngle * 180. / M_PI << " deg sin " << sin(firstAngle) << " sin(pi/4) " << sin(M_PI_4) <<
			"\n angle2 " << secondAngle * 180. / M_PI << " deg sin " << sin(secondAngle) <<
			"\n angle3 " << thirdAngle * 180. / M_PI << " deg sin " << sin(thirdAngle);
		prot << "\n fabs(firstTube.Dia - secondTube.Dia)" << fabs(firstTube.Dia - secondTube.Dia)
			<< " firstTube.Dia " << firstTube.Dia << " secondTube.Dia " << secondTube.Dia;
		prot << "\n fabs(firstTube.Dia - thirdTube.Dia)" << fabs(firstTube.Dia - thirdTube.Dia)
			<< " firstTube.Dia " << firstTube.Dia << " thirdTube.Dia " << thirdTube.Dia;
		prot << "\n fabs(thirdTube.Dia - secondTube.Dia)" << fabs(thirdTube.Dia - secondTube.Dia)
			<< " thirdTube.Dia " << thirdTube.Dia << " secondTube.Dia " << secondTube.Dia;
	}
	if (sin(firstAngle) < 1e-5 &&//Angle between tubes is either 0 or 180 deg, i.e first and second tubes are inline
		sin(secondAngle) >= 0.7071 && //secondAngle, thirdAngle should be >= 45 deg and =<135 deg
		sin(thirdAngle) >= 0.7071 && //  this is case,  if sin(angle) >= sin(45deg = Pi/4)
		// to allow for small differences (rounding of coordinates) for sin(45deg) a value of 0.7071 is used
		fabs(firstTube.Dia - secondTube.Dia) < 1e-3 &&
		thirdTube.Dia <= firstTube.Dia) {
		IsT = true;
		NbTbTStraight[0] = firstTube.Number;  // Tube numbers that are straight (inline) at a T-piece
		NbTbTStraight[1] = secondTube.Number; // Tube numbers that are straight (inline) at a T-piece
		NbTbTOff = thirdTube.Number;          // Tube number that branches off at a T-piece (only one branch)
		NbBrTStraight[0] = firstTube.NbBr;  // Branch numbers that are straight (inline) at a T-piece
		NbBrTStraight[1] = secondTube.NbBr; // Branch numbers that are straight (inline) at a T-piece
		NbBrTOff = thirdTube.NbBr;          // Branch number that branches off at a T-piece (only one branch)
// ------------------------------------------------
// set orientation of Tee piece
// ----------------------------------------------
		TOrientation = TeeOrientation::Undetermined;
		if (Base.showMeshDetail) {
			prot << "\n isT " << IsT;
		}
		if (fabs(firstTube.Height / firstTube.Length) < 0.17365 &&
			fabs(secondTube.Height / secondTube.Length) < 0.17365) { //straight tubes are horizontal (<10 deg)

			if (fabs(thirdTube.Height / thirdTube.Length) < 0.17365) {  // off tube is horizontal
				TOrientation = TeeOrientation::StraightHorOffHor;
				if (Base.showMeshDetail) {
					prot << " orientation : StraightHorOffHor";
				}
			}
			else if (fabs(thirdTube.Height / thirdTube.Length) > 0.5) { // off tube is vertical

				if ((thirdTube.Height > 1e-3 && thirdTube.PointIn == NbPt) || //off tube is vertical up
					(thirdTube.Height < -1e-3 && thirdTube.PointOut == NbPt)) { //off tube is vertical up
					TOrientation = TeeOrientation::StraightHorOffVerUp;
					if (Base.showMeshDetail) {
						prot << " orientation : StraightHorOffVerUp";
					}
				}
				else if ((thirdTube.Height < -1e-3 && thirdTube.PointIn == NbPt) ||
					(thirdTube.Height > 1e-3 && thirdTube.PointOut == NbPt)) { //off tube is vertical down
					TOrientation = TeeOrientation::StraightHorOffVerDown;
					if (Base.showMeshDetail) {
						prot << " orientation :StraigthHorOffVerDown";
					}
				}
			}
		}
		else if ((fabs(firstTube.Height / firstTube.Length) > 0.866) &&
			(fabs(secondTube.Height / secondTube.Length) > 0.866)) { //straight tubes are vertical (>60 deg)
			if (fabs(thirdTube.Height / thirdTube.Length) < 0.17365) {  // off tube is horizontal
				TOrientation = TeeOrientation::StraightVerOffHor;
				if (Base.showMeshDetail) {
					prot << " orientation : StraightVerOffHor";
				}
			}
		}
		else { //straight tubes neither horizontal nor vertical
			if (fabs(thirdTube.Height / thirdTube.Length) < 0.17365) {  // off tube is horizontal
				TOrientation = TeeOrientation::OffHor;
				if (Base.showMeshDetail) {
					prot << " orientation : OffHor";
				}
			}
			else if ((thirdTube.Height > 1e-3 && thirdTube.PointIn == NbPt) || //off tube is vertical up
				(thirdTube.Height < -1e-3 && thirdTube.PointOut == NbPt)) { //off tube is vertical up
				TOrientation = TeeOrientation::OffVerUp;
				if (Base.showMeshDetail) {
					prot << " orientation : OffVerUp";
				}
			}
			else if ((thirdTube.Height < -1e-3 && thirdTube.PointIn == NbPt) ||
				(thirdTube.Height > 1e-3 && thirdTube.PointOut == NbPt)) { //off tube is vertical down
				TOrientation = TeeOrientation::OffVerDown;
				if (Base.showMeshDetail) {
					prot << " orientation : OffVerDown";
				}
			}
		}
		return true;
	}
	else {
		return false;
	}
}

/**
* \fn void determineTee()
* @brief determines all Tees in the mesh
*
* determine which nodes are Tee (T-piece),\n
* which tubes and branches are inline or branching off at this T-piece\n
* as well as spatial orientation of Tee\n
* here only geometry is checked\n
* condition: straight: same diameter and angle between tubes 0 deg or 180 deg,\n
* off: diameter <= straight, angle between 45 deg and 135 deg
 */
void determineTee() {

	_tube Tube1, Tube2, Tube3;
	for (auto& iNode : Nodes) {
		if (iNode.Number == DRUM) continue;
		_point& iPoint = Points[iNode.NbPt];

		if (iPoint.NoTb == 3) { // T-piece can only have 3 tubes that should belong to 3 branches
			Tube1 = Tubes[iPoint.NbTb[0]];
			Tube2 = Tubes[iPoint.NbTb[1]];
			Tube3 = Tubes[iPoint.NbTb[2]];
			if (iNode.isTee(Tube1, Tube2, Tube3)) {
				continue;
			}
			if (iNode.isTee(Tube2, Tube3, Tube1)) {
				continue;
			}
			if (iNode.isTee(Tube1, Tube3, Tube2)) {
				continue;
			}
		}
	}
	return;
}

/**
 * \brief sets up the mesh (nodes, branches)
 *
 * \return int error code
 */
int Mesh() {
	/* Local variables */
	int j;
	size_t PtTbOut, iPt,
		iTb; // , iBr, iNd;
	double hdiff, angle1; // ,

	Nodes.push_back(_node());
	//! per definition point 0 and node 0 is drum
	Nodes[DRUM].Number = DRUM;
	Nodes[DRUM].NbPt = DRUM;
	Nodes[DRUM].Elev = 100.;
	Points[DRUM].NbNd = DRUM;
	hdiff = 100. - Points[DRUM].zCoord / 1e3;
	mNd = 0;
	/*!     --------------------
	 *     determine nodes
		  -------------------- */
		  /*! all points that have more than 2 tubes connected are nodes */
	for (auto& iPoint : Points) {
		if (iPoint.NoTb >= 3 && iPoint.Number > 0) {
			mNd++;
			Nodes.push_back(_node());
			Nodes[mNd].NbPt = iPoint.Number;
			Nodes[mNd].Elev = iPoint.zCoord / 1e3 + hdiff;
			Nodes[mNd].Number = mNd;
			iPoint.NbNd = mNd;
			continue;
		}
		if (iPoint.NoTb == 1) {
			cout << "\nError: Only one tube to point " << iPoint.Number << " see <projectname>TubesPoints.dxf and check input data";
			prot << "\nError: Only one tube to point " << iPoint.Number << " see <projectname>TubesPoints.dxf and check input data";
			Print2dxf(ShowMode::TubesPoints);
			processError(iPoint.Number, "only one tube to point");
			exit(1);
		}
		if (iPoint.NoTb == 0) {// as the points are defined as starting or end points of tubes this should not happen
			cout << "\nPoint " << iPoint.Number << " is not connected to tube. see <projectname>TubesPoints.dxf and check input data";;
			prot << "\nPoint " << iPoint.Number << " is not connected to tube. see <projectname>TubesPoints.dxf and check input data";;
			Print2dxf(ShowMode::TubesPoints);	
			processError(iPoint.Number, "no tube connected to this point");
			exit(1);
		}
	}
	Nodes.shrink_to_fit();

	if (Base.showMeshDetail) {
		prot << "\n    ******************************************";
		prot << "\n    * Mesh checking before renumber of nodes *";
		prot << "\n    ******************************************";
		prot << "\n mPt " << mPt << " mTb " << mTb << " mNd " << mNd << " mBr " << mBr;

		prot << "\n\n tubes start end ";
		for (auto& iTube : Tubes) {
			prot << "\n" << setw(5) << iTube.Number << " " << setw(5) << iTube.PointIn << " " << setw(5) << iTube.PointOut;
		}

		prot << "\n\n Point  Node  No tubes in point ";
		for (auto& iPoint : Points) {
			prot << "\n" << setw(5) << iPoint.Number;
			if (iPoint.NbNd != -1) {
				prot << "  " << setw(5) << iPoint.NbNd;
			}
			else {
				prot << "   --- ";
			}
			prot << " " << setw(5) << iPoint.NoTb << " -> ";
			for (auto& jTb : iPoint.NbTb) {
				prot << setw(8) << jTb;
			}
		}
		prot << "\n\n Node Point";
		for (auto& iNode : Nodes) {
			prot << "\n " << setw(5) << iNode.Number << "  " << setw(5) << iNode.NbPt;
		}
		prot << endl;
	}
	/*!  -------------------------
		 * determination of branches
		  -------------------------
		  a branch is a sequence of connected tubes that connects 2 nodes */
		  /// starting from tubes connected to start node\n
		  /// if end point of this tube is node -> done\n
		  /// finding tube connected to first tube\n
		  /// adjust flow direction and check for node\n
		  /// repeat until node is found\n
	mBr = MINUS1; //=0xffffffffffffffffULL  corresponds to 0 - 1
	for (auto& iNode : Nodes) {
		iPt = iNode.NbPt; 
		for (size_t& nTb : Points[iPt].NbTb) {
			if (Tubes[nTb].NbBr == MINUS1) { // if tube already belong to a branch->skip
				if (iPt == Tubes[nTb].PointOut) {
					if (Base.showMeshDetail) {
						prot << " tube :" << Tubes[nTb].Number << " reverseDirection " << endl;
					}
					Tubes[nTb].reverseDirection();
				}
				PtTbOut = Tubes[nTb].PointOut;
				++mBr;
				Branches.push_back(_branch());
				_branch* mBranch = &Branches[mBr];
				mBranch->Number = mBr;
				mBranch->NbPtIn = Tubes[nTb].PointIn;
				mBranch->NbNdIn = iNode.Number;
				mBranch->NbTbInBr.push_back(Tubes[nTb].Number);
				mBranch->mTbInBr = 0;
				Tubes[nTb].NbBr = mBr;
				mBranch->qSum = Tubes[nTb].q;
				iNode.mBrInNd++;
				iNode.NbBr.push_back(mBr);
				iNode.mBrLeave++;
				iNode.NbBrLeave.push_back(mBr);
				if (Points[PtTbOut].NbNd != MINUS1) { //! end point of tube is at a node 
					mBranch->NbPtOut = PtTbOut;
					mBranch->NbNdOut = Points[PtTbOut].NbNd;
					if (Base.showMeshDetail) {
						prot << "\n end of branch iNd " << iNode.Number << " mBr " << mBr << " nTb " << Tubes[nTb].Number << " mTbInBr " << Branches[mBr].mTbInBr;
					}
				}
				else {
					iTb = nTb;
					while (true) {
						for (size_t& iTb1 : Points[PtTbOut].NbTb) {
							if (iTb1 != iTb) {//! there are only 2 tubes connected to the point: iTb and the other, now we have to handle the other
								if (PtTbOut == Tubes[iTb1].PointOut) {
									//! one end point is connected to endpoint of next tube, this is not correct\n
									//! reverse direction to get end point connected to starting point of next tube\n
									//! we have to deal with direction of branch later\n
									if (Base.showMeshDetail) {
										prot << "\n tube :" << Tubes[iTb1].Number << " reverseDirection ";
									}
									Tubes[iTb1].reverseDirection();
								}
								mBranch->NbTbInBr.push_back(Tubes[iTb1].Number);
								mBranch->mTbInBr++;
								if (Base.showMeshDetail) {
									prot << "\n iNd " << iNode.Number << " mBr " << mBr << " iTb1 " << Tubes[iTb1].Number << " mTbInBr " << mBranch->mTbInBr;
								}
								Tubes[iTb1].NbBr = mBr;
								mBranch->qSum += Tubes[iTb1].q;

								PtTbOut = Tubes[iTb1].PointOut;
								if (Points[PtTbOut].NbNd != MINUS1) {
									mBranch->NbPtOut = PtTbOut;
									mBranch->NbNdOut = Points[PtTbOut].NbNd;
									if (Base.showMeshDetail) {
										prot << "\n end of Branch iNd " << iNode.Number << " mBr " << mBr << " iTb1 " << Tubes[iTb1].Number << " NbNdOut " << mBranch->NbNdOut;
									}
									goto EndOfThisBranch;
								}
								iTb = iTb1;
								break;
							}
						}
					}
					cout << "\nbranch starting with tube # " << mBranch->NbTbInBr[0]
						<< " has not a node at the end";
					prot << "\nbranch starting with tube # " << mBranch->NbTbInBr[0]
						<< " has not a node at the end";
					prot.close();
					exit(1);
				}
			EndOfThisBranch:
				;
				Nodes[mBranch->NbNdOut].mBrInNd++;
				Nodes[mBranch->NbNdOut].NbBr.push_back(mBr);
				Nodes[mBranch->NbNdOut].mBrArrive++;
				Nodes[mBranch->NbNdOut].NbBrArrive.push_back(mBr);
			}
		}
	}
	Branches.shrink_to_fit();

	for (auto& iBranch : Branches) {
		iBranch.deltaH = Nodes[iBranch.NbNdOut].Elev -
			Nodes[iBranch.NbNdIn].Elev;
		iBranch.minArea = 1e30;
		for (auto jTb : iBranch.NbTbInBr) {
			iBranch.minArea = fmin(Tubes[jTb].area * Tubes[jTb].NoParallel, iBranch.minArea);
		}
	}

	StartResistance();

	initFlow();

	// -------------------------------------------------
	// determine which nodes are T-piece
	// and which branches are inline or branching off at this T-piece
	// here only geometry is checked
	// -------------------------------------------------
	// condition: straight: same diameter and angle between tubes 0 or 180 deg,
	// off: diameter <= straight, every angle
	determineTee();

	//	Print2dxf(PathFile, ShowMode::BranchesNodes);

	//------------------------------------------------
	// determine angle between tubes in a branch and if it's an U or S arrangement
	// at this place the flow directions are established
	//-------------------------------------------------
	for (const auto& iBranch : Branches) {
		//      cout<<" iBranch "<<iBranch<<endl;
		if (iBranch.mTbInBr > 0) {
			for (size_t iTbInBr = 1; iTbInBr <= iBranch.mTbInBr; iTbInBr++) {
				_tube* iTube = &Tubes[iBranch.NbTbInBr[iTbInBr]];
				_tube* iTubePrev = &Tubes[iBranch.NbTbInBr[iTbInBr - 1]];
				if (iTube->RadiusBend < 1e-3) { // straight tube
				//------------------------------------------------
				// we have to make exceptions: drum and mud drum
				// there are no bends at drum nozzle
				// it is handled by inlet/outlet resistance factor
				// as well as straight tubes after bends
				//-------------------------------------------------
				   if (iTubePrev->RadiusBend > 1e-3 || // straight tube after bend
						iTubePrev->Dia > MinDrumDiameter || // straight tube after drum
					iTube->Dia > MinDrumDiameter) continue; // tube in drum
					iTube->beta = iTubePrev->Angle3d(*iTube) * 180. / M_PI;
				// cout<<"iBr "<< iBranch.Number<<" mTbInBr "<<Branch[iBranch].mTbInBr  <<" iTbInBr "<<iTbInBr<<" iTb "<<iTb<< " beta "<< Tube[iTb].beta<<endl;
					if (iBranch.mTbInBr > 1 && iTbInBr > 1) {
						_tube* iTube2 = &Tubes[iBranch.NbTbInBr[iTbInBr - 2]];
						if (iTubePrev->Length / iTubePrev->Dia < 20. &&
							fabs(iTubePrev->Dia - iTube->Dia) < 1e-4 &&
							fabs(iTube2->Dia - iTube->Dia) < 1e-4) {
							angle1 = iTube2->Angle3d(*iTube);
							if (angle1 < 0.05) iTubePrev->UorS = USArrangement::S; // S-bend
							if (fabs(angle1 - M_PI) < 0.05) iTubePrev->UorS = USArrangement::U; //U-bend
						}
					}
				}
				else { //elbow 
//angle between iTbInBr-1 and iTbInBr+1 = iTube.beta
					if (iBranch.mTbInBr >= 4 && iTbInBr > 2 && iTbInBr + 1 <= iBranch.mTbInBr){
						_tube* iTube3 = &Tubes[iBranch.NbTbInBr[iTbInBr - 3]];
						_tube* iTube1 = &Tubes[iBranch.NbTbInBr[iTbInBr + 1]];
						if (iTubePrev->Length / iTubePrev->Dia < 20. &&
							fabs(iTubePrev->Dia - iTube->Dia) < 1e-4 &&
							fabs(iTube3->Dia - iTube->Dia) < 1e-4 &&
							fabs(iTube1->Dia - iTube->Dia) < 1e-4) {
							angle1 = iTube1->Angle3d(*iTube3);
							if (angle1 < 0.05) iTube->UorS = USArrangement::S; // S-bend
							if (fabs(angle1 - M_PI) < 0.05) iTube->UorS = USArrangement::U; //U-bend
							iTbInBr++;
						}
					}
				}
			}
		}
	}

	if (Base.showMesh) {
		prot << "\n    *****************";
		prot << "\n    * Mesh checking *";
		prot << "\n    *****************";
		prot << "\n mPt " << mPt << " mTb " << mTb << " mNd " << mNd << " mBr " << mBr;

		prot << "\n\n___________________________________";
		prot << "\n|             tubes               | ";
		prot << "\n__________________________________";
		prot << "\n|        |     Point   | belongs  |";
		prot << "\n| number | start | end | to branch|";
		prot << "\n|________|_______|_____|__________|";
		for (const auto& iTube : Tubes) {
			prot << "\n" << setw(5) << iTube.Number << "    " << setw(5) << iTube.PointIn << "  " << setw(5) << iTube.PointOut << "   " << setw(5) << iTube.NbBr;
		}

		prot << "\n\n Point  Node  No tubes in point ";
		for (const auto& iPoint : Points) {
			prot << "\n" << setw(5) << iPoint.Number;
			if (iPoint.NbNd != MINUS1) {
				prot << "  " << setw(5) << iPoint.NbNd;
			}
			else {
				prot << "   --- ";
			}
			prot << " " << setw(5) << iPoint.NoTb << " -> ";
			size_t index = 0;
			size_t end = iPoint.NbTb.size();
			for (; index < end; index++) {
				if (index > 0 && index % 10 == 0) {
					prot << "\n                      ";
				}
				prot << setw(8) << iPoint.NbTb[index];
			}
		}
		prot << "\n\n Node Point";
		for (const auto& iNode : Nodes) {
		   prot << "\n" << setw(5) << iNode.Number << "  " << setw(5) << iNode.NbPt;
		}
		prot << endl;
		prot << "\n -----------------------------------------------------------" << endl;
		prot << endl;
			//prot << "\nbranch NbPtIn NbPtOut NbNdIn NbNdOut     qSum   tubes  ";
			//for (auto& iBranch : Branches) {
			//	prot << "\n " << setw(5) << iBranch.Number << "  " << setw(5) << iBranch.NbPtIn << "   " << setw(5) <<
			//		iBranch.NbPtOut << " " << setw(5) << iBranch.NbNdIn << "   " << setw(5) <<
			//		iBranch.NbNdOut << " " << setw(8) << iBranch.qSum << " ";
			//	for (auto& element : iBranch.NbTbInBr)
			//		prot << setw(5) << element << " ";
			//}
		if (Base.showMeshDetail) {
			prot << "\nbranch pseudo resistance  ";
			for (auto& iBranch : Branches) {
				prot << "\n " << setw(5) << iBranch.Number << "  " << setw(8) << iBranch.dPConstant;
			}
		}
		prot << endl;
		prot << "\n\n   Node  Point elevation     branches in node";
		for (auto& iNode : Nodes) {
			prot << "\n" << setw(5) << iNode.Number << " " << setw(6) << iNode.NbPt << " " << setw(8) << iNode.Elev << " : ";
			size_t index = 0;
			for (; index <= iNode.mBrInNd; index++) {
				if (index > 0 && index%10 == 0) {
					prot << "\n                        ";
				}
				prot << setw(8) << iNode.NbBr[index];
			}
			prot << "\n " << iNode.mBrArrive + 1 << " arrive :";
			for (index = 0; index <= iNode.mBrArrive; index++) {
				if (index > 0 && index % 10 == 0) {
					prot << "\n              ";
				}
				prot << setw(8) << iNode.NbBrArrive[index];
			}
//			for (auto& element : iNode.NbBrArrive)
//				prot << " " << setw(5) << element;
			prot << "\n " << iNode.mBrLeave + 1 << " leave :";
			for (index = 0; index <= iNode.mBrLeave; index++) {
				if (index > 0 && index % 10 == 0) {
					prot << "\n            ";
				}
				prot << setw(8) << iNode.NbBr[index];
			}
//			for (auto& element : iNode.NbBrLeave)
//				prot << " " << setw(5) << element;
		}
		prot << endl;
		prot << "\nbranch NbPtIn NbPtOut NbNdIn NbNdOut     qSum     g    kind  tubes  ";
		for (const auto& iBranch : Branches) {
			prot << "\n " << setw(5) << iBranch.Number << "  " << setw(5) << iBranch.NbPtIn
				<< "   " << setw(5) << iBranch.NbPtOut << " " << setw(5) << iBranch.NbNdIn
				<< "   " << setw(5) << iBranch.NbNdOut << " " << setw(8) << iBranch.qSum << " "
				<< setw(8) << iBranch.g;
			switch (iBranch.Kind)
			{
			case KindOfBranch::Undefined:
				prot << " undefined        :";
				break;
			case KindOfBranch::Heated:
				prot << " heated           :";
				break;
			case KindOfBranch::Downcomer:
				prot << " downcomer        :";
				break;
			case KindOfBranch::Down2Heated:
				prot << " downcomer->heated:";
				break;
			case KindOfBranch::Heated2Drum:
				prot << " heated->drum     :";
				break;
			case KindOfBranch::Drum2Down:
				prot << " drum->downcomer  :";
				break;
			default:
				prot << " ---------------  :";
				break;
			}
			for (const size_t& element : iBranch.NbTbInBr)
				prot << setw(5) << element << " ";
		}
		prot << endl;
		prot << "\n\n Node is T-piece   straight    off  orientation";
		for (const auto& iNode : Nodes) {
			if (iNode.Number > 0 && iNode.IsT) {
				prot << "\n" << setw(5) << iNode.Number << "     " << iNode.IsT << "    ";
				for (j = 0; j <= 1; ++j) {
					prot << " " << setw(5) << iNode.NbBrTStraight[j];
				}
				prot << "  " << setw(5) << iNode.NbBrTOff << "        ";

				switch (iNode.TOrientation)
				{
				case TeeOrientation::NoTee:
					prot << "NoTee";
					break;
				case TeeOrientation::Undetermined:          // no distinct orientation
					prot << "Undetermined";
					break;
				case	TeeOrientation::StraightHorOffHor:     // 1 Straight horizontal, Off horizontal
					prot << "StraightHorOffHor";     // 1 Straight horizontal, Off horizontal
					break;
				case		TeeOrientation::StraightHorOffVerUp:   // 2 Straight horizontal, Off vertical up
					prot << "StraightHorOffVerUp";
					break;
				case		TeeOrientation::StraightHorOffVerDown: // 3 Straight horizontal, Off vertical down
					prot << "StraightHorOffVerDown";
					break;
				case		TeeOrientation::StraightVerOffHor:     //                             4 : Straight vertical, Off horizontal
					prot << "StraightVerOffHor";
					break;
				case		TeeOrientation::OffHor:                //                             5 : Straight every else, Off horizontal
					prot << "OffHor";
					break;
				case		TeeOrientation::OffVerUp:              //                             6 : Straight every else, Off vertical up
					prot << "OffVerUp";
					break;
				case		TeeOrientation::OffVerDown:             //                             7 : Straight every else, Off vertical down
					prot << "OffVerDown";
					break;
				default:
					break;
				}
			}
		}
		prot << endl;
	}
	Print2dxf(ShowMode::BranchesNodes);
	Base.iterg = 0;
	Print2dxf(ShowMode::Arrows);
	Print2dxf(ShowMode::Flow);
	//  exit(1);
	return 0;
} /* mesh */

/**
 * \fn findMaxNodesConnected()
 * \brief find the maximum number of nodes connected to one (not particular) node (needed for size allocation of sparse matrix)
 *
 * Although the maximum number of branches in a node can be easily determined in mesh-function it is a big difference\n
 * in bidrum boilers where we have ~ 1000 branches connecting steam and mud drum but only one other node connected
 * \return the number (eigen requires int)
 */
int findMaxNodesConnected() {
	size_t maxNodesConnected = 3;
	stVector NodesConnected; ///< vector to hold the numbers of the nodes that are connected to the node 
	for (auto& iNode : Nodes) {
		iNode.NbBr.shrink_to_fit();
		/// as steam drum is not regarded in the size of equation system it will be omitted here also
		if (iNode.Number > DRUM) {
			for (size_t& iBr : iNode.NbBr) {
				size_t otherNd = Branches[iBr].NbNdIn;
				if (iNode.Number == otherNd)  otherNd = Branches[iBr].NbNdOut;
				bool found = false;
				/// look-up if otherNd is already touched before. If yes, break and set flag
				for (size_t& iNC : NodesConnected) {
					if (iNC == otherNd) {
						found = true;
						break;
					}
				}
				/// if the node number is not available in NodesConnected, add it
				if (!found) {
					NodesConnected.push_back(otherNd);
				}
			}
			/// the number of connected nodes is the size of NodesConnected\n 
			/// if the size greater than the already max number -> update
			if (NodesConnected.size() > maxNodesConnected) {
				maxNodesConnected = NodesConnected.size();
			}
			NodesConnected.clear();
		}
	}
	//cout << "\n max Nodes Connected " << maxNodesConnected;
	return static_cast<int> (maxNodesConnected);
}
