/*****************************************************************//**
 * \file   step.cpp
 * \brief  determines the flow and direction for next iteration step
 * 
 * Determination of new flow in branches and reversal of flow\n 
 * the new flow is not the flow as calculated in "CalcPressureNodes" \n
 * it is a damped step except there is a reversal of flow direction 
 * \return bool true if at least one branch changed direction
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/

//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

bool Step(void) {
  /* Local variables */
  stVector reversal;  // branch number of candidate for reversal of direction
  bool IsReverse;
  int istep = 3;
  size_t i;
  double damping = 0.2, 
//     gSumIn, gSumOut, 
     diff1, diff2;
//  bool revDone;

  IsReverse = false;

  if (Base.iterg == 1) {
     // in the first flow iteration the flow direction was guessed
     // The flow in all nodes is now balanced
     // use the new flow directly without any damping (basically, as a new
     // starting value) but not if "vertical" or heated branches should be reversed
     bool ReversalOfHeatedBranches = false;
     for (auto& iBranch : Branches) {
        if (iBranch.gNew < 0. ) {
           if(iBranch.qSum > 1.e-3 || fabs(iBranch.deltaH) > 0.1) {
              ReversalOfHeatedBranches = true;
              break;
           }
        }
     }
     if (!ReversalOfHeatedBranches) {
        for (auto& iBranch : Branches) {
           if (iBranch.gNew < 0.) {
              // guess was wrong -> immediate change of direction
              iBranch.reverseDirection();
//              prot << "\n branch " << iBranch.Number << " reversed in " << iBranch.NbNdIn << " out " << iBranch.NbNdOut;
              IsReverse = false;
              iBranch.g = -iBranch.gNew;
           }
           else {
              iBranch.g = iBranch.gNew;
           }
        }
        for (_node element : Nodes) {
           if (element.NbBrArrive.size() == 0 || element.NbBrLeave.size() == 0) {
              prot << "\n in step: iterg " << Base.iterg << " node " << element.Number << " no arriving or leaving branch";
              cout << "\n in step: iterg " << Base.iterg << " node " << element.Number << " not arriving or leaving branch";
           }
        }
        return IsReverse;
     }
  }

  for (auto& iBranch : Branches) {
    if (iBranch.gNew < 0.) {
      iBranch.neg += 1;
    }
    if (iBranch.neg >= istep) {  // &&
      //			(iBranch.Kind != KindOfBranch::Heated ||
      // iBranch.isHorizontalHeated) && 			iBranch.gNew <
      // 0.)) { //reversal only if negative flow in this step
      IsReverse = true;
      //			reversal.push_back(iBranch.Number);
      //			rev++;
    }
  }
  // to keep mass balance in nodes: if one branch is reversed, use new flow(that
  // is balanced in all nodes) in all branches
  if (IsReverse) {
    for (auto& iBranch : Branches) {
      // using the calculated values from mass balance -> no problems with
      // enthalpy calculation
      iBranch.gPrev3 = iBranch.gPrev2;
      iBranch.gPrev2 = iBranch.gPrev1;
      iBranch.gPrev1 = iBranch.g;
      if (iBranch.gNew < 0.) {
        iBranch.NoChanges++;
        if (iBranch.NoChanges > 5) {
          _node& StartNode = Nodes[iBranch.NbNdIn];
          _node& EndNode = Nodes[iBranch.NbNdOut];
          if (Base.showReverseDetail) {
             prot << "\nbefore reversal branch " << iBranch.Number;
             prot << "\n StartNode number " << StartNode.Number << "  mBrInNd "
                << StartNode.mBrInNd << "  sizeBrInNd "
                << StartNode.NbBr.size() << "  IsT " << StartNode.IsT
                << "  mBrLeave " << StartNode.mBrLeave << "  sizeBrLeave "
                << StartNode.NbBrLeave.size();
             prot << "\n Branch leave ";
             for (i = 0; i <= StartNode.mBrLeave; i++) {
                prot << StartNode.NbBrLeave[i] << " ";
             }

             prot << "\n EndNode number " << EndNode.Number << "  mBrInNd "
                << EndNode.mBrInNd << "  sizeBrInNd " << EndNode.NbBr.size()
                << "  IsT " << EndNode.IsT << "  mBrArrive "
                << EndNode.mBrArrive << "  sizeBrArrive "
                << EndNode.NbBrArrive.size();
             prot << "\n Branch arrive ";
             for (i = 0; i <= EndNode.mBrArrive; i++) {
                prot << EndNode.NbBrArrive[i] << " ";
             }
          } // REVERSE_DEBUG
          if (StartNode.mBrLeave == 0 || EndNode.mBrArrive == 0) {
            iBranch.reverseDirection();
            StartNode = Nodes[iBranch.NbNdIn];
            EndNode = Nodes[iBranch.NbNdOut];
            if (Base.showReverseDetail) {
               prot << "\nafter reversal branch " << iBranch.Number;
               prot << "\n StartNode number " << StartNode.Number << "  mBrInNd "
                  << StartNode.mBrInNd << "  sizeBrInNd "
                  << StartNode.NbBr.size() << "  IsT " << StartNode.IsT
                  << "  mBrLeave " << StartNode.mBrLeave << "  sizeBrLeave "
                  << StartNode.NbBrLeave.size();
               prot << "\n Branch leave ";
               for (i = 0; i <= StartNode.mBrLeave; i++) {
                  prot << StartNode.NbBrLeave[i] << " ";
               }

               prot << "\n EndNode number " << EndNode.Number << "  mBrInNd "
                  << EndNode.mBrInNd << "  sizeBrInNd " << EndNode.NbBr.size()
                  << "  IsT " << EndNode.IsT << "  mBrArrive "
                  << EndNode.mBrArrive << "  sizeBrArrive "
                  << EndNode.NbBrArrive.size();
               prot << "\n Branch arrive ";
               for (i = 0; i <= EndNode.mBrArrive; i++) {
                  prot << EndNode.NbBrArrive[i] << " ";
               }
            } // REVERSE_DEBUG
          }
          if (StartNode.mBrLeave > 0 && EndNode.mBrArrive > 0) {
            cout << "\n branch " << iBranch.Number << " set to zero ";
            prot << "\n branch " << iBranch.Number << " set to zero ";
            iBranch.isFlowSet2zero = true;
            iBranch.g = 0.;
            iBranch.dPdyn = 0.;
            iBranch.dPstat = 0.;
            iBranch.dPIn = 0.;
            iBranch.dPOut = 0.;
            iBranch.dPLinear = 0.;
            iBranch.dPConstant = 0.;
            size_t iStartPoint = StartNode.NbPt;
            _point& StartPoint = Points[iStartPoint];
            size_t iEndPoint = EndNode.NbPt;
            _point& EndPoint = Points[iEndPoint];

            if (Base.showReverseDetail) {
               prot << "\nbefore ";
               prot << "\n StartNode number " << StartNode.Number << "  mBrInNd "
                  << StartNode.mBrInNd << "  sizeBrInNd "
                  << StartNode.NbBr.size() << "  IsT " << StartNode.IsT
                  << "  mBrLeave " << StartNode.mBrLeave << "  sizeBrLeave "
                  << StartNode.NbBrLeave.size() << "  iStartPoint "
                  << iStartPoint << "  StartPoint.NoTb " << StartPoint.NoTb
                  << "  sizeStartPoint.NoTb " << StartPoint.NbTb.size();
               prot << "\n Branch leave ";
               for (i = 0; i <= StartNode.mBrLeave; i++) {
                  prot << StartNode.NbBrLeave[i] << " ";
               }
               prot << "\n tube in point ";
               for (i = 0; i < StartPoint.NoTb; i++) {
                  prot << StartPoint.NbTb[i] << " ";
               }

               prot << "\n EndNode number " << EndNode.Number << "  mBrInNd "
                  << EndNode.mBrInNd << "  sizeBrInNd " << EndNode.NbBr.size()
                  << "  IsT " << EndNode.IsT << "  mBrArrive "
                  << EndNode.mBrArrive << "  sizeBrArrive "
                  << EndNode.NbBrArrive.size() << "  iEndPoint " << iEndPoint
                  << "  EndPoint.NoTb " << EndPoint.NoTb
                  << "  sizeEndPoint.NoTb " << EndPoint.NbTb.size();
               prot << "\n Branch arrive ";
               for (i = 0; i <= EndNode.mBrArrive; i++) {
                  prot << EndNode.NbBrArrive[i] << " ";
               }
               prot << "\n tube in point ";
               for (i = 0; i < EndPoint.NoTb; i++) {
                  prot << EndPoint.NbTb[i] << " ";
               }
            } // REVERSE_DEBUG

            if (StartNode.IsT)
              StartNode.IsT = false;

            for (i = 0; i <= StartNode.mBrInNd; i++) {
               if (StartNode.NbBr[i] == iBranch.Number) {
                  StartNode.NbBr.erase(StartNode.NbBr.begin() + i);
                  break;
               }
            }
            StartNode.mBrInNd--;

            for (i = 0; i <= StartNode.mBrLeave;        i++) {
              if (StartNode.NbBrLeave[i] == iBranch.Number) {
                StartNode.NbBrLeave.erase(StartNode.NbBrLeave.begin() +i);
                break;
              }
            }
            StartNode.mBrLeave--;
            for (i = 0; i < StartPoint.NoTb; i++) {
              size_t iTube = StartPoint.NbTb[i];
              if (Branches[Tubes[iTube].NbBr].Number == iBranch.Number) {
                StartPoint.NbTb.erase(StartPoint.NbTb.begin() + i);
                break;
              }
            }
            StartPoint.NoTb--;
            switch (StartNode.mBrInNd) {
              case 2: {
                 _tube Tube0 = Tubes[StartPoint.NbTb[0]];
                 _tube Tube1 = Tubes[StartPoint.NbTb[1]];
                 _tube Tube2 = Tubes[StartPoint.NbTb[2]];

                if (StartNode.isTee(Tube0, Tube1, Tube2))
                  ;
                else if (StartNode.isTee(Tube1, Tube2, Tube0))
                  ;
                else
                  (StartNode.isTee(Tube0, Tube2, Tube1));
              } break;
              case 1:
                // we have 2 branches connected to the node
                // basically those 2 branches can be connected to one
                // the nodal pressure calculation also works with 2 branches
                // only a deflection need to be taken into account
                _tube TubeArrive = Tubes[Branches[StartNode.NbBrArrive[0]].NbTbInBr[Branches[StartNode.NbBrArrive[0]].mTbInBr]];
                _tube TubeLeave = Tubes[Branches[StartNode.NbBrLeave[0]].NbTbInBr[0]];
                TubeLeave.beta = TubeArrive.Angle3d(TubeLeave) * 180. / M_PI;
                TubeLeave.RadiusBend = 0.;
                break;
            }

            if (EndNode.IsT)
              EndNode.IsT = false;
            for (i = 0; i <= EndNode.mBrInNd;             i++) {
              if (EndNode.NbBr[i] == iBranch.Number) {
                EndNode.NbBr.erase(EndNode.NbBr.begin() + i);
                break;
              }
            }
            EndNode.mBrInNd--;

            for (i = 0; i <= EndNode.mBrArrive; i++) {
              if (EndNode.NbBrArrive[i] == iBranch.Number) {
                EndNode.NbBrArrive.erase(EndNode.NbBrArrive.begin() +i);
                break;
              }
            }
            EndNode.mBrArrive--;
            for (i = 0; i < EndPoint.NoTb; i++) {
              size_t iTube = EndPoint.NbTb[i];
//              if (Branches[Tubes[iTube].NbBr].Number == iBranch.Number) {
                 if (Tubes[iTube].NbBr == iBranch.Number) {
                    EndPoint.NbTb.erase(EndPoint.NbTb.begin() + i);
                break;
              }
            }
            EndPoint.NoTb--;
            switch (EndNode.mBrInNd) {
            case 2: {
               _tube Tube0 = Tubes[EndPoint.NbTb[0]];
               _tube Tube1 = Tubes[EndPoint.NbTb[1]];
               _tube Tube2 = Tubes[EndPoint.NbTb[2]];

               if (EndNode.isTee(Tube0, Tube1, Tube2));
               else if (EndNode.isTee(Tube1, Tube2, Tube0));
               else (EndNode.isTee(Tube0, Tube2, Tube1));
               } 
                  break;
              case 1:
                // we have 2 branches connected to the node
                // basically those 2 branches can be connected to one
                // the nodal pressure calculation also works with 2 branches
                // only a deflection need to be taken into account
                _tube TubeArrive =  Tubes[Branches[EndNode.NbBrArrive[0]].NbTbInBr[Branches[EndNode.NbBrArrive[0]].mTbInBr]];
                _tube TubeLeave = Tubes[Branches[EndNode.NbBrLeave[0]].NbTbInBr[0]];
                TubeLeave.beta = TubeArrive.Angle3d(TubeLeave) * 180. / M_PI;
                TubeLeave.RadiusBend = 0.;
                break;
            }
            if (Base.showReverseDetail) {
               prot << "\nafter ";
               prot << "\n StartNode number " << StartNode.Number << "  mBrInNd "
                  << StartNode.mBrInNd << "  sizeBrInNd "
                  << StartNode.NbBr.size() << "  IsT " << StartNode.IsT
                  << "  mBrLeave " << StartNode.mBrLeave << "  sizeBrLeave "
                  << StartNode.NbBrLeave.size() << "  iStartPoint "
                  << iStartPoint << "  StartPoint.NoTb " << StartPoint.NoTb
                  << "  sizeStartPoint.NoTb " << StartPoint.NbTb.size();
               prot << "\n Branch leave ";
               for (i = 0; i <= StartNode.mBrLeave; i++) {
                  prot << StartNode.NbBrLeave[i] << " ";
               }
               prot << "\n tube in point ";
               for (i = 0; i < StartPoint.NoTb; i++) {
                  prot << StartPoint.NbTb[i] << " ";
               }

               prot << "\n EndNode number " << EndNode.Number << "  mBrInNd "
                  << EndNode.mBrInNd << "  sizeBrInNd " << EndNode.NbBr.size()
                  << "  IsT " << EndNode.IsT << "  mBrArrive " << EndNode.mBrArrive
                  << "  sizeBrArrive " << EndNode.NbBrArrive.size()
                  << "  iEndPoint " << iEndPoint << "  EndPoint.NoTb "
                  << EndPoint.NoTb << "  sizeEndPoint.NoTb "
                  << EndPoint.NbTb.size();
               prot << "\n Branch arrive ";
               for (i = 0; i <= EndNode.mBrArrive; i++) {
                  prot << EndNode.NbBrArrive[i] << " ";
               }
               prot << "\n tube in point ";
               for (i = 0; i < EndPoint.NoTb; i++) {
                  prot << EndPoint.NbTb[i] << " ";
               }
            }// REVERSE_DEBUG

          } else {
            prot << "\n after removing one branch (set to zero) one node has "
                    "no arriving or leaving branch ";
            prot << "\n even with a change of direction";
            prot << "\n Something is rotten in the state of Denmark "
                    "(Hamlet(1.4), Marcellus to Horatio)";
          }
        } else {
        iBranch.reverseDirection();
        if (Base.showReverseDetail) {
           prot << "\n branch " << iBranch.Number << " reversed in " << iBranch.NbNdIn << " out " << iBranch.NbNdOut;
        }
        iBranch.g = -iBranch.gNew;
          iBranch.neg = 0;
        }
      } else {
        iBranch.g = iBranch.gNew;
      }
    }
  } else {  // no reversal: use damping for the new flow
    for (auto& iBranch : Branches) {
      diff1 = fabs((iBranch.g - iBranch.gPrev1) / iBranch.g) * 100.;  // in %
      diff2 = fabs((iBranch.g - iBranch.gPrev2) / iBranch.g) * 100.;  // in %

      iBranch.gPrev3 = iBranch.gPrev2;
      iBranch.gPrev2 = iBranch.gPrev1;
      iBranch.gPrev1 = iBranch.g;
      if (iBranch.gNew > 0.) {
        iBranch.g += (iBranch.gNew - iBranch.g) * damping;
      } else {  // gNew negative, we go towards 0
        //						 if
        //(iBranch.Kind
        //==
        // KindOfBranch::Heated) {
        // iBranch.g = (1. - damping) * iBranch.g;
        //						 }
        //						 else {
        iBranch.g = (1. - damping) * iBranch.g;
        //						 }
      }
      if (fabs(iBranch.g) < 1e-5) {
        iBranch.g = 0.;
      }
    }
  }

  for (_node element : Nodes) {
     if (element.NbBrArrive.size() == 0 || element.NbBrLeave.size() == 0) {
        prot << "\n in step: iterg " << Base.iterg << " node " << element.Number << " no arriving or leaving branch";
        cout << "\n in step: iterg " << Base.iterg << " node " << element.Number << " not arriving or leaving branch";
     }
  }


#if 0
	if (Base.showReverse) {
		if (reversal.size() > 0) {
			prot << "\n reversal ";
			for (auto& iBr : reversal) {
				prot << "  " << iBr;
			}
			prot << "\n end reversal \n";
		}
	}
	// check, if reversal of flow will result in nodes with
	// no arriving or leaving flow
	// do it rev times if only one branch is reversed
	// and this allows another branch to be reversed
	for (i = 1; i <= rev; i++) {
		revDone = false;
		for (auto it = reversal.begin(); it != reversal.end();) {
			size_t iBr = *it;
			if (Base.showReverse) {
				prot << "\n candidate " << iBr;
				prot << "\n NodeIn leave " << Nodes[Branches[iBr].NbNdIn].mBrLeave
					<< " NodeOut arrive " << Nodes[Branches[iBr].NbNdOut].mBrArrive;
			}
			if (Nodes[Branches[iBr].NbNdIn].mBrLeave == 0) {
				//this Branch was only Branch to leave
				if (Base.showReverse) {
					prot << "\n only branch to leave ";
				}
				++it;
			}
			else if (Nodes[Branches[iBr].NbNdOut].mBrArrive == 0) {
				//this Branch was the only one to arrive
				if (Base.showReverse) {
					prot << "\n only branch to arrive ";
				}
				++it;
			}
			else { //reverse
				if (Base.showReverse) {
					prot << "\n reverse branch ";
				}
				Branches[iBr].neg = 0;
				Branches[iBr].NoChanges++;
				///the flows in the other branches are already set
				///if using the calculated gNew the balance in node is far off -> leading to problems in CalcEnthalpyNodes
				size_t iNdIn = Branches[iBr].NbNdIn;
				gSumIn = 0.;
				for (auto& iBr2 : Nodes[iNdIn].NbBrArrive) {
					if (iBr2 != iBr) gSumIn += Branches[iBr2].g;
				}
				for (auto& iBr2 : Nodes[iNdIn].NbBrLeave) {
					if (iBr2 != iBr) gSumIn -= Branches[iBr2].g;
				}
				size_t iNdOut = Branches[iBr].NbNdOut;
				gSumOut = 0.;
				for (auto& iBr2 : Nodes[iNdOut].NbBrArrive) {
					if (iBr2 != iBr) gSumOut += Branches[iBr2].g;
				}
				for (auto& iBr2 : Nodes[iNdOut].NbBrLeave) {
					if (iBr2 != iBr) gSumOut -= Branches[iBr2].g;
				}
				Branches[iBr].g = (fabs(gSumIn) + fabs(gSumOut)) / 2.;

				if (Base.showReverse) {
					prot << "Branch " << iBr << " gSumIn " << gSumIn << " gSumOut " << gSumOut << " gnew " << -Branches[iBr].gNew << " g " << Branches[iBr].g << endl;
				}

				/*           --------------- */
				/*          Flow reversal    */
				/*           --------------- */
				//				prot << "\n start branch reverseDirection :" << iBr << endl;
				//				cout << "\n start branch reverseDirection :" << iBr << endl;

				Branches[iBr].reverseDirection();
				// delete the previous flows, after reverse it is misleading
				Branches[iBr].gPrev1 = 0.;
				Branches[iBr].gPrev2 = 0.;
				Branches[iBr].gPrev3 = 0.;
				it = reversal.erase(it);
				revDone = true;
				IsReverse = true;
			}
		}
		if (Base.showReverse) {
			prot << "\n after one round of reversal \n";
			for (auto& iBr : reversal) {
				prot << "  " << iBr;
			}
			prot << "\nreversal end \n";
		}
		if (reversal.empty() || !revDone) break;
	}

	reversal.clear();
#endif

  return IsReverse;
} /* step */
