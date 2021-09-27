//#include "stdafx.h"

#undef MAINFUNCTION
#include "CommonHeader.h"
/**
 * \fn _branch::reverseDirection()
 *
 */
int _branch::reverseDirection() {
	/* Local variables */
//	size_t i, j;
	/*           --------------- */
	/*          Flow reversal    */
	/*           --------------- */
	if (Base.showReverse) {
		prot << "\n start branch reverseDirection :" << Number;
	}
	swap(NbPtIn, NbPtOut);
	swap(NbNdIn, NbNdOut);
	deltaH = -deltaH;
	if (Base.showReverse) {
		prot << "\nstart branch reverseSequence mtbInBr " << mTbInBr << " NbTbInBr :";
		for (const auto& element : NbTbInBr)
			prot << " " << setw(5) << element;
	}
	/// reversing the sequence of tubes
	reverse(NbTbInBr.begin(), NbTbInBr.end());
	//i = 0;
	//j = mTbInBr;
	//while (i < j) {
	//	swap(NbTbInBr[j], NbTbInBr[i]);
	//	--j;
	//	++i;
	//}
	if (Base.showReverse) {
		prot << "\nend branch reverseSequence NbTbInBr :";
		for (auto& element : NbTbInBr)
			prot << " " << setw(5) << element;
	}
	///reverse all tubes in branch
	for (const auto& NbTb : NbTbInBr) {
		if (Base.showReverse) {
			prot << "\n start tube reverse direction NbTb : " << NbTb << endl;
		}
		::Tubes[NbTb].reverseDirection();
		if (Base.showReverse) {
			prot << "\n  iBr  NbTb Startpoint Endpoint BranchNbPtIn BranchNbPtOut";
			prot << "\n" << setw(5) << Number << " " << " " << setw(5)
				<< NbTb << " " << setw(5) << ::Tubes[NbTb].PointIn << " "
				<< setw(5) << ::Tubes[NbTb].PointOut << " " << setw(5)
				<< NbPtIn << " " << setw(5) << NbPtOut;
		}
	}
	/// Updating NbBrArrive and NbBrLeave for node at start and end
	_node& iNewNodeIn = Nodes[NbNdIn];

	if (Base.showReverse) {
		if (!iNewNodeIn.NbBrLeave.empty()) {
			prot << "\nreverse direction before add nodeIn" << iNewNodeIn.Number;
			prot << "NbBrLeave : ";
			for (auto& element : iNewNodeIn.NbBrLeave) {
				prot << element << " ";
			}
		}
	}
	iNewNodeIn.NbBrLeave.push_back(Number); //extend vector
	iNewNodeIn.NbBrLeave.shrink_to_fit();
	iNewNodeIn.mBrLeave++;
	if (Base.showReverse) {
		prot << "\nreverse direction after add NodeIn" << iNewNodeIn.Number;
		prot << "NbBrLeave : ";
		for (auto& element : iNewNodeIn.NbBrLeave) {
			prot << element << " ";
		}

		prot << "\nbefore erase nodeIn" << iNewNodeIn.Number;
		if (iNewNodeIn.NbBrArrive.empty()) {
			prot << "\nNbBrArrive is empty ";
		}
		else {
			prot << "\nNbBrArrive : ";
			for (auto& element : iNewNodeIn.NbBrArrive) {
				prot << element << " ";
			}
		}
	}
	if (!iNewNodeIn.NbBrArrive.empty()) {
		for (auto iter = iNewNodeIn.NbBrArrive.begin(); iter != iNewNodeIn.NbBrArrive.end(); iter++) {
			if (*iter == Number) {
				iNewNodeIn.NbBrArrive.erase(iter);
				iNewNodeIn.mBrArrive--;
				break;
			}
		}
		if (Base.showReverse) {
			prot << "\nafter erase nodeIn" << iNewNodeIn.Number;
			prot << "\nNbBrArrive : ";
			for (auto& element : iNewNodeIn.NbBrArrive) {
				prot << element << " ";
			}
		}
	}
	_node& iNewNodeOut = Nodes[NbNdOut];

	if (Base.showReverse) {
		prot << "\nbefore add nodeOut" << iNewNodeOut.Number;
		if (iNewNodeOut.NbBrArrive.empty()) {
			prot << "\nNbBrArrive is empty";
		}
		else {
			prot << "\nNbBrArrive : ";
			for (auto& element : iNewNodeOut.NbBrArrive) {
				prot << element << " ";
			}
		}
	}
	iNewNodeOut.NbBrArrive.push_back(Number); //extend vector
	iNewNodeOut.mBrArrive++;
	iNewNodeOut.NbBrArrive.shrink_to_fit();
	if (Base.showReverse) {
		prot << "\nafter add nodeOut" << iNewNodeOut.Number;
		prot << "\nNbBrArrive : ";
		for (auto& element : iNewNodeOut.NbBrArrive) {
			prot << element << " ";
		}

		prot << "\nbefore erase nodeOut" << iNewNodeOut.Number;
		if (iNewNodeOut.NbBrLeave.empty()) {
			prot << "\n NbBrLeave is empty ";
		}
		else {
			prot << "\n NbBrLeave : ";
			for (auto& element : iNewNodeOut.NbBrLeave) {
				prot << element << " ";
			}
		}
	}
	if (!iNewNodeOut.NbBrLeave.empty()) {
		for (auto iter = iNewNodeOut.NbBrLeave.begin(); iter != iNewNodeOut.NbBrLeave.end(); iter++) {
			if (*iter == Number) {
				iNewNodeOut.NbBrLeave.erase(iter);
				iNewNodeOut.mBrLeave--;
				break;
			}
		}
		if (Base.showReverse) {
			prot << "\nafter erase nodeOut" << iNewNodeOut.Number;
			prot << "\nNbBrLeave : ";
			for (auto& element : iNewNodeOut.NbBrLeave) {
				prot << element << " ";
			}
		}
	}
	if (Base.showReverse) {
		prot << "\n arrive iNdIn  " << setw(5) << NbNdIn << " : ";
		for (auto& element : iNewNodeIn.NbBrArrive)
			prot << " " << setw(5) << element;
		prot << "\n leave iNdIn   " << setw(5) << NbNdIn << " : ";
		for (auto& element : iNewNodeIn.NbBrLeave)
			prot << " " << setw(5) << element;
		prot << "\n arrive iNdOut " << setw(5) << NbNdOut << " : ";
		for (auto& element : iNewNodeOut.NbBrArrive)
			prot << " " << setw(5) << element;
		prot << "\n leave iNdOut  " << setw(5) << NbNdOut << " : ";
		for (auto& element : iNewNodeOut.NbBrLeave)
			prot << " " << setw(5) << element;
		prot << "\n  iBr \t j \t NbTb st Startpoint \t Endpoint \t NbPtIn \t NbPtOut" << endl;
		prot << "\n end branch reverse" << endl;
	}
	return 0;
} /* reverse Branch */

/** \fn _tube::reverseDirection
* @brief reverses flow direction in tube
*
*/
void _tube::reverseDirection() {
	double ddummy;
	swap(PointIn, PointOut);
	FactHeat = 2. - FactHeat;
	Height = -Height;
	ddummy = DiaOrificeIn;
	DiaOrificeIn = DiaOrificeOut;
	DiaOrificeOut = ddummy;
	return;
}
