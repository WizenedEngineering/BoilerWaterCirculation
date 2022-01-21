/*! \file   main.cpp */
// starting date August 1, 2016, 8:50 PM


#include "stdafx.h"
#include <cstdlib>
#define MAINFUNCTION
#include "CommonHeader.h"
#undef MAINFUNCTION

using namespace std;
extern int findMaxNodesConnected();

/*!
 * \brief calculates the water/steam flow in a closed network of tubes in natural circulation
 *
 * \param argc number of arguments
 * \param argv program arguments, the name of the data project
 * \return int
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/
 /*!
 * ---------
 *
 *    Copyright 2021 Rainer Jordan
 *
 *    Licensed under the European Union Public Licence (EUPL), Version 1.2 or - as soon they
 *    will be approved by the European Commission - subsequent
 *    versions of the EUPL (the "Licence");\n
 *    You may not use this work except in compliance with the Licence.
 *
 *    You may obtain a copy of the Licence at:
 *
 *    https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12 (and chose your language)
 *
 *    Unless required by applicable law or agreed to in
 *    writing, software distributed under the Licence is
 *    distributed on an "AS IS" basis,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *    express or implied.
 *    See the Licence for the specific language governing
 *    permissions and limitations under the Licence.
 *
 * -----------------
 */

int main(int argc, char** argv) {
	/* Local variables */
	string protname, Project, readG, iname, rname, cname, line, TubeName;
	double diff, dissipation;
	bool maxg = false;
	double diffsum;
	char answer;

	ifstream gReadData;
	ofstream dchange;
	//      istringstream iss;
	/* ---------- */
	/* data input */
	/* ---------- */
	if (argc > 1) {
		Project = argv[1];
	}
	else {
		cout << "\nProject name : ";
		cin >> Project;
	}
//	PathFile = "..\\..\\WSCDATA\\" + Project + "\\" + Project;
	PathFile = "WSCDATA\\" + Project + "\\" + Project;
	protname = PathFile + ".pro";
	cout << protname;
	// Set exceptions to be thrown on failure
	prot.exceptions(ofstream::failbit | ofstream::badbit);

	try {
		prot.open(protname.c_str());
	}
	catch (system_error& e) {
		cerr << protname << " : " << e.code().message() << endl;
		exit(1);
	}

	/**
	 * This is the main program to call the different functions.
	 *
	 * 1) input data is read from file\n
	* function: readData() */
	/* ------------------- */
	readData();

	/**
	* 2) Setting up Mesh (network) and the initial flow in Branches\n
	* function: Mesh()   */
	/* ------------------------------------------------- */
	Mesh();
	//size the matrix and vectors for equation system, mNd is known now
	//allocated here once,
	SparseMatrix <double> S(mNd, mNd); /* System sparse matrix */
	S.reserve(VectorXi::Constant(mNd, findMaxNodesConnected() + 2));
	VectorXd B(mNd); /* right hand side vector............*/
	VectorXd X(mNd); /* solution vector ...................*/

  /* --------------------------------- */
  /** 3) although the initial flow in the branches is set, it can be overwritten by data from file "*projectname*StartFlow.dat" */
  /* --------------------------------- */
	readG = PathFile + "startFlow.dat";

	gReadData.open(readG.c_str());
	if (!gReadData.good()) {
		cout << "\n no start values for flow available, using build-in values" << endl;
	}
	else {
		size_t BranchNumber;
		for (size_t
			i = 0; i <= mBr; ++i) {
			gReadData >> BranchNumber;
			gReadData >> Branches[BranchNumber].g;
			if (Branches[BranchNumber].g < 0.) {
				Branches[BranchNumber].reverseDirection();
				Branches[BranchNumber].g = -Branches[BranchNumber].g;
			}
		}
		gReadData.close();
	}//file opened correctly
	 /*     --------------------------------- */
	/** 4) Initial guess for Node Pressure   */
	/*     --------------------------------- */
	for (auto& iNode : Nodes) {
		iNode.pNode = 0.9 * Drum.rhoW * 9.80665 * (100. - iNode.Elev); //drum elevation is 100m
	}
	/* -------------------------------- */
	 /** 5) start of main loop to calculate the flow in branches  */
	 /* -------------------------------- */
 //   Base.showDPTube = true;
 //   Base.showDPBranch = true;
 //   Base.showNodePressure = true;
	for (Base.iterg = 1; Base.iterg <= Base.maxit + 2; ++Base.iterg) {
		if (Base.iterg == Base.maxit) { //show all
			Base.showEnth = true;
			Base.showDPTube = true;
			Base.showDPBranch = true;
			Base.showNodePressure = true;
			Base.showFlow = true;
			Base.showReverse = true;
		}
		/* ----------------- */
		/** 6) save flow to file, renaming this file (...StartFlow.dat) can give a start point for subsequent iterations\n
		 * function: saveFlow() */
		 /* ----------------- */
		saveFlow();
		if (Base.showFlow) {
			Print2dxf(ShowMode::Flow);
			Print2dxf(ShowMode::Arrows);
		}
		/*     ------------------------------------------------------------------ */
		/** 7) Iterative calculation of enthalpies in nodes and enthalpy or steam flow at branch inlets\n
		* function: CalcEnthalpyNodes()  */
		/*     ------------------------------------------------------------------ */
		if (CalcEnthalpyNodes()) break;
		/*     ------------------------------------------------ */
		/**8) Calculation of pressure difference in branches and factors for characteristic curves\n
		*      Model: dp = dPLinear*g + dPConstant\n
		* function: _branch.dpBranch() */
		/*     ------------------------------------------------ */
		maxg = false;
		if (Base.showDPBranch) {
			prot << "\n\n   *************";
			prot << "\n   * dp-Branch *";
			prot << "\n   *************" << endl;
		}

		for (auto& iBranch : Branches) {
			if (Tubes[iBranch.NbTbInBr[0]].xIn > 1.) {
				cout << " XIn in branch " << iBranch.Number << "  > 1." << endl;
			}
			if (Base.showDPBranch) {
				prot << "\n\n   *************";
				prot << "\n   * Branch # " << iBranch.Number << endl;
			}
			if (iBranch.isFlowSet2zero) {
				if (Base.showDPBranch) {
					prot << "\n Flow is set to 0 " << endl;
				}
			}
			else {
				iBranch.dpBranch();
			}
			//         prot << "\n no " << iBranch.Number << " g " << iBranch.g << " dPlin " << iBranch.dPLinear << " dPconst " << iBranch.dPConstant;
			if (Tubes[iBranch.NbTbInBr[0]].xIn > 1.) {
				cout << " XIn in branch " << iBranch.Number << "  > 1." << endl;
			}
		}

		/** 9) summing of dissipation energy (will be needed in later versions of the program) */
		dissipation = TotalDissipation();
		//      saveFlow();

				/* ----------------------------- */
			  /** 10) calculation of nodal pressure for given characteristic curves\n
			  * and calculate new flow from nodal pressure differences and characteristic curves\n
			  * stop program if failure\n
			  * function: ::CalcPressureNodes() */
			  /* ----------------------------- */
		CalcPressureNodes(S, B, X);

		/* ------------------------------ */
		/** 11) comparison of flow in branches to previous iteration step and break loop if result is close enough*/
		/* ------------------------------ */
		if (Base.showFlow) {
			prot << "\n New flow in iteration step " << Base.iterg;
			prot << "\n iBr\tg\tgnew\tdiff";
		}

		maxg = true;
		diffsum = 0.;
		double diffmax = 0.;
		size_t Branchmax = MINUS1;
		for (const auto& iBranch : Branches) {
			diff = fabs(iBranch.gNew - iBranch.g);
			if (fabs(diff / iBranch.gNew) * 100. > diffmax) {
				Branchmax = iBranch.Number;
				diffmax = fabs(diff / iBranch.gNew) * 100.;
			}

			diffsum += diff;

			if (fabs(diff / iBranch.gNew) * 100. > Base.tol) { // Same tolerance in % for each branch flow
				if (iBranch.g / iBranch.minArea > 0.1) { //disregard branches with 0 flow
					//(very small differences lead to very big relative difference)
					maxg = false;
				}
			}
			if (Base.showFlow) {
				prot << "\n" << iBranch.Number << "\t" << iBranch.g << "\t" << iBranch.gNew << "\t" << diff << "\t" << fabs(diff / iBranch.g) * 100. << " maxg " << maxg;
			}
		}
		if (Base.showFlow) {
			prot << "\n Iteration " << Base.iterg << " Flow in branches: sum of differences " << diffsum << " flow in drum " << Nodes[0].gSum;
			prot << "\n" << Branches[Branchmax].Number << "\t" << Branches[Branchmax].g << "\t" << Branches[Branchmax].gNew << "\t" << fabs(Branches[Branchmax].g - Branches[Branchmax].gNew) << "\t" << diffmax << endl;
		}
		cout << "\n Iteration " << Base.iterg << " Flow in branches: sum of differences " << diffsum << " max diff " << diffmax << "% in Branch " << Branchmax << endl;
		//      cout << "\n pnode1 " << Nodes[1].pNode << " full static " << Drum.rhoW * 9.80665 * (100. - Nodes[1].Elev);
		if ((maxg && (fabs(Nodes[0].gSum) < Base.tol)) || diffsum < Base.tol) {
			maxg = true;
			break;
		}
		if (Base.iterg >= Base.maxit) {
			//        Base.showNodePressure = true;
			if (Base.iterg == Base.maxit + 2) {
				cout << "\n\nMaximum number of iterations reached. Is result close enough (y/n)? ";
				cin >> answer;
				if (answer != 'n') {
					maxg = true;
				}
				break;
			}
		}
		/* ------------------------------------------- */
		/** 12) determination of flow for next iteration step as well as reverse of flow direction if needed\n
		* function: Step() */
		/* ------------------------------------------- */
		if (Step()) Print2dxf(ShowMode::Arrows); //if at least one branch reversed show arrows
	}

	/* -------------------------------------------- */
	/** 13) end of iteration loop                  */
	/** 14) once more updating tube data with final flow */
	/* -------------------------------------------- */
	for (auto& iBranch : Branches) {
		if (Base.showDPBranch) {
			prot << "\n\n   *************";
			prot << "\n   * Branch # " << iBranch.Number << endl;
		}
		iBranch.dpBranch();
	}
	if (maxg) {
		/* -------------- */
		/** 16)save results as text file as well as different .dxf files\n
		* function : SaveResults() */
		/* -------------- */
		SaveResults();
	}

	cname = PathFile + ".bcd";
	// Set exceptions to be thrown on failure
	dchange.exceptions(std::ofstream::failbit | std::ofstream::badbit);

	try {
		dchange.open(cname.c_str());
	}
	catch (std::system_error& e) {
		std::cerr << e.code().message() << std::endl;
		exit(1);
	}
	dchange << "\n branch\tdirection changes" << endl;
	for (const auto& iBranch : Branches) {
		if (iBranch.NoChanges > 0) {
			dchange << iBranch.Number << "\t" << iBranch.NoChanges << endl;
		}
	}
	dchange.close();
	if (maxg) {
		cout << "\n successful iteration" << endl;
		prot << "\n successful iteration" << endl;
	}
	else {
		cout << "\n iteration not converging" << endl;
		prot << "\n iteration not converging max. iterations: " << Base.maxit << " iterations:  " << Base.iterg << endl;
	}
	prot.close();

	return 0;
} /* MAIN__ */



