/*****************************************************************//**
 * \file   LES.cpp
 * \brief sets up and solves the linear equation system
 *
 * in each node the flow should be balanced (mass balance)\n
 * as the flow sum in the node depends on the nodal pressure the pressure in all nodes can be determined\n
 * without drum: number of unknown = number of equations
 * \param S reference to system sparse matrix
 * \param B reference to right hand side vector
 * \param X reference to solution vector
 * \return true if successful
 *
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/
//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

bool LES(SparseMatrix<double>& S, /* System matrix as n*? sized vector.*/
	VectorXd& B, /* right hand side vector............*/
	VectorXd& X /* solution vector ...................*/
) // true, if successful
{
	/* Local variables */
	double direction;
	size_t iNd, jNd;
	//  prot << "\n in LES";
	  /**
		* setting of equation system.\n
		* Based on mass balance in all nodes\n
		* As steam drum pressure is fixed, the steam drum is not included\n
		* the size of the square matrix is "number of nodes" - 1
		*/
	S.setZero();
	for (iNd = 1; iNd <= mNd; ++iNd) {
		bool isSeti = false;
		bool isSetj = false;
		// index for Eigen has to be long long, iNd is unsigned long long 
		// in the cast the -1 is also applied for offset to exclude steam drum 
		long long iNdLL = static_cast<long long> (iNd - 1ULL);
		B.coeffRef(iNdLL) = 0.;
		for (const auto& iBr : Nodes[iNd].NbBr) {
			_branch iBranch = Branches[iBr];
			if (fabs(iBranch.dPLinear) > 1e-6) {
				jNd = iBranch.NbNdOut;
				direction = 1.;
				if (jNd == iNd) {
					jNd = iBranch.NbNdIn;
					direction = -1.;
				}
				if (isSeti) {
					S.coeffRef(iNdLL, iNdLL) += 1. / iBranch.dPLinear;
				}
				else {
					S.insert(iNdLL, iNdLL) = 1. / iBranch.dPLinear;
					isSeti = true;
				}
				//            A[iNd-1][iNd-1] += 1./Branch[iBr].dynLin;
				if (jNd > 0) {
					long long jNdLL = static_cast<long long> (jNd - 1ULL);
					if (isSetj) {
						S.coeffRef(iNdLL, jNdLL) -= 1. / iBranch.dPLinear;
					}
					else {
						S.insert(iNdLL, jNdLL) = -1. / iBranch.dPLinear;
						isSetj = true;
					}
					//               A[iNd-1][jNd-1]-=1./Branch[iBr].dynLin;
				}
				B.coeffRef(iNdLL) += direction * iBranch.dPConstant / iBranch.dPLinear;
			}
		}
	}
	if (Base.showNodePressure) {
		prot << "\n Test data (sparse Matrix ):\n";
		prot << S << endl;
		prot << "\n Test data (right hand side):\n";
		prot << B << endl;
	}
	/// calling solver from Eigen library
 //   ConjugateGradient<SparseMatrix<double>, Eigen::Lower | Eigen::Upper  > solver;doesn't work on all models
 //   ConjugateGradient < SparseMatrix<double>, Eigen::Lower|Eigen::Upper, IncompleteCholesky<SparseMatrix<double> > solver; doesn't work at all
 //   BiCGSTAB<SparseMatrix<double>> solver;doesn't work
	SimplicialLDLT<SparseMatrix<double>> solver;

	solver.compute(S);
	if (solver.info() != Success) {
		// decomposition failed
		prot << "\n problem solving equation system: decomposition failed ";
		cout << "\n problem solving equation system: decomposition failed ";
		return false;
	}
	X = solver.solve(B);
	if (solver.info() != Success) {
		// solving failed
		prot << "\n problem solving equation system: solving failed ";
		cout << "\n problem solving equation system: solving failed ";
		return false;
	}
	else {
		/// if calculation successful, copy data of result vector to pressure of nodes "pNode"
		for (iNd = 1; iNd <= mNd; ++iNd) {
			long long iNdLL = static_cast<long long> (iNd - 1ULL);
			Nodes[iNd].pNode = X.coeffRef(iNdLL);
		}
		Nodes[0].pNode = 0.; // pressure in drum set to 0
	}
	return true;
} /* LES */

bool SingleEquation() // true, if successful
{
	/* Local variables */
	double numerator = 0.,
		denominator = 0.;
	  /**
		* setting of equation.\n
		* Based on mass balance in node 1\n
		* As steam drum pressure is fixed, the steam drum is not included\n
		*/
	for (const auto& iBr : Nodes[1].NbBr) {
		_branch iBranch = Branches[iBr];
		if (fabs(iBranch.dPLinear) > 1e-6) {
			if (iBranch.NbNdOut == 1) {
				numerator -= iBranch.dPConstant / iBranch.dPLinear;
			}
			else {
				numerator += iBranch.dPConstant / iBranch.dPLinear;
			}
			denominator += 1. / iBranch.dPLinear;
		}
	}
	Nodes[1].pNode = numerator / denominator;
	Nodes[0].pNode = 0.; // pressure in drum set to 0
	return true;
} /* LES */

