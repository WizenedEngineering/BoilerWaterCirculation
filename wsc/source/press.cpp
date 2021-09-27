/**
* \file press.cpp
* @brief calculation of pressure in nodes and a resulting new flow in the branches
*
* @param S  reference to System sparse matrix
* @param B  reference to right hand side vector
* @param X  reference to solution vector
* @return int error code
*/
//#include "stdafx.h"

#undef MAINFUNCTION
#include "CommonHeader.h"

int CalcPressureNodes(SparseMatrix<double>& S, /* System sparse matrix */
	VectorXd& B, /* right hand side vector............*/
	VectorXd& X /* solution vector ...................*/
) // error code as return value
{
	/// to reduce numerical errors (small differences of big numbers)
	/// the drum pressure in Nodes[0].pNode is set to 0. ;\n
	/// the absolute pressure doesn't matter because at this place we only
	/// deal with differences to calculate the flow

	/// calling function LES to set up and solve the linear equation system or 
	/// calling function SingleEquation (for bidrum boilers with only both drums as nodes)

	if (mNd == 1) {
		if (!SingleEquation()) {
			prot << "\nerror in solver of equation " << endl;
			std::cout << "\nerror in solver of equation " << endl;
			exit(1);
		}
	}
	else if (!LES(S, B, X)) {
		prot << "\nerror in solver of equation system " << endl;
		std::cout << "\nerror in solver of equation system " << endl;
		exit(1);
	}

	if (Base.showNodePressure) {
		prot << "\n   *****************";
		prot << "\n   * Node Pressure *";
		prot << "\n   *****************\n";
		for (auto& iNode : Nodes) {
			prot << "\n iNd " << iNode.Number << " press " << iNode.pNode << " previous " <<
				iNode.pPrev << " diff " << iNode.pNode - iNode.pPrev;
		}
	}
	for (auto& iNode : Nodes) {
		iNode.pPrev = iNode.pNode;
		iNode.gSum = 0.;
	}
	/// with the new nodal pressures new flows in the branches are calculated 
	for (auto& iBranch : Branches) {
		if (fabs(iBranch.dPLinear) > 1e-6) {
			iBranch.gNew = ((Nodes[iBranch.NbNdIn].pNode - Nodes[iBranch.NbNdOut].pNode) - iBranch.dPConstant) / iBranch.dPLinear;
			Nodes[iBranch.NbNdIn].gSum += iBranch.gNew;
			Nodes[iBranch.NbNdOut].gSum -= iBranch.gNew;
		}
		else {
			iBranch.gNew = 0.;
		}
		if (Base.showNodePressure) {
			prot << "\niBr " << iBranch.Number << " g " << iBranch.g << " gNew " << iBranch.gNew << " deltap " << (Nodes[iBranch.NbNdIn].pNode - Nodes[iBranch.NbNdOut].pNode);
		}
	}
	/// checking if the flow in all nodes is balanced
	for (auto& iNode : Nodes) {
		if (fabs(iNode.gSum) > 0.01) {
			prot << "\nerror in solver of equation system: Flow in node " << iNode.Number << " not balanced" << endl;
			cout << "\nerror in solver of equation system: Flow in node " << iNode.Number << " not balanced" << endl;
			prot << "\n Node " << iNode.Number << " gsum " << iNode.gSum;
			exit(1);
		}
	}
	//if ( Base.showNodePressure ) {
	//   for ( auto& iNode : Nodes ) {
	//      if ( fabs ( iNode.gSum ) > 0.01 ) {
	//         prot << "\n Node " << iNode.Number << " gsum " << iNode.gSum;
	//      }
	//   }
	//}

	return 0;
} /* druck_ */

