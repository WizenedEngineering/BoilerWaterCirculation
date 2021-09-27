/*****************************************************************//**
 * \file initFlow.cpp
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/
//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"
#include <list>

 /**
 * \struct VertexInfo
 * \brief structure to hold one vertex of shortest distance tree
 *
 */
struct VertexInfo {
	size_t NbNdPrev = MINUS1; ///< index of previous node
	size_t NbBr = MINUS1; ///< Number of Branch between NbNd and previous Node
	double weight = 1e30; ///< Pseudo resistance factor from start node to NbNd
};

/**
 * \fn void printBranch(const _branch& iBranch, const size_t& iNd)
 * \brief prints data of branch to protocol
 *
 * \param iBranch reference to branch iBranch
 * \param iNd reference to node iNd
 */
void printBranch(const _branch& iBranch, const size_t& iNd) {
	prot << "\n iNd " << iNd << " iBr " << iBranch.Number << " deltaH " << iBranch.deltaH;
	switch (iBranch.Kind) {
	case KindOfBranch::Undefined:
		prot << " kind : Undefined";
		break;
	case KindOfBranch::Heated:
		prot << " kind : Heated";
		break;
	case KindOfBranch::Downcomer:
		prot << " kind : Downcomer";
		break;
	case KindOfBranch::Down2Heated:
		prot << " kind : Down2Heated";
		break;
	case KindOfBranch::Heated2Drum:
		prot << " kind : Heated2Drum";
		break;
	case KindOfBranch::Drum2Down:
		prot << " kind : Drum2Down";
		break;
	default:
		prot << " kind : unknown (should not happen)";
		break;
	}
	prot << " node in " << iBranch.NbNdIn << " node out " << iBranch.NbNdOut << endl;

	return;
}

/**
 * \fn void renumberNodes(list<size_t>& NodesHeatedIn, list<size_t>& NodesHeatedOut, list<size_t>& NodesDowncomerIn, list<size_t>& NodesDowncomerOut,
		  vector<bool>& isNdDownIn, vector<bool>& isNdDownOut, vector<bool>& isBrDirectionSet)
 * \brief renumber the nodes
 *
 * renumber the nodes that the numbering tries to follows an estimated flow from drum\n
 * the purpose is to more easily iterate through the Nodes vector for the iterative determination of the enthalpies in the nodes (enth.cpp)\n
 * the new numbering reduces significantly the number of iterations in CalcEnthapyNodes\n
 * lower part of boiler close to start of Nodes vector\n
 * upper part of boiler at end of Nodes vector\n
 * middle part of boiler according elevation\n
 * \param NodesHeatedIn list of numbers of nodes at inlet of heated branches
 * \param NodesHeatedOut list of numbers of nodes at outlet of heated branches
 * \param NodesDowncomerIn list of numbers of nodes at inlet of downcomers
 * \param NodesDowncomerOut list of numbers of nodes at outlet of downcomers
 * \param isNdDownIn vector same size as Nodes, true if node is inlet of downcomer
 * \param isNdDownOut vector same size as Nodes, true if node is outlet of downcomer
 * \param isBrDirectionSet vector same size as Branches, true if direction of branch is set
 */
void renumberNodes(list<size_t>& NodesHeatedIn, list<size_t>& NodesHeatedOut, list<size_t>& NodesDowncomerIn, list<size_t>& NodesDowncomerOut,
	vector<bool>& isNdDownIn, vector<bool>& isNdDownOut, vector<bool>& isBrDirectionSet) {
	stVector NewNodes(mNd + 1, MINUS1), //contains the new node numbers, index corresponds to old number
		OldNodes(mNd + 1, MINUS1); // contains the old node numbers, index corresponds to new number
	stVector Vertex,
		oldVertex,
		lowerStart,
		upperEnd;
	vector<stVector> Vertices;
	size_t otherNode;
	bool isHeated;
	/*for (auto& iNode : Nodes) {
		prot << "\n  renumber 42 " << setw(5) << iNode.Number << "  " << setw(5) << iNode.NbPt;
	}*/
	/**
	 * the new numbers follow the water flow via connection to downcomer -> downcomer -> lower headers->\n
	 * heating surfaces -> upper headers -> connection to drum
	 *
	 * Each branch has 2 nodes, NbNdIn and NbNdOut. We are going through the Nodes-vector, one is the node, we are looking at, "iNd" or "iNode" and the node at *opposite side of branch is called "otherNode".
	 * -# drum keeps drum and is starting point
	 * -# from drum to inlet of downcomer
	 * -# Branches from drum to downcomer are already set( KindOfBranch::Drum2Down), only follow them
	 *    - if other node is inlet node of downcomer it should not be followed because it will be handled in the next step
	 * -# all downcomers\n
	 *    "NodesDowncomerIn" and "NodesDowncomerOut" are determined in calling function. They simply be taken.
	 * -# lower headers from downcomer in until last heatedIn
	 *    - by now "Vertices" is filled with drum, nodes between drum , downcomerIn and downcomerOut\n
	 *    now all nodes from back-end are followed (basically downcomerOut)
	 *    - the last element of "Vertices" is taken and all branches in that node followed
	 *    - if the branch is heated or downcomer we don't follow
	 *    - if the branch is heated we save the node number in "lowerStart"
	 * -# copy what's done so far to new node numbers vector
	 * -# upper headers, overflow etc. from drum to heatedOut going backward
	 *    - have to clear Vertices first
	 *    - starting point is drum again
	 *    - only follow Branches that are not set already
	 *    - if node is not handled (visited) already, push node number to vertex and set KindOfBranch::Heated2Drum
	 * -# copy what's done to new node numbers vector, starting from end of vector
	 * -# all branches in-between, starting from lower header (lowerStart)
	 *    - have to clear Vertices vector first
	 *    - follow the branches from this node only if is not KindOfBranch::Downcomer, KindOfBranch::Down2Heated or KindOfBranch::Drum2Down
	 *    - if node is not handled (visited) already, push node number to vertex
	 * -# copy what's done to new node numbers vector
	 *  - The new numbers are settled now
	 * -# update isNdDownOut\n
	 *    a temporary vector is used to store the values
	 * -# update isNdDownIn\n
	 *    a temporary vector is used to store the values
	 * -# overwrite NodesHeatedIn with new numbers\n
	 *    a temporary vector is used to store the values\n
	 *    a complicated construct\n
	 *    my main point is, that I don't want to lose the pointer to NodesHeatedIn because I refer to it in the calling program
	 * -# overwrite NodesHeatedOut with new numbers\n
	 *    a temporary vector is used to store the values
	 * -# overwrite NodesDowncomerIn (if available) with new numbers\n
	 *    a temporary vector is used to store the values
	 * -# overwrite NodesDowncomerOut with new numbers\n
	 *    a temporary vector is used to store the values
	 * -# update NbNdIn and NbNdOut for all branches
	 * -# update the node number of each point (if it's a node)
	 * -# update the node numbers itself
	 *    - copy at right position in temporary vector
	 *    - copy temporary vector back to Nodes
	 */
	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (auto& iNode : Nodes) {
		iNode.isVisited = false;
	}
	/*
	 * -# drum keeps drum and is starting point
	 */
	Nodes[DRUM].isVisited = true;
	Vertex.push_back(DRUM);
	Vertices.push_back(Vertex);
	/*
	 * -# from drum to inlet of downcomer
	 */
	if (!NodesDowncomerIn.empty()) {
		while (true) {
			Vertex.clear();
			for (const size_t& iNd : Vertices.back()) {
				for (const size_t& iBr : Nodes[iNd].NbBr) {
					if (Base.showMeshDetail) {
						printBranch(Branches[iBr], iNd);
					}
					/*
					 *-# Branches from drum to downcomer are already set( KindOfBranch::Drum2Down), only follow them
					 */
					if (Branches[iBr].Kind == KindOfBranch::Drum2Down) {
						otherNode = Branches[iBr].NbNdOut;
						/*
						 * - if other node is inlet node of downcomer it should not be followed because it will be handled in the next step
						 */
						if (isNdDownIn[otherNode]) {
							if (Base.showMeshDetail) {
								prot << "\n don't follow downcomer ";
							}
							continue;
						}
						else {
							Vertex.push_back(otherNode);
							Nodes[otherNode].isVisited = true;
						}
					}
				}
			}
			if (Vertex.empty()) {
				break;
			}
			else {
				Vertices.push_back(Vertex);
			}
		}; // while (true)

/* -# all downcomers\n
*  "NodesDowncomerIn" and "NodesDowncomerOut" are determined in calling function. They simply be taken.
*/
		Vertex.clear();
		for (const size_t& iNd : NodesDowncomerIn) {
			Vertex.push_back(iNd);
			Nodes[iNd].isVisited = true;
		}
		Vertices.push_back(Vertex);
	}
	if (!NodesDowncomerOut.empty()) {
		Vertex.clear();
		for (const size_t& iNd : NodesDowncomerOut) {
			Vertex.push_back(iNd);
			Nodes[iNd].isVisited = true;
		}
		Vertices.push_back(Vertex);
	}

	/*
	 * -# lower headers from downcomer in until last heatedIn
	 *    - by now "Vertices" is filled with drum, nodes between drum , downcomerIn and downcomerOut\n
	 *    now all nodes from back-end are followed (basically downcomerOut)
	 */
	while (true) {
		oldVertex.clear();
		oldVertex = Vertex;
		Vertex.clear();
		for (const size_t& iNd : oldVertex) {
			isHeated = false;
			/*
			 *    - the last element of "Vertices" is taken and all branches in that node followed
			 */
			for (const size_t& iBr : Nodes[iNd].NbBr) {
				if (Base.showMeshDetail) {
					printBranch(Branches[iBr], iNd);
				}
				/*
				 *   - if the branch is heated or downcomer we don't follow
				 *   - if the branch is heated we save the node number in "lowerStart"
				 */
				if (Branches[iBr].Kind == KindOfBranch::Heated) {//don't follow heated Branches
					isHeated = true;
				}
				else if (Branches[iBr].Kind != KindOfBranch::Downcomer) { // don't follow downcomers
					otherNode = Branches[iBr].NbNdOut;
					if (iNd == otherNode) {
						if (Nodes[Branches[iBr].NbNdIn].isVisited) continue;
						Branches[iBr].reverseDirection();
						otherNode = Branches[iBr].NbNdOut;
					}
					if (!Nodes[otherNode].isVisited) {
						isBrDirectionSet[iBr] = true;
						Vertex.push_back(otherNode);
						//						prot << " push back otherNode " << otherNode << endl;
						Nodes[otherNode].isVisited = true;
						// should be set in initflow
						if (Branches[iBr].Kind == KindOfBranch::Undefined) {
							Branches[iBr].Kind = KindOfBranch::Down2Heated;
						}
					}
				}
			}
			if (isHeated) {
				lowerStart.push_back(iNd);
			}
		}
		if (Vertex.empty()) {
			break;
		}
		else {
			Vertices.push_back(Vertex);
		}
	}; // while (true)
/*
 * -# copy what's done so far to new node numbers vector
 */
	size_t i = 0; // i is needed later on
	for (stVector iV : Vertices) {
		for (size_t iNd : iV) {
			NewNodes[iNd] = i;
			OldNodes[i] = iNd;
			i++;
		}
	}
	/* NodesHeatedIn already contains all the content of lowerStart
	The elements of lowerStart should be at end of NodesHeatedIn where the nodes closest to downcomer should be at end. for this we have to reverse lowerStart*/

	//   reverse(lowerStart.begin(), lowerStart.end());  
	std::vector<size_t>::reverse_iterator rit = lowerStart.rbegin();
	for (; rit != lowerStart.rend(); ++rit) {
		NodesHeatedIn.remove(*rit);
		NodesHeatedIn.push_back(*rit);
	}
	/*
	 * -# upper headers, overflow etc. from drum to heatedOut going backward
	 *  - have to clear Vertices first
	 */
	for (stVector iV : Vertices) {
		iV.clear();
	}
	Vertices.clear();
	/*
	 * - starting point is drum again
	 */
	Vertex.push_back(DRUM);
	Vertices.push_back(Vertex);
	while (true) {
		oldVertex.clear();
		oldVertex = Vertex;
		Vertex.clear();
		for (size_t iNd : oldVertex) {
			isHeated = false;
			for (size_t& iBr : Nodes[iNd].NbBr) {
				if (Base.showMeshDetail) {
					printBranch(Branches[iBr], iNd);
				}
				/*
				 *   - only follow Branches that are not set already
				 *   - if node is not handled (visited) already, push node number to vertex and set KindOfBranch::Heated2Drum
				 */
				if (Branches[iBr].Kind == KindOfBranch::Heated) { // don't follow heated branches 
					isHeated = true;
				}
				else if (Branches[iBr].Kind != KindOfBranch::Downcomer && // don't follow downcomers 
					Branches[iBr].Kind != KindOfBranch::Drum2Down && // don't follow Drum->downcomers 
					Branches[iBr].Kind != KindOfBranch::Down2Heated) { // don't follow downcomers->heated 
					otherNode = Branches[iBr].NbNdIn;
					if (iNd == otherNode) {
						if (Nodes[Branches[iBr].NbNdOut].isVisited) continue;
						Branches[iBr].reverseDirection();
						otherNode = Branches[iBr].NbNdIn;
					}
					if (!Nodes[otherNode].isVisited) {
						isBrDirectionSet[iBr] = true;
						Vertex.push_back(otherNode);
						Branches[iBr].Kind = KindOfBranch::Heated2Drum;
						//						prot << " push back otherNode " << otherNode << endl;
						Nodes[otherNode].isVisited = true;
					}
				}
			}
			if (isHeated && iNd != DRUM) {
				upperEnd.push_back(iNd);
			}
		}
		if (Vertex.empty()) {
			break;
		}
		else {
			Vertices.push_back(Vertex);
		}
	}; // while (true)
/*
 *  -# copy what's done to new node numbers vector, starting from end of vector
 */
	size_t j = mNd;
	for (stVector iV : Vertices) {
		for (size_t iNd : iV) {
			if (iNd != 0) {
				NewNodes[iNd] = j;
				OldNodes[j] = iNd;
				j--;
			}
		}
	}
	/* NodesHeatedIn already contains all the content of lowerStart
	The elements of lowerStart should be at end of NodesHeatedIn where the nodes closest to downcomer should be at end*/
	std::vector<size_t>::reverse_iterator riter = upperEnd.rbegin();
	for (; riter != upperEnd.rend(); ++riter) {
		NodesHeatedOut.remove(*riter);
		NodesHeatedOut.push_back(*riter);
	}
	/*
 *   -# all branches in-between, starting from lower header (lowerStart)
 *   - have to clear Vertices vector first
 */
	for (stVector iV : Vertices) {
		iV.clear();
	}
	Vertices.clear();

	Vertices.push_back(lowerStart);
	while (true) {
		Vertex.clear();
		for (size_t iNd : Vertices.back()) {
			for (size_t& iBr : Nodes[iNd].NbBr) {
				if (Base.showMeshDetail) {
					printBranch(Branches[iBr], iNd);
				}
				/*
				 *  - follow the branches from this node only if is not KindOfBranch::Downcomer, KindOfBranch::Down2Heated or KindOfBranch::Drum2Down
				 *   - if node is not handled (visited) already, push node number to vertex
				 */
				if (Branches[iBr].Kind != KindOfBranch::Downcomer && // no downcomer
					Branches[iBr].Kind != KindOfBranch::Down2Heated && // no lower header
					Branches[iBr].Kind != KindOfBranch::Drum2Down) { // no connection drum ->downcomer
					otherNode = Branches[iBr].NbNdIn;
					if (iNd == otherNode) otherNode = Branches[iBr].NbNdOut;
					if (Nodes[otherNode].isVisited) continue;
					Vertex.push_back(otherNode);
					//					prot << " push back otherNode " << otherNode << endl;
					Nodes[otherNode].isVisited = true;
				}
			}
		}
		if (Vertex.empty()) {
			break;
		}
		else {
			Vertices.push_back(Vertex);
		}
	} // while (true)
/*
  *   -# copy what's done to new node numbers vector
  */
	if (Vertices.size() > 1) {
		for (auto iV = Vertices.begin(); iV != Vertices.end(); ++iV) {
			if (iV != Vertices.begin()) {
				for (size_t iNd : *iV) {
					NewNodes[iNd] = i;
					OldNodes[i] = iNd;
					i++;
				}
			}
		}
	}
	//prot << " NewNodes ";
	//for (int iNd : NewNodes) {
	//	prot << "  " << iNd;
	//}
	//prot << endl;
	//prot << " OldNodes ";
	//for (int iNd : OldNodes) {
	//	prot << "  " << iNd;
	//}
	//prot << endl;

	//prot << "\n\n before renumbering isNdDownOut ";
	//for (bool iisNdDownOut : isNdDownOut)
	//	prot << " " << setw(5) << iisNdDownOut;
	//prot << endl;
/*
 *  - The new numbers are settled now
 *  - update isNdDownOut\n
 *    a temporary vector is used to store the values
 */

	size_t Size = isNdDownOut.size();
	vector<bool> temp1(Size);

	for (i = 0; i < Size; i++) {
		temp1[i] = isNdDownOut[OldNodes[i]];
	}
	for (i = 0; i < Size; i++) {
		isNdDownOut[i] = temp1[i];
	}
	//prot << "\n\n after renumbering isNdDownOut ";
	//for (bool iisNdDownOut : isNdDownOut)
	//	prot << " " << setw(5) << iisNdDownOut;
	//prot << endl;
/*
  *   - update isNdDownIn\n
  *    a temporary vector is used to store the values
  */

  //		Size is same as isNdDownOut 
	for (i = 0; i < Size; i++) {
		temp1[i] = isNdDownIn[OldNodes[i]];
	}
	for (i = 0; i < Size; i++) {
		isNdDownIn[i] = temp1[i];
	}
	temp1.clear();

	//prot << "\n\n before renumbering nodeHeatedIn ";
	//for (int &iNodeHeatedIn : NodesHeatedIn)
	//	prot << " " << setw(5) << iNodeHeatedIn;
	//prot << endl;
/*
  *   - overwrite NodesHeatedIn with new numbers\n
  *    a temporary vector is used to store the values\n
  *    a complicated construct\n
  *     my main point is, that I don't want to lose the pointer to NodesHeatedIn because I refer to it in the calling program
  */
	stVector temp2;
	for (auto it = NodesHeatedIn.begin(); it != NodesHeatedIn.end(); ++it) {
		temp2.push_back(NewNodes[*it]);
	}
	NodesHeatedIn.clear();
	for (auto it = temp2.begin(); it != temp2.end(); ++it) {
		NodesHeatedIn.push_back(*it);
	}
	//   prot << "\n\n after renumbering nodeHeatedIn ";
	//   for (size_t &iNodeHeatedIn : NodesHeatedIn)
	//   	prot << " " << setw(5) << iNodeHeatedIn;
	//   prot << endl;
		//prot << "\n\n before renumbering nodeHeatedOut ";
		//for (int &iNodeHeatedOut : NodesHeatedOut)
		//	prot << " " << setw(5) << iNodeHeatedOut;
		//prot << endl;
	/*
	  *   - overwrite NodesHeatedOut with new numbers
	  *    a temporary vector is used to store the values
	  */

	temp2.clear();
	for (auto it = NodesHeatedOut.begin(); it != NodesHeatedOut.end(); ++it) {
		temp2.push_back(NewNodes[*it]);
	}
	NodesHeatedOut.clear();
	for (auto it = temp2.begin(); it != temp2.end(); ++it) {
		NodesHeatedOut.push_back(*it);
	}
	//  prot << "\n\n after renumbering nodeHeatedOut ";
	//  for (size_t &iNodeHeatedOut : NodesHeatedOut)
	//  	prot << " " << setw(5) << iNodeHeatedOut;
  //   prot << endl;
  /*
	*  - overwrite NodesDowncomerIn (if available) with new numbers
	*    a temporary vector is used to store the values
	*/
	if (!NodesDowncomerIn.empty()) {

		//prot << "\n\n before renumbering nodeDowncomerIn ";
		//for (int &iNodeDowncomerIn : NodesDowncomerIn)
		//	prot << " " << setw(5) << iNodeDowncomerIn;
		//prot << endl;

		temp2.clear();
		for (auto it = NodesDowncomerIn.begin(); it != NodesDowncomerIn.end(); ++it) {
			temp2.push_back(NewNodes[*it]);
		}
		NodesDowncomerIn.clear();
		for (auto it = temp2.begin(); it != temp2.end(); ++it) {
			NodesDowncomerIn.push_back(*it);
		}
		//prot << "\n\n after renumbering NodesDowncomerIn ";
		//for (int &iNodeDowncomerIn : NodesDowncomerIn)
		//	prot << " " << setw(5) << iNodeDowncomerIn;
		//prot << endl;
	}
	//prot << "\n\n before renumbering NodesDowncomerOut ";
	//for (int &iNodeDowncomerOut : NodesDowncomerOut)
	//	prot << " " << setw(5) << iNodeDowncomerOut;
	//prot << endl;
/*
  *   - overwrite NodesDowncomerOut with new numbers
  *    a temporary vector is used to store the values
  */

	temp2.clear();
	for (auto it = NodesDowncomerOut.begin(); it != NodesDowncomerOut.end(); ++it) {
		temp2.push_back(NewNodes[*it]);
	}
	NodesDowncomerOut.clear();
	for (auto it = temp2.begin(); it != temp2.end(); ++it) {
		NodesDowncomerOut.push_back(*it);
	}
	//prot << "\n\n after renumbering NodesDowncomerOut ";
	//for (int &iNodeDowncomerOut : NodesDowncomerOut)
	//	prot << " " << setw(5) << iNodeDowncomerOut;
	//prot << endl;
	//	prot << "\n Branch  node in  node out" << endl;
/*
 *   - update NbNdIn and NbNdOut for all branches
 */
	for (auto& iBranch : Branches) {
		iBranch.NbNdIn = NewNodes[iBranch.NbNdIn];
		iBranch.NbNdOut = NewNodes[iBranch.NbNdOut];
		//	prot << setw(5) << iBranch.Number << setw(5) << iBranch.NbNdIn << setw(5) << iBranch.NbNdOut << endl;
	}

	//prot << "\nPoint   NbNd   " << endl;
/*
  *   - update the node number of each point (if it's a node)
  */

	for (auto& iPoint : Points) {
		if (iPoint.NbNd != MINUS1) {
			iPoint.NbNd = NewNodes[iPoint.NbNd];
		}
		//		prot << setw(5) << iPoint.Number << setw(5) << iPoint.NbNd << endl;
	}

	//prot << "\n Number   Point" << endl;
	//for (auto& iNode : Nodes) {
	//	prot << setw(5) << iNode.Number << setw(5) << iNode.NbPt << endl;
	//}
/*
 *   - copy at right position in temporary vector
 */
	vector <_node> tempNodes(mNd + 1);
	for (auto& iNode : Nodes) {
		tempNodes[NewNodes[iNode.Number]] = iNode;
		tempNodes[NewNodes[iNode.Number]].Number = NewNodes[iNode.Number];
	}
	/*
	  *   - copy temporary vector back to Nodes
	  */
	Nodes = tempNodes;
	tempNodes.clear();
	//prot << "\n  Number   Point" << endl;
	//for (auto& iNode : Nodes) {
	//	prot << setw(5) << iNode.Number << setw(5) << iNode.NbPt << endl;
	//}
	return;
}

// comparison of elevation, highest first

/**
 * \fn bool compare_elevation_high(const size_t& first, const size_t& second)
 * \brief comparison of node elevation (high)
 *
 * true if elevation of first node is higher than second node
 * \param first: index of first node
 * \param second: index of second node
 * \return bool true if elevation of first node is higher than second node
 */
bool compare_elevation_high(const size_t& first, const size_t& second) {
	if (Nodes[first].Elev > Nodes[second].Elev) {
		return true;
	}
	else {
		return false;
	}
}

/**
 * \fn bool compare_elevation_low(const size_t& first, const size_t& second)
 * \brief comparison of node elevation (low)
 *
 * true if elevation of first node is lower than second node
 * \param first: index of first node
 * \param second: index of second node
 * \return bool true if elevation of first node is lower than second node
 */
bool compare_elevation_low(const size_t& first, const size_t& second) {
	if (Nodes[first].Elev < Nodes[second].Elev) {
		return true;
	}
	else {
		return false;
	}
}

/**
 * \fn void printVisitedVertices(const vector<VertexInfo>& Vertices, size_t iNode)
 * \brief writes info of visited vertices (nodes) to protocol file
 *
 * for debugging
 * \param Vertices: vector of vertices
 * \param iNode: number of starting node (root of tree)
 */
void printVisitedVertices(const vector<VertexInfo>& Vertices, size_t iNode) {
	prot << "\n Starting Node " << iNode << "\n";
	for (size_t i = 0; i <= mNd; i++) {
		if (Vertices[i].weight < 1e20) {
			prot << std::setw(5) << i << " , " << Vertices[i].NbNdPrev << "," << Vertices[i].NbBr << "," << Vertices[i].weight << " | ";
		}
		if (i % 6 == 0) prot << "->\n ";
	}
	prot << "\n---------------------------------------------------------------------------------------\n";
	return;
}

/**
 * \fn void put_into(list<size_t>& sorted, const vector<VertexInfo>& Vertices, size_t iVx, double weight)
 * \brief puts iVx (node number) into sorted list and keep it sorted according weight with lowest weight at beginning of list
 *
 * \param sorted: list of indexes (iVx, node number)
 * \param Vertices: vector of all Vertices
 * \param iVx: index = node number
 * \param weight: weight of Vertex iVx
 */
void put_into(list<size_t>& sorted, const vector<VertexInfo>& Vertices, size_t iVx, double weight) {
	/**
	 * -# if list empty put it in
	 * -# if iVx is already in list, it might be at the wrong position -> remove it
	 * -# iterate through list\n
	 *    if weight of list element is higher than weight, insert the iVx at this place
	 * -# if no place found to insert, put at end of list
	 */
	if (sorted.empty()) {
		sorted.push_back(iVx);
		return;
	}
	sorted.remove(iVx);
	for (auto iter = sorted.begin(); iter != sorted.end(); iter++) {
		if (Vertices[*iter].weight > weight) {
			iter = sorted.insert(iter, iVx);
			return;
		}
	}
	sorted.push_back(iVx);
	return;
}
#if 0
/*
* \fn void remove_one(list<size_t>& theList, const size_t number)
 * \brief removes one element from list
 *
 * list element with content >number< should be removed
 * \param theList: the list
 * \param number: list element with content "number" should be removed
 */
void remove_one(list<size_t>& theList, const size_t number) {
	for (auto iter = theList.begin(); iter != theList.end(); iter++) {
		if (*iter == number) {
			theList.erase(iter);
			return;
		}
	}
}
#endif
/**
 *\fn void initialize_Vertices(vector<VertexInfo>& Vertices)
 *  \brief initialize vector of vertices
 *
 * \param Vertices: vector of Vertices
 */
void initialize_Vertices(vector<VertexInfo>& Vertices) {
	for (auto& iVertex : Vertices) {
		iVertex.NbNdPrev = MINUS1;
		iVertex.NbBr = MINUS1;
		iVertex.weight = 1e30;
	}
	return;
}

/**
 * \brief process node at other end of branch
 *
 * \param Vertices: vector of Vertices
 * \param tempVx: list of temporary Vertices
 * \param iBr: branch number
 * \param iNd: node number of starting node
 * \param NodalFlow: Flow at root of tree
 * \param R: two-phase factor for pressure drop
 * \param minWeightTarget: minimum weight of target
 * \return size_t index of other node
 */
size_t process_OtherNode(vector<VertexInfo>& Vertices, list<size_t>& tempVx,
	const size_t& iBr, const size_t& iNd, const double NodalFlow,
	const double R, const double minWeightTarget) {
	/**
	 * otherNode is the node at other end of branch no matter if it's inlet or outlet node\n
	 * the flow direction has not been settled yet and can change later on
	 *
	 * the new weight of other node is weight of this node + "weight of branch"
	 *
	 * the new weight is compared to already minimum weight of other node\n
	 * if it is lower the data of the node (weight, NbNdPrev, NbBr) is replaced\n
	 * and it is put at the right place in the sorted list tempVx
	 */
	size_t otherNode = Branches[iBr].NbNdOut;
	if (iNd == otherNode) otherNode = Branches[iBr].NbNdIn;
	double Weight = Vertices[iNd].weight + Branches[iBr].dPConstant * (Branches[iBr].g * Branches[iBr].g + R * NodalFlow * NodalFlow);

	if (Base.showMeshDetail) {
		prot << "\notherNode " << otherNode << " iNd Weight " << Vertices[iNd].weight << " otherNodeWeight " <<
			Vertices[otherNode].weight << " Weight   " << Weight << " minweightTarget " << minWeightTarget;
		prot << "\n dPconst " << Branches[iBr].dPConstant << " g " << Branches[iBr].g << " NodalFlow " << NodalFlow;
	}
	if (Weight < minWeightTarget) {
		if (Vertices[otherNode].weight > Weight) {
			Vertices[otherNode].NbNdPrev = iNd;
			Vertices[otherNode].NbBr = iBr;
			put_into(tempVx, Vertices, otherNode, Weight);
			Vertices[otherNode].weight = Weight;
			return otherNode;
		}
	}
	return MINUS1;
}

int setDowncomer(_branch& iBranch, vector<bool>& isBrDirectionSet, list<size_t>& NodesDowncomerOut,
	vector<bool>& isNdDownOut, list<size_t>& NodesDowncomerIn, list<size_t>& NodesDowncomerIn2) {

	iBranch.Kind = KindOfBranch::Downcomer;
	if (Base.showMeshDetail) {
		prot << "\n set downcomer iBranch " << setw(5) << iBranch.Number;
	}
	if (iBranch.deltaH > 0.) {
		iBranch.reverseDirection();
	}
	isBrDirectionSet[iBranch.Number] = true;
	//typically there is a check if a node is downcomer outlet
	//an array same size as Node can be handled easier than several times traversing a vector
	//vector is needed for renumbering
	isNdDownOut[iBranch.NbNdOut] = true;

	//   if (NodesDowncomerOut.empty()) {
	NodesDowncomerOut.push_back(iBranch.NbNdOut);
	//   }
	//   else {
	//      bool already_exist = false;
	//      for (const auto& iNodeDowncomerOut : NodesDowncomerOut) {
	//         if (iNodeDowncomerOut == iBranch.NbNdOut) {
	//            already_exist = true;
	//            break;
	//         }
	//      }
	//      if (!already_exist) {
	//         NodesDowncomerOut.push_back(iBranch.NbNdOut);
	//      }
	//   }

	Nodes[iBranch.NbNdOut].gSum = 0.;
	if (iBranch.NbNdIn != DRUM) {//if start point of downcomer already drum
		// we don't need to look for connections between drum and downcomer
		NodesDowncomerIn.push_back(iBranch.NbNdIn);
		//      NodesDowncomerIn2.push_back(iBranch.NbNdIn);
	}
	// if several downcomers end in same node, this node should only appear once in NodesDowncomerOut and NodesDowncomerIn
	if (!NodesDowncomerIn.empty()) {
		NodesDowncomerIn.sort();
		NodesDowncomerIn.unique();
		NodesDowncomerIn2 = NodesDowncomerIn;
	}
	NodesDowncomerOut.sort();
	NodesDowncomerOut.unique();

	return 0;
}


/**
 * \brief function to determine the kind of branches and initial flow
 *
 *  the weight takes into account the flow from previously processed paths
*/
void initFlow() {
	/* Local variables */
	size_t iNd, iNdIn, iNdOut;
	size_t iBr;
	size_t otherNode;
	size_t iNdSource;
	double Flow = 0.;
	list<size_t> NodesHeatedIn;
	list<size_t> NodesHeatedOut;
	list<size_t> NodesDowncomerIn, //it is consumed during the search for connections between drum and downcomer
		NodesDowncomerIn2; // to be used in renumber
	list<size_t> NodesDowncomerOut;

	vector<bool> isNdDownOut(mNd + 1, false);
	vector<bool> isNdDownIn(mNd + 1, false);
	vector<bool> isBrDirectionSet(mBr + 1, false);

	vector<VertexInfo> Vertices(mNd + 1);
	list<size_t> tempVx;
	//	VertexInfo PaVx;
	size_t iNdTarget;
	list<int> PathNodes;
	vector<size_t> Drum2Down;
	vector<vector < size_t>> Drum2AllDown;
	double minWeightTarget = 1e30;

	/** The flow direction at least for first iteration step should be physically possible
	 * (each node needs to have branches arriving and leaving
	 * as well as there should exist paths between outlet of heated branches and drum,
	 * inlet of heated branches and downcomer, drum and inlet of downcomers)
	 * also the mass balance at the nodes should be satisfied. Otherwise problems in calculation of branch inlet enthalpies
	 *
	 * as paths are determined they can be used for the initial guess of the flow
	 *
	 * 1. determination of downcomers (are not heated and pass the given level or\n
	 *    have a given inlet enthalpy below 0 for heated downcomers) kind = KindOfBranch::Downcomer\n
	 *    Flow in downcomer is determined at the end, after all flows at NodesDowncomerOut are known\n
	 *    in further calculation there is a check if a node is downcomer outlet\n
	 *    an bool vector same size as Nodes can be handled easier than several times traversing the Nodes vector\n
	 *    this vector is also needed for renumbering\n
	 *    if several downcomers end in same node, this node should only appear once in NodesDowncomerOut
	 * 2. connection between drum and downcomer (if downcomer is not directly connected to drum)) kind = KindOfBranch::Drum2Downcomer\n
	 *    for connection between NodesDowncomerIn and drum it's a straight Dijkstra algorithm <a href="https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm">Dijkstra algorithm on Wikipedia</a> \n
	 *    Flow in downcomer to be determined later. Flow in Branches passed by different paths is summed up.
	 * 3. all heated branches except heated downcomers should have upward flow (except horizontal heated branches)  kind = KindOfBranch::Heated
	 *    - flow is assumed at a circulation ratio that is input value
	 *    - flow direction has to be settled first (upward) before adding to inlet and outlet node. Direction of flow in horizontal heated branches will be handled later
	 *    - if inlet of heated branch already drum we don't need to look for connections between downcomer and heated branch
	 *    - if endpoint of heated branch already drum we don't need to look for connections between heated branch and drum
	 *    - NodesHeatedIn and NodesHeatedOut can contain multiple entries of same node\n
	 *      got problems during sorting according elevation -> have to take out double entries
	 *    - the lists need to be sorted: for NodesHeatedIn the highest elevation first, for NodesHeatedOut the lowest elevation first
	 * 4. renumbering of nodes this makes life easier later on
	 * 5. connection between downcomers and heated branches (flow downcomer -> heated branches) kind = KindOfBranch::Down2Heated\n
	 *    A path has to be found between downcomer(s) and inlet of heated branches.
	 *    It is not the shortest path from Rotterdam to Groningen, but some kind of: Dijkstra algorithm is used in this case.\n
	 *   (it's more like the shortest path from Nijmegen to the next fishing port, no matter which particular one, at least for the connection between NodesHeatedIn and one downcomer)
	 *    - check if NodesHeatedIn is already NodesDowncomerOut -> no need look for path
	 *    - have to find a path from iNodeHeatedIn to one iNodeDowncomerOut
	 *       - don't follow heated branches up
	 *       - don't follow Branches that flow is already set against
	 *       - don't follow downcomers
	 *       - search tree ends when a downcomer outlet node is hit
	 *       - walk back. Flow in Branches passed by different paths is summed up. If flow direction is wrong, change direction
	 * 6. we summed up all flows in each downcomer out\n
	 *    it can happen that multiple downcomers end at the same iNodeDowncomerOut -> the flow will be distributed equally to the downcomers
	 * 7. branches connected to outlet of downcomers are set\n
	 *    also flow in downcomers are set\n
	 *    now the flow from drum to downcomers can be set
	 * 8. connection between heated branch -> drum  kind = KindOfBranch::Heated2Drum\n
	 *    connection between NodesHeatedOut (single point) and drum (also single point) it's a straight Dijkstra algorithm\n
	 *       - as the flow is mixture of water and steam a two-phase factor is added to the pseudo resistance
	 *       - don't follow heated branches down
	 *       - don't follow downcomers
	 *       - don't follow Branches that flow is already set against the momentary flow direction
	 *       - when drum is found -> walk back and set flows\n
	 * 9. set upward flow in all remaining non-horizontal branches\n
	 *    this also covers connection pipes that will not be touched by the paths
	 * 10. check for nodes that have no arriving or no leaving branches
	 *    - reversing direction of one branch that is able to change without violation at the other end
	 *    - at least for the start we can't change the flow direction in heated upward branches and downcomers
	 * 11. looking for branches that are parallel, i.e. same start node and same end node\n
	 *    we can improve the "path" method on those places\n
	 *    basically the flow in each step is quite random, this can lead to one branch having "heavy" weight\n
	 *    but in the next step a big flow is coming that goes to the "lighter" branch making it far more "heavy" than the other one\n
	 *    in the parallel branches the total flow should be kept but distributed according pseudo pressure drop
	 */

	if (Base.showMeshDetail) {
		prot << "\n\n before downcomer ?????????????????????????????????????????????????????????????????????????????????????????\n\n";
		std::cout << "\n\n before downcomer";
	}

	// 1) determination of downcomers (are not heated and pass center level of heated branches or 
	//    have a given inlet enthalpy below 0 for heated downcomers) kind = KindOfBranch::Downcomer
	for (auto& iBranch : Branches) {
		//      prot << "\n Branch " << iBranch.Number << " elevIn " << Nodes[iBranch.NbNdIn].Elev <<" elevOut "<< Nodes[iBranch.NbNdOut].Elev << " level " << Base.DowncomerLevel;
		if (iBranch.qSum < 1e-3 &&
			((Nodes[iBranch.NbNdOut].Elev < Base.DowncomerLevel &&
				Nodes[iBranch.NbNdIn].Elev > Base.DowncomerLevel) ||
				(Nodes[iBranch.NbNdOut].Elev > Base.DowncomerLevel &&
					Nodes[iBranch.NbNdIn].Elev < Base.DowncomerLevel))) {
			// flow to be handled later
			setDowncomer(iBranch, isBrDirectionSet, NodesDowncomerOut, isNdDownOut,
				NodesDowncomerIn, NodesDowncomerIn2);
		}
		else { //check for heated downcomers (f.i. bi-drum) criterion: inlet enthalpy in downcomer has to be lower than drum water enthalpy
			bool isDown = false;
			for (auto& iTb : iBranch.NbTbInBr) {
				if (Tubes[iTb].EnthInGiven < 0.) {
					isDown = true;
					break;
				}
			}
			if (isDown) {
				setDowncomer(iBranch, isBrDirectionSet, NodesDowncomerOut, isNdDownOut,
					NodesDowncomerIn, NodesDowncomerIn2);
			}
		}
	}
	if (Base.showMeshDetail) {
		std::cout << "\n after downcomers ";
		prot << "\n\n after downcomers nodeDowncomerIn ";
		for (auto& iNodeDowncomerIn : NodesDowncomerIn)
			prot << " " << setw(5) << iNodeDowncomerIn;
		prot << endl;
		prot << "\n\n nodeDowncomerOut ";
		for (int i = 1; i <= mNd; i++)
			if (isNdDownOut[i]) prot << " " << setw(5) << i;
		prot << endl;
	}
	if (Base.showMeshDetail) {
		prot << "\n\n before drum ->downcomer  ?????????????????????????????????????????????????????????????????????????????????????????\n\n";
		std::cout << "\n\n before drum ->downcomer ";
	}

	// 2) connection between drum and inlet of downcomer (if not direct connected to drum)) kind = KindOfBranch::Drum2Down
	/* for connection between NodesDowncomerIn and drum it's a straight Dijkstra algorithm */
	size_t iNodeDowncomerIn;
	while (!NodesDowncomerIn.empty()) {
		iNodeDowncomerIn = NodesDowncomerIn.front();
		if (Base.showMeshDetail) {
			prot << "\nroot " << iNodeDowncomerIn;
		}
		initialize_Vertices(Vertices);
		tempVx.push_back(iNodeDowncomerIn);
		Vertices[iNodeDowncomerIn].weight = 0.;
		while (!tempVx.empty()) {
			iNd = tempVx.front();
			tempVx.pop_front();
			if (Base.showMeshDetail) {
				prot << "\n---------------------------------------------------------------------------------------\n";
			}
			if (iNd != DRUM) { //iNd is drum -> already target
				for (size_t& jBr : Nodes[iNd].NbBr) {
					_branch* jBranch = &Branches[jBr];
					if (Base.showMeshDetail) {
						printBranch(Branches[jBr], iNd);
					}
					if (jBranch->Kind == KindOfBranch::Downcomer && iNd == jBranch->NbNdIn) { //don't follow downcomers
						if (Base.showMeshDetail) {
							prot << "\n don't follow downcomer ";
						}
						continue;
					}
					otherNode = process_OtherNode(Vertices, tempVx, jBr, iNd, Flow, 1., minWeightTarget);
					if (otherNode == DRUM) {//->drum
						tempVx.clear();
						break;
					}
				}
			}
		}
		if (Base.showMeshDetail) {
			printVisitedVertices(Vertices, iNodeDowncomerIn);
			prot << endl;
		}
		iNdSource = iNodeDowncomerIn;
		iNd = DRUM;
		Drum2Down.push_back(iNodeDowncomerIn);
		while (iNd != iNdSource) {
			otherNode = Vertices[iNd].NbNdPrev;
			iBr = Vertices[iNd].NbBr;
			if (Base.showMeshDetail) {
				prot << "start iBr" << iBr << endl;
			}
			if (Branches[iBr].NbNdOut == otherNode) {
				Branches[iBr].reverseDirection();
			}
			Drum2Down.push_back(iBr);
			if (Branches[iBr].Kind == KindOfBranch::Undefined) {
				Branches[iBr].Kind = KindOfBranch::Drum2Down;
			}
			iNd = otherNode;
		}
		Drum2AllDown.push_back(Drum2Down);
		Drum2Down.clear();
		NodesDowncomerIn.remove(iNdSource);
	}//end of NodesDowncomerIn
	if (Base.showMeshDetail) {
		prot << "\n\n before heated branches ?????????????????????????????????????????????????????????????????????????????????????????\n\n";
		std::cout << "\n\n before heated branches ";
	}
	// 3) :all heated branches should have upward flow kind = heated (except horizontal heated branches ) 
	// heated downcomers are already handled
	for (auto& iBranch : Branches) {
		if (iBranch.qSum > 1e-3 && !(iBranch.Kind == KindOfBranch::Downcomer)) {
			iBranch.Kind = KindOfBranch::Heated;
			iBranch.g = iBranch.qSum / Drum.enthEvap * Base.CRStart; //Circulation ratio as input value
			if (Base.showMeshDetail) {
				prot << "\n iBranch " << setw(5) << iBranch.Number;
			}
			if (fabs(iBranch.deltaH) > 1e-3) {
				if (iBranch.deltaH < 0.) {
					iBranch.reverseDirection();
				}
				// direction of flow in horizontal heated branches will be handled later
				isBrDirectionSet[iBranch.Number] = true;
			}
			else {
				iBranch.isHorizontalHeated = true;
			}
			//flow direction has to be settled first before adding to inlet and outlet node
			Nodes[iBranch.NbNdIn].gSumLeave += iBranch.g;
			Nodes[iBranch.NbNdOut].gSumArrive += iBranch.g;

			if (iBranch.NbNdIn != DRUM) { // add to list
				NodesHeatedIn.push_back(iBranch.NbNdIn);
			}

			if (iBranch.NbNdOut != DRUM) {//if endpoint of heated branch already drum
				// we don't need to look for connections between heated branch and drum
				NodesHeatedOut.push_back(iBranch.NbNdOut);
			}
		}
	}
	// NodesHeatedIn and NodesHeatedOut can contain multiple entries of same node
	// got problems during sorting according elevation ->
	// have to take out double entries
	NodesHeatedIn.sort();
	NodesHeatedIn.unique();
	NodesHeatedOut.sort();
	NodesHeatedOut.unique();

	// the lists need to be sorted: for NodesHeatedIn the highest elevation first, for NodesHeatedOut the lowest elevation first
	try {
		NodesHeatedIn.sort(compare_elevation_high);
		NodesHeatedOut.sort(compare_elevation_low);
	}
	catch (std::exception& e) {
		std::cout << e.what() << '\n';
	}


	if (Base.showMeshDetail) {
		std::cout << "\nafter heated branches nodeHeatedIn ";
		prot << "\n\n after heated branches nodeHeatedIn ";
		for (auto& iNodeHeatedIn : NodesHeatedIn)
			prot << " " << setw(5) << iNodeHeatedIn;
		prot << endl;
		prot << "\n nodeHeatedOut ";
		for (auto& iNodeHeatedOut : NodesHeatedOut)
			prot << " " << setw(5) << iNodeHeatedOut;
		prot << endl;
	}
	//4)here re-numbering of nodes
	renumberNodes(NodesHeatedIn, NodesHeatedOut, NodesDowncomerIn2, NodesDowncomerOut, isNdDownIn, isNdDownOut, isBrDirectionSet);

	//prot << "\n\n after back nodeHeatedIn ";
	//for (auto &iNodeHeatedIn : NodesHeatedIn)
	//	prot << " " << setw(5) << iNodeHeatedIn;
	//prot << endl;
	// 5)connection between downcomers and heated branches (flow downcomer -> heated branches) kind = KindOfBranch::Down2Heated
	if (Base.showMeshDetail) {
		prot << "\n\n before downcomer ->heated branches 5 ?????????????????????????????????????????????????????????????????????????????????????????\n\n";
		std::cout << "\n\n before downcomer ->heated branches 5";
	}
	/* A path has to be found between downcomer(s) and inlet of heated branches

	It is not the shortest path from Rotterdam to Groningen, but some kind of: Dijkstra algorithm is used in this case
	(it's more like the shortest path from Nijmegen to the next fishing port, no matter which particular one,
	at least for the connection between NodesHeatedIn and one downcomer)
	 */

	while (!NodesHeatedIn.empty()) {
		iNdTarget = MINUS1;
		minWeightTarget = 1e30;
		size_t iNodeHeatedIn = NodesHeatedIn.front();
		NodesHeatedIn.pop_front();
		if (Base.showMeshDetail) {
			prot << "\nroot " << iNodeHeatedIn;
		}
		Flow = Nodes[iNodeHeatedIn].gSumLeave;
		initialize_Vertices(Vertices);
		// check if iNodeHeatedIn is already NodesDowncomerOut -> no need look for path
		if (isNdDownOut[iNodeHeatedIn]) { // iNodeHeatedIn = NodeDowncomerOut
			if (Base.showMeshDetail) {
				prot << "\niNodeHeatedIn " << iNodeHeatedIn << " is equal downcomerOut";
			}
			// the flow in heated tubes has to be added to flow in NodeDowncomerOut for the flow in downcomer
			Nodes[iNodeHeatedIn].gSum += Nodes[iNodeHeatedIn].gSumLeave; //gSum = flow at downcomer out
		}
		else {
			// have to look for a path from iNodeHeatedIn to one iNodeDowncomerOut
			tempVx.push_back(iNodeHeatedIn);
			Vertices[iNodeHeatedIn].weight = 0.;
			while (!tempVx.empty()) {
				iNd = tempVx.front();
				tempVx.pop_front();
				if (Base.showMeshDetail) {
					prot << "\n---------------------------------------------------------------------------------------\n";
				}
				if (isNdDownOut[iNd]) continue;  //NdDownOut is already a target, don't try to search for other connections from here
				for (size_t& jBr : Nodes[iNd].NbBr) {
					_branch* jBranch = &Branches[jBr];
					if (Base.showMeshDetail) {
						printBranch(Branches[jBr], iNd);
					}
					if (jBranch->Kind == KindOfBranch::Heated &&
						iNd == jBranch->NbNdIn &&
						!jBranch->isHorizontalHeated) {
						if (Base.showMeshDetail) {
							prot << "\n don't follow heated up";
						}
						continue; //don't follow heated branches up
					}
					if (isBrDirectionSet[jBr] && iNd == jBranch->NbNdIn) { // don't follow Branches that flow is already set against
						if (Base.showMeshDetail) {
							prot << "\n don't follow against set flow ";
						}
						continue; // don't follow Branches that flow is already set against
					}
					if (jBranch->Kind == KindOfBranch::Downcomer &&
						iNd == jBranch->NbNdIn) { //don't follow downcomers
						if (Base.showMeshDetail) {
							prot << "\n don't follow downcomers ";
						}
						continue;
					}
					otherNode = process_OtherNode(Vertices, tempVx, jBr, iNd, Flow, 1., minWeightTarget);
					if (otherNode < MINUS1 && isNdDownOut[otherNode]) {
						iNdTarget = otherNode;
						minWeightTarget = Vertices[iNdTarget].weight;
					}
				}
			}
			if (Base.showMeshDetail) {
				printVisitedVertices(Vertices, iNodeHeatedIn);
				prot << endl;
			}
			if (iNdTarget == MINUS1) { // no downcomer found
				processError(Nodes[iNodeHeatedIn].NbPt, "downcomer to heated branch");
			}
			//walk_back
			Flow = Nodes[iNodeHeatedIn].gSumLeave;

			Nodes[iNdTarget].gSum += Flow;
			iNd = iNdTarget;
			//         prot << "\n walk back iNdTarget " << iNdTarget << " heatedIn " << iNodeHeatedIn << endl;
			while (iNd != iNodeHeatedIn) {
				otherNode = Vertices[iNd].NbNdPrev;
				std::reverse(Nodes[otherNode].NbBr.begin(), Nodes[otherNode].NbBr.end());
				// if there are 2 parallel branches between downcomer and last node before, 
				// the branch that comes in the vector NbBr first will be preferred
				// therefore after being in final path the vector will be reversed  
				iBr = Vertices[iNd].NbBr;
				if (Base.showMeshDetail) {
					prot << "\nstart iBr" << iBr << " initial " << Branches[iBr].g;
				}
				if (Branches[iBr].NbNdIn == otherNode) { //flow from iNd -> otherNode
					Branches[iBr].reverseDirection();
					if (Branches[iBr].Kind == KindOfBranch::Heated) {
						//                     remove_one(NodesHeatedIn, otherNode);
						NodesHeatedIn.remove(otherNode);
						NodesHeatedOut.push_back(otherNode);
						NodesHeatedOut.sort(compare_elevation_low);
						double g = Branches[iBr].qSum / Drum.enthEvap * Base.CRStart;
						Nodes[iNd].gSumLeave += g;
						Nodes[iNd].gSumArrive -= g;
						Nodes[otherNode].gSumLeave -= g;
						Nodes[otherNode].gSumArrive += g;
					}
				}
				isBrDirectionSet[iBr] = true;
				Branches[iBr].g += Flow;
				if (Branches[iBr].Kind == KindOfBranch::Undefined) {
					Branches[iBr].Kind = KindOfBranch::Down2Heated;
				}
				//           prot << " new " << Branches[iBr].g ;

				iNd = otherNode;
			}
		} //end iNodeHeatedIn != NDownOut
	}//end of NodesHeatedIn
	if (Base.showMeshDetail) {
		prot << "\n after  nodeHeatedIn ";
		std::cout << "\n after  nodeHeatedIn ";
		for (auto& jNodeHeatedIn : NodesHeatedIn)
			prot << " " << setw(5) << jNodeHeatedIn;
		prot << " should be empty " << endl;
	}

	// we summed up all flows in each downcomer out
	// it can happen that multiple downcomers end at the same iNodeDowncomerOut
	// the flow will be distributed equally to the downcomers
//   int NoDowncomerInNode = 0;
	for (iNd = 1; iNd <= mNd; iNd++) {
		if (isNdDownOut[iNd]) {
			int NoDowncomerInNode = 0;
			for (size_t& jBr : Nodes[iNd].NbBrArrive) {
				if (Branches[jBr].Kind == KindOfBranch::Downcomer) {
					NoDowncomerInNode++;
				}
			}
			for (size_t& jBr : Nodes[iNd].NbBrArrive) {
				_branch* jBranch = &Branches[jBr];
				if (jBranch->Kind == KindOfBranch::Downcomer) {
					jBranch->g += Nodes[iNd].gSum / NoDowncomerInNode;
					Nodes[jBranch->NbNdIn].gSumLeave += jBranch->g;
				}
			}
		}
	}
	if (Base.showMeshDetail) {
		prot << "\n after nodeDowncomerOut ";
		std::cout << "\n after nodeDowncomerOut ";
		prot << endl;
		prot << "\n direction not set \n";
		for (iBr = 0; iBr <= mBr; iBr++) {
			if (!isBrDirectionSet[iBr]) prot << "  " << iBr;
		}
		prot << "\n\n before Heated -> drum  ?????????????????????????????????????????????????????????????????????????????????????????\n";
		std::cout << "\n\n before Heated -> drum ";
		prot << "\n NodesHeatedOut ";
		for (auto& iNodeHeatedOut : NodesHeatedOut)
			prot << " " << setw(5) << iNodeHeatedOut;
		prot << endl;
		Print2dxf(ShowMode::BranchesNodes);
		//std::exit(1234);
	}
	// branches connected to lower side of downcomers are set
	// also flow in downcomers set
	//	now the flow from drum to downcomers can be set
	for (auto& conn : Drum2AllDown) {
		Flow = Nodes[conn[0]].gSumLeave;
		for (auto it = conn.begin() + 1; it != conn.end(); ++it) {
			Branches[*it].g += Flow;
		}
		conn.clear();
	}
	Drum2AllDown.clear();

	// 4) : connection between heated branch -> drum  kind = Heated2Drum

	/* for connection between NodesHeatedOut and drum it's a straight Dijkstra algorithm
	as it is the connection between start and target A* algorithm can be used also
	 */
	double R = 1. + 2400. * pow(1. / (Base.CRStart * Drum.pMPa / 0.0980665), 0.96); //two-phase-factor according to Becker
	// we don't need to treat the case of NodeHEatedOut.empty because all heated branches already end in drum-> no need to search for connecting branches
	while (!NodesHeatedOut.empty()) {
		size_t iNodeHeatedOut = NodesHeatedOut.front();
		NodesHeatedOut.pop_front();
		if (Base.showMeshDetail) {
			prot << "\nroot " << iNodeHeatedOut;
		}
		Flow = Nodes[iNodeHeatedOut].gSumArrive;
		initialize_Vertices(Vertices);
		// have to look for a path or paths from iNodeHeatedOut to drum
		tempVx.push_back(iNodeHeatedOut);
		Vertices[iNodeHeatedOut].weight = 0.;
		minWeightTarget = 1e30;
		while (!tempVx.empty()) {
			iNd = tempVx.front();
			tempVx.pop_front();
			if (iNd == DRUM) continue;   // don't go further from drum 
			if (Base.showMeshDetail) {
				prot << "\n---------------------------------------------------------------------------------------\n";
			}
			for (size_t& jBr : Nodes[iNd].NbBr) {
				_branch* jBranch = &Branches[jBr];
				if (Base.showMeshDetail) {
					printBranch(Branches[jBr], iNd);
				}
				if (jBranch->Kind == KindOfBranch::Heated &&
					iNd == jBranch->NbNdOut &&
					!jBranch->isHorizontalHeated) {
					if (Base.showMeshDetail) {
						prot << "\n don't follow heated down";
					}
					continue; //don't follow heated branches down
				}
				if (jBranch->Kind == KindOfBranch::Downcomer) { //don't follow downcomers
					if (Base.showMeshDetail) {
						prot << "\n don't follow downcomer ";
					}
					continue;
				}
				if (isBrDirectionSet[jBr] && iNd == jBranch->NbNdOut) {
					if (Base.showMeshDetail) {
						prot << "\n don't follow against set flow ";
					}
					continue; // don't follow Branches that flow is already set against
				}
				otherNode = process_OtherNode(Vertices, tempVx, jBr, iNd, Flow, R, minWeightTarget);
				if (!otherNode) {//otherNode == 0 ->drum
					minWeightTarget = Vertices[0].weight;
					//              break;
				}
			}
			//	}
		}
		tempVx.clear();
		if (Base.showMeshDetail) {
			printVisitedVertices(Vertices, iNodeHeatedOut);
			prot << endl;
		}
		Flow = Nodes[iNodeHeatedOut].gSumArrive;
		iNd = DRUM;
		// we are starting from drum and want to go back 
				  // if drum has no previous node on this path there must be an error 
		if (Vertices[DRUM].NbNdPrev == MINUS1) {
			processError(Nodes[iNodeHeatedOut].NbPt, "heated branch to drum");
		}
		while (iNd != iNodeHeatedOut) {
			otherNode = Vertices[iNd].NbNdPrev;
			iBr = Vertices[iNd].NbBr;
			_branch* iBranch = &Branches[iBr];
			if (Base.showMeshDetail) {
				prot << "start iBr" << iBr << endl;
			}
			if (iBranch->NbNdOut == otherNode) {
				iBranch->reverseDirection();
				if (iBranch->Kind == KindOfBranch::Heated) {
					//              remove_one(NodesHeatedOut, otherNode);
					NodesHeatedOut.remove(otherNode);
					NodesHeatedIn.push_back(otherNode);
					NodesHeatedIn.sort(compare_elevation_high);
					double g = iBranch->qSum / Drum.enthEvap * Base.CRStart;
					Nodes[iNd].gSumLeave -= g;
					Nodes[iNd].gSumArrive += g;
					Nodes[otherNode].gSumLeave += g;
					Nodes[otherNode].gSumArrive -= g;
				}
				//					remove_one(NodesHeatedOut, otherNode);
			}
			isBrDirectionSet[iBr] = true;
			iBranch->g += Flow;
			//         prot << "\n iBr " << iBr << " g " << iBranch->g;
			if (iBranch->Kind == KindOfBranch::Undefined) {
				iBranch->Kind = KindOfBranch::Heated2Drum;
			}
			iNd = otherNode;
		}
		//      prot << "\n=============================================================";
	}//end of NodesHeatedOut
	if (Base.showMeshDetail) {
		std::cout << "\n after nodeHeatedOut ";
		prot << "\n after nodeHeatedOut ";
		for (auto& jNodeHeatedOut : NodesHeatedOut)
			prot << " " << setw(5) << jNodeHeatedOut;
		prot << " should be empty " << endl;
	}

	// set upward flow in all remaining non-horizontal branches
	// this also covers connection pipes that will not be touched by the paths

	for (auto& iBranch : Branches) {
		if (iBranch.Kind == KindOfBranch::Undefined) {
			if (iBranch.deltaH < -1e-3) {
				iBranch.reverseDirection();
			}
		}
	}

	// check for nodes that have no arriving or no leaving branches
	// at least for the start we can't change the flow direction in heated upward branches
	// and downcomers
	for (auto& iBranch : Branches) {
		if ((iBranch.Kind != KindOfBranch::Heated &&
			iBranch.Kind != KindOfBranch::Downcomer) ||
			iBranch.isHorizontalHeated) {
			if ((Nodes[iBranch.NbNdIn].NbBrArrive.size() == 0 &&
				Nodes[iBranch.NbNdOut].mBrArrive > 0) ||
				(Nodes[iBranch.NbNdOut].NbBrLeave.size() == 0 &&
					Nodes[iBranch.NbNdIn].mBrLeave > 0)) {
				if (Base.showMeshDetail) {
					prot << "\n reverseDirection iBranch due to no arriving or no leaving branches : " << iBranch.Number << endl;
				}
				iBranch.reverseDirection();
			}
		}
	}

	// looking for branches that are parallel, i.e. same start node and same end node
	// we can improve the "path" method on those places
	// basically the flow in each step is quite random, this can lead to one branch having "heavy" weight
	// but in the next step a big flow is coming that goes to the "lighter" branch making it far more "heavy" than the other one
	// in the parallel branches the total flow should be kept but distributed according pseudo pressure drop

	for (iBr = 0; iBr < mBr; iBr++) {
		_branch* iBranch = &Branches[iBr];
		iNdIn = iBranch->NbNdIn;
		iNdOut = iBranch->NbNdOut;
		for (size_t jBr = iBr + 1ULL; jBr <= mBr; jBr++) {
			_branch* jBranch = &Branches[jBr];
			if ((iNdIn == jBranch->NbNdIn) &&
				(iNdOut == jBranch->NbNdOut)) {
				double gTotal = iBranch->g + jBranch->g;
				double root = sqrt(jBranch->dPConstant / iBranch->dPConstant);
				iBranch->g = root * gTotal / (1. + root);
				jBranch->g = gTotal - iBranch->g;
				//                        prot << "\n parallel iBr" << iBr << " jBr " << jBr;
			}
		}
	}
	if (Base.showMeshDetail) {
		prot << "\n  nodes mass balance ";
		for (auto& iNode : Nodes) {
			iNode.gSumArrive = 0.;
			iNode.gSumLeave = 0.;
			for (const auto& Arrive : iNode.NbBrArrive) {
				iNode.gSumArrive += Branches[Arrive].g;
			}
			for (const auto& Leave : iNode.NbBrLeave) {
				iNode.gSumLeave += Branches[Leave].g;
			}
			prot << "\n iNd " << iNode.Number << " diff " << iNode.gSumArrive - iNode.gSumLeave;
		}
	}
	//	Base.iterg = 2;
	//	Print2dxf(PathFile, 1);
	//	std::exit(1234);

	return;
}
