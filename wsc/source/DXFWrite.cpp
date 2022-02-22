/*!
 * \file DXFWrite.cpp
 * \author Rainer_Jordan@<very, very warm>mail.com
 * ## Licence
 * Licensed under the EUPL, Version 1.2 or later
 * \fn int DXFWrite ( ostream& outData, ShowMode mode )
 * \brief writing mesh geometry to .dxf file according mode
 *
 * ShowMode::TubesPoints,     : showing tubes and points\n
 * ShowMode::BranchesNodes,   : showing tubes and branch numbers, coloring according kind of branch\n
 * ShowMode::Flow,            : showing tubes, coloring according kind of branch, flow in layer\n
 * ShowMode::Arrows           : showing arrows of flow direction only\n
 * ShowMode::SafetyFactor     : showing tubes, coloring according SafetyFactor\n
 *          red = below 1, yellow = 1.0 ... 1.2, green = good\n
 * ShowMode::SteamQuality     : showing tubes, coloring according steam quality\n
 * red = above 0.2, orange = between 0.2 ... 0.1, yellow = 0.1 ... 1/15, light green = 1/15...0.05, dark green < 0.05\n
 * ShowMode::VoidFraction     : showing tubes, coloring according void fraction\n
 *            red = above 0.95, yellow = 0.85..0.95, green < 0.85\n
 * ShowMode::MassVel          : showing tubes, coloring according mass velocity\n
 *            rainbow from red (lowest) to blue (highest)\n
 * ShowMode::FlowPattern      : showing tubes, coloring according flow pattern (Steiner)\n
 *            red = stratified, wavy or mist, green = slug,plug, bubble or annular\n
 * ShowMode::DPdynPerLength       : showing tubes, coloring according dpdyn / length (to find bottle neck)\n
 *            rainbow from blue (lowest) to red (highest)\n
 * ShowMode::ResistanceFactor : showing tubes, coloring according dpdyn / length / mass velocity^2 (to find bottle neck)\n
 *            rainbow from blue (lowest) to red (highest)\n

 * \param outData outstream of .dxf file
 * \param mode what should be shown in .dxf file
 * \return int error
 */
#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

using namespace std;

extern int DXFWriteTube(ostream& outData, const _tube& iTube, unsigned long long& handle,
	const string& LayerName, unsigned short color, bool writeArrow);

extern int DXFWritePoint(ostream& outData, unsigned long long& handle,
	const string& LayerName, unsigned short color, double xPos, double yPos, double zPos);

extern int DXFWriteText(ostream& outData, unsigned long long& handle, const string& LayerName, unsigned short color,
	double PointMidX, double PointMidY, double PointMidZ, double height, const string& content);

extern int DXFWriteArrow(ostream& outData, unsigned long long& handle,
	unsigned short color, double PointMidX, double PointMidY,
	double PointMidZ, double dx, double dy, double dz);

extern _point OCS2WCS(double OCSCenterX, double OCSCenterY, double OCSCenterZ,
	double Nx, double Ny, double Nz,
	double OCSRadius, double OCSAngle);

	int DXFWrite(ostream& outData,
	ShowMode mode
	// ShowMode::TubesPoints, : showing tubes and points
	// ShowMode::BranchesNodes: showing tubes and coloring according kind of branch
	// ShowMode::Flow, : showing tubes, coloring according kind of branch, flow in layer
	// ShowMode::Arrows       : showing arrows of flow direction only
	// ShowMode::SafetyFactor : showing tubes, coloring according SafetyFactor
			 //    red = below 1, yellow = 1.0 ... 1.2, green = good
	// ShowMode::SteamQuality : showing tubes, coloring according steam quality
		 //    red = above 0.2, orange = between 0.2 ... 0.1, yellow = 0.1 ... 1/15,
			//    light green = 1/15...0.05, dark green < 0.05
	// ShowMode::VoidFraction : showing tubes, coloring according void fraction
		  //    red = above 0.95, yellow = 0.85..0.95, green < 0.85
	// ShowMode::MassVel      : showing tubes, coloring according mass velocity
		 //    rainbow from red (lowest) to blue (highest)
	// ShowMode::FlowPattern  : showing tubes, coloring according flow pattern(Steiner)
		  //    red = stratified, wavy or mist, green = slug,plug, bubble or annular
	// ShowMode::DPdynPerLength   : showing tubes, coloring according dpdyn / length(to find bottle neck)
		 //    rainbow from blue (lowest) to red (highest)
	// ShowMode::ResistanceFactor : showing tubes, coloring according dpdyn / length / mass velocity^2(to find bottle neck)
		 //    rainbow from blue (lowest) to red (highest)
) {
	string LayerName, content;
	int error = 0;
	//  int i;
	double ValueMin = 1e30, ValueMax = 0., ValueStep, Value;
	unsigned short color = 0, TEN = 10;
	//center point of line
	double PointMidX, PointMidY, PointMidZ, dx, dy, dz;
	double OCSMidAngle;
	_point* PtIn;
	_point* PtOut;
	// font height
	double height;
	//using new handles
	unsigned long long handle = 100ULL;
	//determine max and min extension
	double xExtMax = -200000., yExtMax = -200000., zExtMax = -200000.,
		xExtMin = 200000., yExtMin = 200000., zExtMin = 200000.;
	for (const auto& iPoint : Points) {
		xExtMax = fmax(iPoint.xCoord, xExtMax);
		yExtMax = fmax(iPoint.yCoord, yExtMax);
		zExtMax = fmax(iPoint.zCoord, zExtMax);
		xExtMin = fmin(iPoint.xCoord, xExtMin);
		yExtMin = fmin(iPoint.yCoord, yExtMin);
		zExtMin = fmin(iPoint.zCoord, zExtMin);
	}
	height = (zExtMax - zExtMin) / 100.;
	outData << "999\nWSC\n  0\nSECTION\n  2\nCLASSES\n  0\nENDSEC\n  0\nSECTION\n  2\nENTITIES\n";
	switch (mode) {
	case ShowMode::TubesPoints:
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				content = to_string(iTube.Number);
				LayerName = "Tube No " + content;
				error = DXFWriteTube(outData, iTube, handle,
					LayerName, color, true);
				if (error) return error;
			}
		}
		break;
	case ShowMode::BranchesNodes:
		for (const auto& iTube : Tubes) {
			content = to_string(iTube.NbBr);
			LayerName = "Br " + content;
			switch (Branches[iTube.NbBr].Kind)
			{
			case KindOfBranch::Undefined:
				color = 0;
				break;
			case KindOfBranch::Heated:
				color = 1;
				break;
			case KindOfBranch::Downcomer:
				color = 2;
				break;
			case KindOfBranch::Down2Heated:
				color = 3;
				break;
			case KindOfBranch::Heated2Drum:
				color = 4;
				break;
			case KindOfBranch::Drum2Down:
				color = 5;
				break;
			default:
				color = 0;
				break;
			}
			error = DXFWriteTube(outData, iTube, handle,
				LayerName, color,  true);
			if (error) return error;
		}
		for (const auto& iNode : Nodes) {
			LayerName = "Nd " + to_string(iNode.Number);
			error = DXFWritePoint(outData, handle, LayerName, 0,
				Points[iNode.NbPt].xCoord,
				Points[iNode.NbPt].yCoord,
				Points[iNode.NbPt].zCoord);
			if (error) return error;
		}
		break;
	case ShowMode::Flow:
		for (const auto& iTube : Tubes) {
			if (!Branches[iTube.NbBr].isFlowSet2zero) {
				if (iTube.Flow > ValueMax)ValueMax = iTube.Flow;
				if (iTube.Flow < ValueMin)ValueMin = iTube.Flow;
			}
		}
		ValueStep = ceil((ValueMax - ValueMin) / 17.);
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				content = to_string(iTube.Number);
				LayerName = "Tube " + content + " " + to_string(Branches[iTube.NbBr].g / iTube.NoParallel);
				//         prot << "msdmin " << ValueMin << " msdmax " <<ValueMax << " step " << ValueStep << endl;
				color = 0;
				if (Branches[iTube.NbBr].g / iTube.NoParallel > 1e-5) {
					for (unsigned short ic = 1; ic <= 17; ic++) {
						if (Branches[iTube.NbBr].g / iTube.NoParallel < (ValueMin + ic * ValueStep)) {
							color = ic * TEN;
							break;
						}
					}
				}
				error = DXFWriteTube(outData, iTube, handle,
					LayerName, color,(!Branches[iTube.NbBr].isFlowSet2zero) );
				if (error) return error;
			}
		}
		break;
	case ShowMode::Arrows:
		for (const auto& iTube : Tubes) {
			if (!Branches[iTube.NbBr].isFlowSet2zero) {
			   PtIn = &Points[iTube.PointIn];
				PtOut = &Points[iTube.PointOut];
				if (iTube.RadiusBend < 1.) { // straight tube
					PointMidX = (PtIn->xCoord + PtOut->xCoord) / 2.;
					PointMidY = (PtIn->yCoord + PtOut->yCoord) / 2.;
					PointMidZ = (PtIn->zCoord + PtOut->zCoord) / 2.;
				}
				else { // bend
					OCSMidAngle = (iTube.OCSStartAngle + iTube.OCSEndAngle) / 2.;
					_point WCSMidPoint = OCS2WCS(iTube.OCSCenterX, iTube.OCSCenterY, iTube.OCSCenterZ,
							iTube.Nx, iTube.Ny, iTube.Nz, iTube.RadiusBend, OCSMidAngle);
					PointMidX = WCSMidPoint.xCoord;
					PointMidY = WCSMidPoint.yCoord;
					PointMidZ = WCSMidPoint.zCoord;
				}
				dx = PtIn->xCoord - PtOut->xCoord;
				dy = PtIn->yCoord - PtOut->yCoord;
				dz = PtIn->zCoord - PtOut->zCoord;

				error = DXFWriteArrow(outData, handle, color,
						PointMidX, PointMidY, PointMidZ, dx, dy, dz);
				if (error) return error;
			}
		}
		error = DXFWritePoint(outData, handle, "0", 0,
			100000., 100000., 100000.);
		break;
	case ShowMode::SafetyFactor:
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				content = to_string(iTube.Number);
				LayerName = "Tb " + content + " " + to_string(iTube.SafetyFactor);
	         color = 0;
	         if (iTube.q > 0.) {
					if (iTube.SafetyFactor < 1.) color = 1;
					else if (iTube.SafetyFactor < 1.2) color = 2;
					else color = 3;
				}
				//            prot << "\n" << LayerName << " color " << color;
				error = DXFWriteTube(outData, iTube, handle,
					LayerName, color, (!Branches[iTube.NbBr].isFlowSet2zero));
				if (error) return error;
			}
		}
		break;
	case ShowMode::SteamQuality:
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				content = to_string(iTube.Number);
				LayerName = "Tb " + content + " " + to_string(round(iTube.xOut * 1e3) / 1e3)
					+ " CR " + to_string(round(10./iTube.xOut)/ 10.);
				color = 0;
				if (!Branches[iTube.NbBr].isFlowSet2zero) {
					if (iTube.xOut > 0.2) color = 1;
					else if (iTube.xOut > 0.1) color = 30;
					else if (iTube.xOut > 1. / 15) color = 2;
					else if (iTube.xOut > 0.05) color = 3;
					else if (iTube.xOut > 0.0) color = 150;
				}
				// red = above 0.2 (CR: <5), orange = between 0.2 ... 0.1 (CR: 5...10), yellow = 0.1 ... 1/15 (CR: 10...15),
				// light green = 1/15...0.05 (CR: 15...20), dark green < 0.05 (CR: >20)
				error = DXFWriteTube(outData, iTube, handle,
					LayerName, color,(!Branches[iTube.NbBr].isFlowSet2zero) );
				if (error) return error;
			}
		}
		break;
	case ShowMode::VoidFraction:
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				content = to_string(iTube.Number);
				LayerName = "Tb " + content + " " + to_string(round(iTube.VoidFractionOut * 1e3) / 1e3);
				color = 0;
				if (!Branches[iTube.NbBr].isFlowSet2zero) {
					if (iTube.VoidFractionOut > 0.95) color = 1;
					else if (iTube.VoidFractionOut > 0.85) color = 2;
					else if (iTube.VoidFractionOut > 0.6) color = 3;
					else if (iTube.VoidFractionOut > 1e-6) color = 150;
					else color = 0;
					// red > 0.95, yellow > 0.85, light green > 0.6, dark green > 1e-6
				}
				error = DXFWriteTube(outData, iTube, handle,
					LayerName, color,(!Branches[iTube.NbBr].isFlowSet2zero));
				if (error) return error;
			}
		}
		break;
	case ShowMode::MassVel:
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				if (!Branches[iTube.NbBr].isFlowSet2zero) {
					if (iTube.MassVel > ValueMax)ValueMax = iTube.MassVel;
					if (iTube.MassVel < ValueMin)ValueMin = iTube.MassVel;
				}
			}
		}
		ValueStep = ceil((ValueMax - ValueMin) / 17.);
		//         prot << "msdmin " << ValueMin << " msdmax " <<ValueMax << " step " << ValueStep << endl;
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				content = to_string(iTube.Number);
				LayerName = "Tb " + content + " " + to_string(round(iTube.MassVel)) ;
		   	color = 0;
				if (!Branches[iTube.NbBr].isFlowSet2zero) {
					for (unsigned short ic = 1; ic <= 17; ic++) {
						if (iTube.MassVel < (ValueMin + ic * ValueStep)) {
							color = ic * TEN;
							break;
						}
					}
				}
				//                        prot << "\n Itb " << iTube.Number << " MassVel " << iTube.MassVel << " color " << color << endl;
		// rainbow from red (lowest) to blue (highest)
				error = DXFWriteTube(outData, iTube, handle,
					LayerName, color, (!Branches[iTube.NbBr].isFlowSet2zero));
				if (error) return error;
			}
		}
		break;
	case ShowMode::Velocity:
		//for (const auto& iTube : Tubes) {
		//	if (!Branches[iTube.NbBr].isFlowSet2zero) {
		//		if (iTube.velOut > ValueMax)ValueMax = iTube.velOut;
		//		if (iTube.velOut < ValueMin)ValueMin = iTube.velOut;
		//	}
		//}
		ValueMax = 5.;
		ValueMin = 0.;
		ValueStep = (ValueMax - ValueMin) / 17.;
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				content = to_string(iTube.Number);
				LayerName = "Tb " + content + " " + to_string(round(iTube.velOut * 100.) / 100.);
				if (!Branches[iTube.NbBr].isFlowSet2zero) {
					color = 10;
					for (unsigned short ic = 1; ic <= 17; ic++) {
						if (iTube.velOut > (ValueMax - ic * ValueStep)) {
							color = ic * TEN;
							break;
						}
					}
				}
				else {
					color = 0;
				}
				//                        prot << "\n Itb " << iTube.Number << " MassVel " << iTube.MassVel << " color " << color << endl;
								// rainbow from red (lowest) to blue (highest)
				error = DXFWriteTube(outData, iTube, handle,
					LayerName, color, (!Branches[iTube.NbBr].isFlowSet2zero));
				if (error) return error;
			}
		}
		break;
	case ShowMode::FlowPattern:
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				content = to_string(iTube.Number);
				if (!Branches[iTube.NbBr].isFlowSet2zero) {
					if (iTube.xOut > 1e-3 && fabs(iTube.Height / iTube.Length) < 0.174) {
						switch (iTube.steiner) {
						case FlowPattern::stratified:
							LayerName = "Tb " + content + " stratified";
							color = 1;
							break;
						case FlowPattern::wavy:
							LayerName = "Tb " + content + " wavy";
							color = 1;
							break;
						case FlowPattern::bubble:
							LayerName = "Tb " + content + " bubble";
							color = 3;
							break;
						case FlowPattern::plugSlug:
							LayerName = "Tb " + content + " plug or slug";
							color = 3;
							break;
						case FlowPattern::mist:
							LayerName = "Tb " + content + " mist";
							color = 1;
							break;
						case FlowPattern::annular:
							LayerName = "Tb " + content + " annular";
							color = 3;
							break;
						default:
							LayerName = "Tb " + content + " ---";
							color = 0;
						}
					}
					else {
						LayerName = "Tb " + content + " ---";
						color = 0;
					}
				}
				else {
					LayerName = "Tb " + content + " no flow";
					color = 1;
				}
				error = DXFWriteTube(outData, iTube, handle,
					LayerName, color,(!Branches[iTube.NbBr].isFlowSet2zero));
				if (error) return error;
			}
		}
		break;
	case ShowMode::DPdynPerLength:
	//
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				if (!(iTube.DiaOrificeIn > 0. || iTube.DiaOrificeOut > 0.)) { //disregard orifices
					Value = iTube.dpdyn / iTube.Length;
					if (Value > ValueMax)ValueMax = Value;
					if (Value < ValueMin)ValueMin = Value;
				}
			}
		}
		ValueStep = ceil((ValueMax - ValueMin) / 17.);
		//         prot << "\n\nmsdmin " << ValueMin << " msdmax " << ValueMax << " step " << ValueStep << endl;
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				content = to_string(iTube.Number);
				LayerName = "Tb " + content + " " + to_string(round(iTube.dpdyn / iTube.Length));
				if (iTube.DiaOrificeIn > 0. || iTube.DiaOrificeOut > 0.) {
					color = 210;
				}
				else if (iTube.PointIn == 0 || iTube.PointOut == 0 || Branches[iTube.NbBr].isFlowSet2zero) {
					color = 0;
				}
				else {
					color = 10;
					Value = iTube.dpdyn / iTube.Length;
					for (unsigned short ic = 1; ic <= 17; ic++) {
						//                  prot << "\n ic " << ic << " value " << (ValueMin + ic * ValueStep) << " dp/l " << iTube.dpdyn / iTube.Length;
						if (Value < (ValueMin + ic * ValueStep)) {
							color = (18 - ic) * TEN;
							//                     prot << " color " << color;
							break;
						}
					}
				}
				//            prot << "\n Itb " << iTube.Number << " MassVel " << iTube.MassVel << " color " << color << endl;
				// rainbow from blue (lowest) to red (highest)
				error = DXFWriteTube(outData, iTube, handle,
					LayerName, color, (!Branches[iTube.NbBr].isFlowSet2zero));
				if (error) return error;
			}
		}
		break;
	case ShowMode::ResistanceFactor:
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				if (!(iTube.DiaOrificeIn > 0. || iTube.DiaOrificeOut > 0.)) { // disregard orifices
					Value = iTube.dpdyn / iTube.Length / iTube.MassVel / iTube.MassVel;
					if (Value > ValueMax)ValueMax = Value;
					if (Value < ValueMin)ValueMin = Value;
				}
			}
		}
		ValueStep = (ValueMax - ValueMin) / 17.;
		//		prot << "\n\nvaluemin " << ValueMin << " valuemax " << ValueMax << " step " << ValueStep << endl;
		for (const auto& iTube : Tubes) {
			if (iTube.Dia < MinDrumDiameter) {
				Value = iTube.dpdyn / iTube.Length / iTube.MassVel / iTube.MassVel;
				content = to_string(iTube.Number);
				LayerName = "Tb " + content + " " + to_string(Value);
				if (iTube.DiaOrificeIn > 0. || iTube.DiaOrificeOut > 0.) {
					color = 210;
				}
				else if (Branches[iTube.NbBr].isFlowSet2zero) {
					color = 0;
				}
				else {
					color = 10;
					for (unsigned short ic = 1; ic <= 17; ic++) {
						//						prot << "\n number " << iTube.Number<<" ic " << ic << " value " << (ValueMin + ic * ValueStep) << " dp/l/msd^2 " << Value;
						if (Value < (ValueMin + ic * ValueStep)) {
							color = (18 - ic) * TEN;
							break;
						}
					}
					//					prot << " color " << color;
				}
				// rainbow from blue (lowest, color= 180 ) to red (highest, color = 10)
				error = DXFWriteTube(outData, iTube, handle,
					LayerName, color,(!Branches[iTube.NbBr].isFlowSet2zero) );
				if (error) return error;
			}
		}
		break;
		}
		outData << "  0\nENDSEC\n  0\nSECTION\n  2\nOBJECTS\n  0\nENDSEC\n  0\nEOF" << std::endl;
		return error;
	}
