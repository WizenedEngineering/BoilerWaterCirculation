/*! \file CommonHeader.h
   \brief contains some class definitions, global functions and global data 
*/

#ifndef _NCOMMON_H
#define _NCOMMON_H
//#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <system_error>
//#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include<Eigen/IterativeLinearSolvers>
//#include <Eigen/SparseCholesky>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include "propH2O.h"
using namespace std;
using namespace Eigen;

#ifdef MAINFUNCTION
#define EXTERN
#else
#define EXTERN extern
#endif

/*!
\def MINUS1
\brief a -1 for size_t variable 
  gives 0 when incremented
*/
#define MINUS1 0xffffffffffffffffULL

static const size_t DRUM = 0;
static const double MinDrumDiameter = 0.5;
//using Matrix = std::vector<std::vector<double>>;
//using dVector = std::vector<double>;
//using iVector = std::vector<int>;
using stVector = std::vector<size_t>;

/// No.. number of (count), how many
/// Nb.. specific number of tube, branch etc.

/*!
\var size_t mBr; 
\brief max index of branches (number of branches - 1)
*/
EXTERN size_t mBr; // max index of branches (number of branches - 1)

/*!
\var size_t mTb; 
\brief max index of tubes (number of tubes - 1)
*/
EXTERN size_t mTb; // max index of tubes (number of tubes - 1)
/*!
\var size_t mPt; 
\brief max index of Points (number of points - 1)
*/
EXTERN size_t mPt; // max index of Points (number of points - 1)

/*!
\var size_t mNd; 
\brief max index of nodes (number of nodes - 1)
*/
EXTERN size_t mNd; // max index of nodes (number of nodes - 1)

/*!
\var string PathFile
\brief full path to working directory as well as project name
*/
EXTERN string PathFile;

/*!
\var ofstream prot;
\brief file handle for protocol file 

used to trace data during calculation 
*/
EXTERN ofstream prot;

/*!
* \enum FlowPattern
* \brief enum FlowPattern gives the flow patterns 
*
* Based on Steiner 
*/
enum class FlowPattern {
	undetermined, ///< flow pattern is undetermined
	stratified,   ///< stratified flow: to be avoided
	wavy,         ///< wavy flow: to be avoided
	bubble,       ///< bubble flow
	plugSlug,     ///< plug or slug flow
	mist,         ///< mist flow: to be avoided
	annular       ///< annular flow 
};

/*!
* \enum KindOfBranch
* \brief enum KindOfBranch stores the different kind of branches like heated, downcomer etc.
*
* just to get rid of index numbers
*/
enum class KindOfBranch
{
   Undefined,   ///< not of the other, maybe not in any initial path
   Heated,      ///< branch is heated
   Downcomer,   ///< downcomer (not heated and passing downcomer level)
   Down2Heated, ///< connection between downcomer and heated branch (f.i. lower header)
   Heated2Drum, ///< connection between heated branch and drum (f.i. upper header, overflow)
   Drum2Down    ///< connection between drum and downcomer inlet
};

/*!
* \enum ShowMode
* \brief What should be shown by dxf file 
* 
* just to get rid of index numbers
*/
enum class ShowMode
{
	TubesPoints,     ///< showing tubes and points 
   BranchesNodes,   ///< showing tubes and coloring according kind of branch
	Flow,            ///< showing tubes and coloring according kind of branch, Flow [kg/s] in layer
   Velocity,        ///< showing tubes, coloring according velocity
   Arrows,          ///< showing arrows of flow direction only
	SafetyFactor,    ///< showing tubes, coloring according SafetyFactor
	SteamQuality,    ///< showing tubes, coloring according steam quality
	VoidFraction,    ///< showing tubes, coloring according void fraction
   MassVel,         ///< showing tubes, coloring according mass velocity
	FlowPattern,     ///< showing tubes, coloring according flow pattern(Steiner)
	DPdynPerLength,  ///< showing tubes, coloring according dpdyn / length(to find bottle neck)
	ResistanceFactor ///< showing tubes, coloring according dpdyn / length / mass velocity^2(to find bottle neck)
};

/*!
* \enum TFlow
* \brief indicates the flow patten in Tee
*
*/
enum class TFlow
{
   StraightSeparationInletOff,       ///< flow in straight is separated: inlet of off branch
   StraightSeparationInletStraight,  ///< flow in straight is separated: inlet of straight branch
   StraightSeparationOutletStraight, ///< flow in straight is separated: outlet of straight branch
   OffSeparationOutletOff,           ///< flow in off is separated to 2 straight: outlet of off branch
   OffSeparationInletStraight,       ///< flow in off is separated to 2 straight: inlet of straight branch
   StraightUnionOutletOff,           ///< flow in straight is union of straight and off: outlet of off branch
   StraightUnionOutletStraight,      ///< flow in straight is union of straight and off: outlet of straight branch
   StraightUnionInletStraight,       ///< flow in straight is union of straight and off: inlet of straight branch
   OffUnionOutletStraight,           ///< flow in off is union from 2 straight: outlet of straight branch
   OffUnionInletOff                  ///< flow in off is union from 2 straight: inlet of off branch
};

/*!
 * \enum TeeOrientation
 * \brief enum TeeOrient gives the orientation of a Tee
 *
 * is needed for steam distribution in "enth"
 */
enum class TeeOrientation
{
	NoTee,
	Undetermined,          ///< no distinct orientation
   StraightHorOffHor,     ///< Straight horizontal, Off horizontal
   StraightHorOffVerUp,   ///< Straight horizontal, Off vertical up
   StraightHorOffVerDown, ///< Straight horizontal, Off vertical down
   StraightVerOffHor,     ///< Straight vertical, Off horizontal
   OffHor,                ///< Straight every else, Off horizontal
   OffVerUp,              ///< Straight every else, Off vertical up
   OffVerDown             ///< Straight every else, Off vertical down
};

/*!
 * \enum TbRegion
 * \brief enum TbRegion gives which part of a tube should be calculated
 *
 */
enum class TbRegion
{
   Tube,        ///< plain tube, resistance factor to be calculated  
   SinglePlace, ///< single place, resistance determined by given resistance factor
   Orifice,     ///< orifice, resistance factor to be calculated 
   Bend         ///< bend or elbow, resistance factor to be calculated  
};

/*!
 * \enum USArrangement
 * \brief enum USArrangement indicates if successive bends are in U or S arrangement
 *
 */
enum class USArrangement
{
   No, ///< no U or S 
   S,  ///< S arrangement
   U,  ///< U arrangement
};

/*!
* \class _base
* \brief Base holds some basic values
*
*/
class _base {
public:
    double tol;            ///< tolerance
    double Rough;          ///< tube roughness in mm (valid for all tubes)
    double DowncomerLevel; ///< level that downcomers pass, determined by highest and lowest point and LevelFact (input)
    double CRStart;        ///< initial guess of circulation ratio for the first branch flow estimation initFlow() 
    bool showMesh;         ///< switch for printing of protocol showing mesh data
    bool showEnth;         ///< switch for printing of protocol showing enthalpy in nodes data
    bool showDPTube;       ///< switch for printing of protocol showing pressure difference in tubes
    bool showDPBranch;     ///< switch for printing of protocol showing pressure difference in branches
    bool showNodePressure; ///< switch for printing of protocol showing nodal pressure
    bool showFlow;         ///< switch for printing of protocol showing flow
    bool showReverse;      ///< switch for printing of protocol showing details on reversing branches
    bool showMeshDetail;   ///< switch for printing of protocol showing more details on mesh data and initial flow
    bool showReverseDetail;///< switch for printing of protocol showing more details on reversing branches
    bool showDPTubeDetail; ///< switch for printing of protocol showing more details on pressure difference in tubes
    int iterg;             ///< counter for flow iterations
    int maxit;             ///< maximum flow iterations
    char Method;           ///< calculation method (J for Jirous/Jirous, W for Welzer/Welzer, R for Rouhani/Becker, C for Chexal,Lellouche/VDIHA,G for Woldesemayat/VDIHA
};
/*!
\var _base Base
\brief single invocation of structure _base 

Base is holding some base values
*/
EXTERN _base Base;

/*!
* \struct _drum
* \brief data of the steam drum
*
*/
struct _drum {
    double pMPa;      ///< drum pressure in MPa(absolute)
    double enthW;     ///< spec. enthalpy of saturated water in drum [kJ/kg]
    double enthS;     ///< spec. enthalpy of saturated water in drum [kJ/kg]
    double enthEvap;  ///< evaporation enthalpy at drum pressure[kJ/kg]
    double rhoW;      ///< density of saturated water in drum [kg/m3]
    double LevelW;    ///< Water level in drum, difference to drum centerline [mm]
    // + level is above center, - level is below center
    double dpDyn;     ///< pressure drop of drum internals [Pa]
    double qSum;      ///< sum of all heat in the boiler in kW
    double FlowSteam; ///< steam flow from drum [kg/s]
};
/*
\var  _drum Drum;
\brief single invocation of structure -drum
*/
EXTERN _drum Drum;

/*!
* \class _point 
* \brief data related to one point
*
*/
class _point {
public:
    size_t Number;
    size_t NbNd;   ///< number of node in this point, preset to -1
    stVector NbTb; ///<Tube number of tube connected to this point
    size_t NoTb;   ///< number of tubes connected to the point
    double xCoord; ///< x Coordinate of point in mm
    double yCoord; ///< y Coordinate of point in mm
    double zCoord; ///< z Coordinate of point in mm
    //   constructor

    _point() {
        Number = MINUS1;
        NbNd = MINUS1; // number of node in this point, preset to -1
        NoTb = 0;      // number of tubes connected to the point
        NbTb.clear(); //Tube number of tube connected to this point
        xCoord = -1.; // x Coordinate of point in mm
        yCoord = -1.; // y Coordinate of point in mm
        zCoord = -1.; // z Coordinate of point in mm
    }
};

/*!
\var vector <_point> Points;
\brief vector holding all points (of class _point)
*/
EXTERN vector <_point> Points;

/*!
* \class _tube 
* \brief class handling one tube
*
*/
class _tube {
public:
    size_t Number;         ///< index number in Tubes-vector
    size_t PointIn;        ///< number of point at tube start (inlet)
    size_t PointOut;       ///< number of point at tube end (outlet)
    size_t NbBr;           ///< number of the branch this tube belongs to; preset to max. possible number
    double Dia;            ///< inside diameter of tube [m]
    double area;           ///< cross section area [m2]
    double Length;         ///< length [m]
    double Height;         ///< height [m]
    double q;              ///< heat absorption of all parallel tubes [kW]
    double NoParallel;     ///< number of parallel tubes
    double ksiAdd;         ///< resistance factor for additional resistance
    double ksiIn;          ///< resistance factor at inlet, needed for determination of energy dissipation
    double ksiOut;         ///< resistance factor at outlet, needed for determination of energy dissipation
    double ksiTube;        ///< resistance factor of tube, needed for determination of energy dissipation
    double beta;           ///< angle between tubes in deg
    USArrangement UorS;    ///< 3 sequential tubes in U or S arrangement
    double RadiusBend;     ///< radius of bend [m]
    double DiaOrificeIn;   ///< diameter of orifice at tube inlet [m], if <=0 not available
    double DiaOrificeOut;  ///< diameter of orifice at tube outlet [m], if <=0 not available
    int NoSections;        ///< number of sections
    double LengthSection;  ///< Length of section [m]
    double HeightSection;  ///< Height of section [m]
    double HeatSection;    ///< heat absorption of one section of one tube [kW]
    double pPaIn;          ///< absolute pressure at tube inlet [Pa]
    double EnthIn;         ///< enthalpy at tube inlet [kJ/kg]
    double EnthInGiven;    ///< given enthalpy at tube inlet [kJ/kg], below drum saturation enthalpy to indicate distinct (even heated) downcomers 
    double xIn;            ///< Steam quality at tube inlet [-]
    double VoidFractionIn; ///< void fraction at inlet (volume of steam / total volume) [-]
    double velIn;          ///< mixture velocity at tube inlet [m/s]
    double pPaOut;         ///< absolute pressure at tube outlet [Pa]
    double EnthOut;        ///< enthalpy at tube outlet [kJ/kg]
    double xOut;           ///< SteamQuality at tube outlet [-]
    double VoidFractionOut;///< void fraction at outlet  (volume of steam / total volume) [-]
    double velOut;         ///< mixture velocity at tube outlet [m/s]
    double FactHeat;       ///< factor for position of heating
    double HeatFlux;       ///< heat flux in tube [kW/m2]
    double dpIn;           ///< pressure drop at inlet [Pa]
    double dpOut;          ///< pressure drop at outlet [Pa]
    double dpdyn;          ///< dynamic pressure difference [Pa]
    double dpstat;         ///< static pressure difference [Pa]
    double rhoIn;          ///< density at inlet [kg/m3]
    double rhoOut;         ///< density at outlet [kg/m3]
    double rhoMean;        ///< mean density [kg/m3]
    double Flow;           ///< Flow in one single tube [kg/s]
    double MassVel;        ///< mass velocity in tube [kg/m2 s]
    //   double Reynolds; //Reynolds number
    double velWater;       ///< velocity if density is water density ("superficial" water velocity)[m/s]
    double SafetyFactor;   ///< safety factor against flow separation or overheating (should be higher than 1)
    double SafetyKorneev;   ///< safety factor against overheating according to Korneev
    double SafetyTaitel_Dukler;  ///< safety factor against flow separation according to Taitel-Dukler
    double SafetySteiner;   ///< safety factor against flow separation according to Steiner
    double SafetyKonkov;   ///< safety factor against overheating according to Kon'kov
    double SafetyDoroshchuk;   ///< safety factor against overheating according to Doroshchuk
    double SafetyKatto_Ohno;   ///< safety factor against overheating according to Katto-Ohno
    double SafetyGroeneveld;   ///< safety factor against overheating according to Groeneveld
    FlowPattern steiner;   ///< flow pattern according Steiner
    double OCSStartAngle;  ///< start angle of arc (bend) in OCS, needed for .dxf file
    double OCSEndAngle;    ///< start angle of arc (bend) in OCS, needed for .dxf file
    double OCSCenterX;     ///< x-coordinate of arc (bend) center point in OCS, needed for .dxf output file
    double OCSCenterY;     ///< y-coordinate of arc (bend) center point in OCS, needed for .dxf output file
    double OCSCenterZ;     ///< z-coordinate of arc (bend) center point in OCS, needed for .dxf output file
    double Nx;             ///< x-component of arc (bend) extrusion vector, needed for .dxf output file
    double Ny;             ///< y-component of arc (bend) extrusion vector, needed for .dxf output file
    double Nz;             ///< z-component of arc (bend) extrusion vector, needed for .dxf output file
    std::string Name;      ///< name of the tube

    /// constructor
    _tube() {
        Number = MINUS1;
        FactHeat = 1.;
        PointIn = MINUS1;
        PointOut = MINUS1;
        UorS = USArrangement::No;
        NbBr = MINUS1;
        dpIn = 0.;
        dpOut = 0.;
        DiaOrificeIn = 0.;
        DiaOrificeOut = 0.;
        NoSections = 1; // number of sections
        LengthSection = 0.; // Length of section [m]
        HeightSection = 0.; // Height of section [m]
        HeatSection = 0.; // heat absorption of one section of one tube [kW]
        VoidFractionOut = 0.;
        Dia = 0.;
        area = 0.;
        Length = 0.;
        Height = 0.;
        q = 0.;
        NoParallel = 0.;
        ksiAdd = 0.;
        ksiIn=0.; // resistance factor at inlet
        ksiOut=0.; //resistance factor at outlet
        ksiTube=0.; // resistance factor of tube
        RadiusBend = 0.;
        pPaIn = 0.;
        EnthIn = 0.;
        xIn = 0.;
        VoidFractionIn = 0.;
        velIn = 0.;
        pPaOut = 0.;
        EnthOut = 0.;
        xOut = 0.;
        velOut = 0.;
        beta = 0.;
        HeatFlux = 0.;
        dpdyn = 0.;
        dpstat = 0.;
        rhoIn = 0.;
        rhoOut = 0.;
        rhoMean = 0.;
        Flow = 0.;
        MassVel = 0.;
        velWater = 0.;
        SafetyFactor = -10.;
        SafetyKorneev = -10.;
        SafetyTaitel_Dukler = -10.;  
        SafetySteiner = -10.;   
        SafetyKonkov = -10.;   
        SafetyDoroshchuk = -10.;   
        SafetyKatto_Ohno = -10.;   
        SafetyGroeneveld = -10.;  
        steiner = FlowPattern::undetermined;
        OCSStartAngle = 0.;
        OCSEndAngle = 0.;
        OCSCenterX = 0.;
        OCSCenterY = 0.;
        OCSCenterZ = 0.;
        Nx = 0.;
        Ny = 0.;
        Nz = 0.;
        Name = " ";
    }
    //   functions
 /*!
*  \brief determines spatial angle between 2 tubes *this* and *other*
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
* \param [in] other reference to other tube
* \return spatial angle in radians
*/
    double Angle3d(const _tube& other);

/*!
 * \brief Determines the number of sections of this tube
 *
 * The calculation of pressure difference (pressure drop) in a long tube can lead to\n
 * quite different pressures at inlet and outlet and therefore different properties of water/steam.\n
 * To keep the properties in a small range the tube is split into different sections where (hopefully) the properties are almost constant.  
 * The number of sections has to be determined for the lowest flow in each flow iteration step because the steam quality depends on flow\n 
 * if number of sections from previous steps is higher than calculated one -> the higher number is used
 * \param [in] g flow in the tube [kg/s]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
     */
    void DetermineNoSections(double g);

/*!
* \brief calculates the pressure difference in a tube
*
* the tube has different regions:\n
* Inlet\n
* inlet orifice\n
* Bend or elbow\n
* straight tube divided into sections\n
* outlet orifice\n
* outlet\n
* those regions are handled sequentially
* 
* 
* 
* \param [in] ResFactIn resistance factor at inlet [-]
* \param [in] ResFactOut resistance factor at outlet [-]
* \param [in] g Flow in all parallel tubes [kg/s]
* \return int error code
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*/
    int dpTube ( double ResFactIn, double ResFactOut, double g );

/*!
* \brief calculation of pressure difference in tube section
*
* function determines inlet condition of the section and branches to dpSinglePhase(water or steam) or dpTwoPhase\n
* within the section there is no switch from single phase to two phase\n
* even if after some length the conditions call for 2 phase the whole section is calculated as single phase. 
* The next section is calculated as two phase.\n
* The same principle applies to the transition between 2 phase and steam.\n 
* Even if the tubes are cooled or 2 phase downward flow this result in a conservative approach
* \param [in] region single place (inlet, outlet = 1,  orifice=2, bend = 3) or tube = 0
* \param [in] FrictCoeff friction coefficient [-]
*         if single place: friction coefficient of this particular place
*         if tube: additional friction coefficient
* \param [in] enthInSect inlet enthalpy [kJ/kg]
* \param [in] enthOutSect outlet enthalpy [kJ/kg]
* \param [in] LengthOrifice length (thickness) of orifice [m]
* \param [out] PhaseChange bool indication if there was a phase change in the section
* \param [in,out] tSatOutSect saturation temperature at outlet [K]
* \param [in,out] volWOutSect spec. volume of saturated water at outlet  [m3/kg]
* \param [in,out] volSOutSect spec. volume of saturated steam at outlet [m3/kg]
* \param [in,out] SurfTensOutSect surface tension at outlet [N/m]
* \param [in,out] dynVisSOutSect dynamic viscosity of saturated steam at outlet [Pa s]
* \param [in,out] dynVisWOutSect dynamic viscosity of saturated water at outlet [Pa s]
* \param [in,out] dynVisOutSect dynamic viscosity at outlet [Pa s], single phase
* \param [in,out] enthWSatOutSect spec. enthalpy of saturated water at outlet [kJ/kg]
* \param [in,out] enthSSatOutSect spec. enthalpy of saturated steam at outlet [kJ/kg]
* \param [in,out] pPaOutSect absolute pressure at outlet [Pa]
* \param [in,out] dpDynSect dynamic pressure drop [Pa]
* \param [in,out] dpStatSect  static pressure difference [Pa]
* \param [in,out] xInSect steam quality at inlet
* \param [in,out] xOutSect steam quality at outlet
* \param [in,out] rhoInSect density at inlet [kg/m3]
* \param [in,out] rhoOutSect  density at outlet [kg/m3]
*         if no previous section, f.i. tube inlet, it has to be 0.
* \param [in, out] VoidInSect  void fraction at inlet
* \param [in, out] VoidOutSect void fraction at outlet
* \param [out] rhoMeanSect mean velocity [m/s].
* \param [out] VelSect mean velocity [m/s]
* \return int error code: 0 all ok, 1 phase change, 2 out of turbulent range 3 pressure drop too high, outlet pressure below atmosphere
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
   */
    int dPSection ( // returns error code: 0 all ok, 1 phase change, 2 out of turbulent range
        //3 pressure drop too high, outlet pressure below atmosphere
        // no phase change in tube or single place
        TbRegion region, // indicating the region of the tube: single place (inlet, outlet etc.) is default but FrictCoeffAdd must be given, orifice, bend, tube
        double FrictCoeff,      //friction coefficient [-]
        // if single place: friction coeff of this particular place
        // if tube: additional friction coeff
        //           double l, // length of tube section[m]
        //           double h, // height of tube section(+ upward flow. - downward flow) [m]
        double enthInSect,       // inlet enthalpy [kJ/kg]
        double enthOutSect,      // outlet enthalpy [kJ/kg]
        double LengthOrifice,    // length (thickness) of orifice [m]
        bool &PhaseChange,       // a phase change occurred
        double& tSatOutSect,     //saturation temperature at outlet [K]
        double& volWOutSect,     //spec. volume of saturated water at outlet  [m3/kg]
        double& volSOutSect,     //spec. volume of saturated steam at outlet [m3/kg]
        double& SurfTensOutSect, //surface tension at outlet [N/m]
        double& dynVisSOutSect,  // dynamic viscosity of saturated steam at outlet [Pa s]
        double& dynVisWOutSect,  //dynamic viscosity of saturated water at outlet [Pa s]
        double& dynVisOutSect,   //dynamic viscosity at outlet [Pa s], single phase
        double& enthWSatOutSect, //spec. enthalpy of saturated water at outlet [kJ/kg]
        double& enthSSatOutSect, //spec. enthalpy of saturated steam at outlet [kJ/kg]
        double& pPaOutSect,      // absolute pressure at outlet [Pa]
        double& dpDynSect,       // dynamic pressure drop [Pa]
        double& dpStatSect,      // static pressure difference [Pa]
        double& xInSect,         //steam quality at inlet
        double& xOutSect,        //steam quality at outlet
        double& rhoInSect,       // density at inlet [kg/m3]
        double& rhoOutSect,      // density at outlet [kg/m3]
        // if no previous section, it has to be 0.
        double& VoidInSect , // void fraction at inlet
        double& VoidOutSect, // void fraction at outlet 
        double& rhoMeanSect,     // mean density [kg/m3]
        double& VelSect );       // mean velocity [m/s]

/*!
* \brief  dPSinglePhase calculates the static and dynamic pressure difference of a tube section for water or steam phase
*
* \param [in] region indicating the region of the tube: single place (inlet, outlet etc.) is default but FrictCoeffAdd must be given, orifice, bend, tube
* \param [in] FrictCoeff friction coefficient [-]
* \param [in] enthInSect inlet enthalpy [kJ/kg]
* \param [in] enthOutSect outlet enthalpy [kJ/kg]
* \param [in] LengthOrifice Length of orifice [m]
* \param [in] index index for water or steam
* \param [out] PhaseChange  =0, if outlet above saturation = 1
* \param [in,out] dynViscOutSect dyn. viscosity at inlet [Pa s]..
* \param [in,out] pPaOutSect absolute pressure at outlet [Pa]
* \param [in,out] dpDynSect dynamic pressure drop [Pa]
* \param [in,out] dpStatSect static pressure difference [Pa]
* \param[in,out]  xOutSect steam quality at outlet
* \param [in,out] rhoInSect density at inlet [kg/m3]
* \param [in,out] rhoOutSect density at outlet [kg/m3]
* \param [out] rhoMeanSect mean density [kg/m3]
* \param [out] VelSect mean velocity [m/s]
* \return int error code 
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
   */
    int dPSinglePhase ( // returns error code: 0 all ok,
        //3 pressure drop too high, outlet pressure below atmosphere
        // no phase change in tube or single place
        TbRegion region,        // single place (inlet, outlet, valve :1, orifice, valve) or tube = 0
        double FrictCoeff,      //friction coefficient [-]
        // if single place: friction coeff of this particular place
        // if tube: additional friction coeff not covered by deflection, orifice etc. 
        double enthInSect,      // inlet enthalpy [kJ/kg]
        double enthOutSect,     // outlet enthalpy [kJ/kg]
        double LengthOrifice,   // Length of orifice [m]
        int index,              // index for water or steam
        bool& PhaseChange,      // =0, if outlet above saturation = 1
        double& dynViscOutSect, // dyn. viscosity at inlet [Pa s]
        double& pPaOutSect,     // absolute pressure at outlet [Pa]
        double& dpDynSect,      // dynamic pressure drop [Pa]
        double& dpStatSect,     // static pressure difference [Pa]
        double& xOutSect,       // steam quality at outlet
        double& rhoInSect,      // density at inlet [kg/m3]
        double& rhoOutSect,     // density at outlet [kg/m3]
        // if no previous section, it has to be 0.
        double& rhoMeanSect,    // mean density [kg/m3]
        double& VelSect );      // mean velocity [m/s]

/*!
* \brief dPTwoPhase calculates the static and dynamic pressure difference of a tube section in 2 phase flow
*
* It is assumed that the tube sections are small enough that the properties are mean value between inlet and outlet
* 
* \param  [in] region indicating the region of the tube: single place (inlet, outlet etc.) is default but FrictCoeffAdd must be given, orifice, bend, tube
* \param  [in] FrictCoeffAdd friction coefficient [-]if single place: friction coeff of this particular place, if tube: additional friction coeff
* \param  [in] enthIn inlet enthalpy [kJ/kg]
* \param  [in] enthOut outlet enthalpy [kJ/kg]
* \param  [in] LengthOrifice length (thickness) of orifice [m]
* \param  [out] PhaseChange did a phase change happen 2 phase to single phase
* \param  [in,out] tSatOut saturation temperature at outlet [degC]
* \param  [in,out] volWOut spec. volume of saturated water at outlet  [m3/kg]
* \param  [in,out] volSOut spec. volume of saturated steam at outlet [m3/kg]
* \param  [in,out] SurfTensOut surface tension at outlet [N/m]
* \param  [in,out] dynVisSOut dynamic viscosity of saturated steam at outlet [Pa s]
* \param  [in,out] dynVisWOut dynamic viscosity of saturated water at outlet [Pa s]
* \param  [in,out] enthWSatOut spec. enthalpy of saturated water at outlet [kJ/kg]
* \param  [in,out] enthSSatOut spec. enthalpy of saturated steam at outlet [kJ/kg]
* \param  [in,out] pPaOutSect absolute pressure at outlet [Pa]
* \param  [in,out] dpDynSect dynamic pressure drop [Pa]
* \param  [in,out] dpStatSect static pressure difference [Pa]
* \param  [in,out] xInSect steam quality at inlet
* \param  [in,out] xOutSect steam quality at outlet
* \param  [in,out] rhoInSect density at inlet [kg/m3] 
* \param  [in,out] rhoOutSect density at outlet [kg/m3]  if no previous section, it has to be 0.
* \param  [out] VoidInSect void fraction at inlet
* \param  [out] VoidOutSect void fraction at outlet 
* \param  [out] rhoMeanSect mean density [kg/m3]
* \param  [out] VelSect mean velocity [m/s]
* \return error code
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
   */
    int dPTwoPhase (
        TbRegion region, // indicating the region of the tube: single place (inlet, outlet etc.) is default but FrictCoeffAdd must be given, orifice, bend, tube
        double FrictCoeffAdd, //friction coefficient [-]
        // if single place: friction coeff of this particular place
        // if tube: additional friction coeff
        //           double l, // length of tube section[m]
        //           double h, // height of tube section(+ upward flow. - downward flow) [m]
        double enthIn, // inlet enthalpy [kJ/kg]
        double enthOut, // outlet enthalpy [kJ/kg]
        double LengthOrifice, // length (thickness) of orifice [m]
        bool& PhaseChange,
        double& tSatOut, //saturation temperature at outlet [degC]
        double& volWOut, //spec. volume of saturated water at outlet  [m3/kg]
        double& volSOut, //spec. volume of saturated steam at outlet [m3/kg]
        double& SurfTensOut, //surface tension at outlet [N/m]
        double& dynVisSOut, // dynamic viscosity of saturated steam at outlet [Pa s]
        double& dynVisWOut, //dynamic viscosity of saturated water at outlet [Pa s]
        double& enthWSatOut, //spec. enthalpy of saturated water at outlet [kJ/kg]
        double& enthSSatOut, //spec. enthalpy of saturated steam at outlet [kJ/kg]
        double& pPaOutSect, // absolute pressure at outlet [Pa]
        double& dpDynSect, // dynamic pressure drop [Pa]
        double& dpStatSect, // static pressure difference [Pa]
        double& xInSect, //steam quality at inlet
        double& xOutSect, //steam quality at outlet
        double& rhoInSect, // density at inlet [kg/m3]
        double& rhoOutSect, // density at outlet [kg/m3]
        // if no previous section, it has to be 0.
        double& VoidInSect, // void fraction at inlet
        double& VoidOutSect, // void fraction at outlet 
        double& rhoMeanSect, // mean density [kg/m3]
        double& VelSect ); // mean velocity [m/s]

/*!
* \brief calculates the mixture density according to Chexal-Lellouche
*
* Source: EPRI Report TR-106326, Void Fraction Technology for Design \n
* and Analysis, December 1996. 
* 
* \param [in] x steam quality [-]
* \param [in] pMPa pressure [MPa]
* \param [in] volW spec.volume of saturated water [m3/kg]
* \param [in] volS spec.volume of saturated steam [m3/kg]
* \param [in] dynVisW dynamic viscosity of saturated water [Pa s]
* \param [in] dynVisS dynamic viscosity of saturated steam [Pa s]
* \param [in] SurfTens surface tension [N/m]
* \param [out] VoidFraction void fraction [-]
* \return mixture density [kg/m3]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
     */
    double Density_chexal(double x, double pMPa, double volW,
        double volS, double dynVisW, double dynVisS, double SurfTens, double& VoidFraction);

/*!
* \brief calculates the dynamic pressure difference according to Heat Atlas
* 
* Source:VDI Heat Atlas, Second Edition L2.2
*
* \param [in] pPaInSect pressure at inlet [MPa]
* \param [in] pPaOutSect pressure at outlet [MPa]
* \param [in] volSSatIn spec.volume of saturated steam at inlet [m3/kg]
* \param [in] volSSatOut spec.volume of saturated steam at outlet [m3/kg]
* \param [in] volWSatIn spec.volume of saturated water at inlet [m3/kg]
* \param [in] volWSatOut spec.volume of saturated water at outlet [m3/kg]
* \param [in] dynVisWSatIn dynamic viscosity of saturated water at inlet [Pa s]
* \param [in] dynVisWSatOut dynamic viscosity of saturated water at outlet [Pa s]
* \param [in] dynVisSSatIn dynamic viscosity of saturated steam at inlet [Pa s]
* \param [in] dynVisSSatOut dynamic viscosity of saturated water at outlet [Pa s]
* \param [in] xMean steam quality [-]
* \return dynamic pressure difference [Pa]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
     */
    double dpdyn_HeatAtlas(double pPaInSect, double pPaOutSect, double volSSatIn,
        double volSSatOut, double volWSatIn, double volWSatOut, double dynVisWSatIn,
        double dynVisWSatOut, double dynVisSSatIn, double dynVisSSatOut, double xMean
        );

/*!
* \brief calculates the mixture density according to Rouhani
*
* mixture density for upward flow according Rouhani/Axelson and Gomez for downward flow
* 
* Source: Rouhani, S. Z., and Axelsson, E. Calculation of Void Volume Fraction in the Subcooled and Quality Boiling Regions \n
* International Journal of Heat and Mass Transfer, vol. 13, no. 2, pp. 383–393, 1970\n
* and \n
* Gomez, L. E, Shoham, O., Schmidt, Z., Chokshi, R. N., and Northug, T. Unified Mechanistic Model for Steady State Two Phase Flow:\n
* Horizontal to Vertical Upward Flow. Society of Petroleum Engineers Journal, vol. 5, pp. 339–350, 2000.
* 
* \param [in] x steam quality [-]
* \param [in] rhoW density of saturated water [kg/m3]
* \param [in] rhoS density of saturated steam [kg/m3]
* \param [in] SurfTens surface tension at outlet [N/m]
* \param [out] VoidFraction void fraction [-]
* \return mixture density [kg/m3]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
    */
    double Density_Rouhani(double x, double rhoW, double rhoS, double SurfTens, double& VoidFraction );
    double Gomez(double x, double rhoW, double rhoS, double SurfTens, double VoidFractionInput);

/*!
* \brief calculates the dynamic pressure difference according to Becker
* 
* Source: Becker K.M., Hernborg O, and Bode M. "An experimental study of Pressure Gradients for
* Flow of Boiling Water in Vertical Round Ducts (Part 4)", AB-86, 1962
* \param [in] zeta friction coefficient [-]
* \param [in] xIn steam quality at inlet [-]
* \param [in] xOut steam quality at outlet [-]
* \param [in] pPaIn pressure at inlet [MPa]
* \param [in] pPaOut pressure at outlet [MPa]
* \param [in] volWSatIn spec.volume of saturated water at inlet [m3/kg]
* \param [in] volWSatOut spec.volume of saturated water at outlet [m3/kg]
* \return dynamic pressure difference [Pa]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
     */
    double dpdyn_Becker(double zeta, double xIn, double xOut, double pPaIn, double pPaOut,
        double volWSatIn, double volWSatOut);

/*!
* \brief calculates the dynamic pressure difference according to Jirous
* 
* Source:F.Jirous, Analytische Methode der Berechnung des Naturumlaufs
* bei Dampferzeugern, VGB Kraftwerkstechnik 58, Heft 5
*
* \param [in] zeta friction coefficient [-]
* \param [in] pMPa pressure [MPa]
* \param [in] volSIn spec.volume of saturated steam at inlet [m3/kg]
* \param [in] volSOut spec.volume of saturated steam at outlet [m3/kg]
* \param [in] volWIn spec.volume of saturated water at inlet [m3/kg]
* \param [in] volWOut spec.volume of saturated water at outlet [m3/kg]
* \param [in] x steam quality [-]
* \return dynamic pressure difference [Pa]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
    */
    double dpdyn_Jirous(double zeta, double pMPa, double volSIn,double volSOut, double volWIn, double volWOut, double x);

/*!
* \brief calculates the mixture density according to Jirous
*
* Source:F.Jirous, Analytische Methode der Berechnung des Naturumlaufs
* bei Dampferzeugern, VGB Kraftwerkstechnik 58, Heft 5
*
* \param [in] x steam quality [-]
* \param [in] volW spec.volume of saturated water [m3/kg]
* \param [in] volS spec.volume of saturated steam [m3/kg]
* \param [out] VoidFraction void fraction [-]
* \return mixture density [kg/m3]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
    */
    double Density_Jirous(double x, double volW, double volS, double& VoidFraction);

/*!
* \brief calculates the dynamic pressure difference according to Welzer
* 
* Source:Herbert Welzer, Probleme der Zweiphasenstroemung in Siederohren von Verdampfersystemen, Dissertation, 1969
*
* \param [in] zeta friction coefficient [-]
* \param [in] volW spec.volume of saturated water [m3/kg]
* \param [in] volS spec.volume of saturated steam [m3/kg]
* \param [in] enthW spec. enthalpy of saturated water [m3/kg]
* \param [in] enthS spec. enthalpy of saturated steam [m3/kg]
* \param [in] visW dynamic viscosity of saturated water [Pa s]
* \param [in] xInSect steam quality [-]
* \param [in] xOutSect steam quality [-]
* \return dynamic pressure difference [Pa]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
    */
    double dpdyn_Welzer(double zeta, double volW, double volS, double enthW, double enthS,
       double visW, double xInSect, double xOutSect);

/*!
* \brief calculates the mixture density according to Welzer
*
* Source:Herbert Welzer, Probleme der Zweiphasenstroemung in Siederohren von Verdampfersystemen, Dissertation, 1969
*
* \param [in] x steam quality [-]
* \param [in] volW spec.volume of saturated water [m3/kg]
* \param [in] volS spec.volume of saturated steam [m3/kg]
* \param [out] VoidFraction void fraction [-]
* \return mixture density [kg/m3]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
    */
    double Density_Welzer(double x, double volW, double volS, double& VoidFraction);

/*!
* \brief calculates the mixture density according to Woldesemayat
* 
* Source : Woldesemayat, M. A., and Ghajar, A. J. Comparison of Void Fraction Correlations for Different Flow Patterns\n
* in Horizontal and Upward Inclined Pipes. International Journal of Multiphase Flow, vol. 33, no. 4, pp. 347–370, 2007.
*
* \param [in] x steam quality [-]
* \param [in] rhoW density of saturated water [m3/kg]
* \param [in] rhoS density of saturated steam [m3/kg]
* \param [in] SurfTens surface tension at outlet [N/m]
* \param [in] pMPa pressure [MPa]
* \param [out] VoidFraction void fraction [-]
* \return mixture density [kg/m3]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
    */
    double Density_Woldesemayat(double x, double rhoW, double rhoS, double SurfTens, double pMPa, double& VoidFraction);

/*!
* \brief friction factor according Moody's diagram
* 
* Source: I don't recall, where I found the formula. It is in **good** agreement with the diagram and no iteration needed.  
*
* \param [in] Reynolds Reynolds number
* \param [in] RelRough relative Roughness (rough/dia) [-]
* \return friction factor
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
    */
    double FrictFact ( double Reynolds, double RelRough );

/*!
* \brief calculation of resistance factor for bends, elbows and sharp edged deflection
*
* Source: Idel'chik, Handbook of hydraulic resistance, Coefficients of local resistance and of friction
* 
* \param [in] Reynolds Reynolds number
* \param [in] AngleDeg Bend angle [degrees] 
* \return resistance factor
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
   */
    double ksiBend ( double Reynolds, double AngleDeg );

/*!
* \brief calculates resistance factor of orifice
*
* Source: Brandt, F., Dampferzeuger: Kesselsysteme, Energiebilanz, Stroemungstechnik. FDBR Fachbuchreihe Band 3, Vulkan-Verlag Essen, 1992  
*
* \param [in] LengthOrifice length (thickness) of orifice
* \param [in] visc dynamic viscosity of fluid [Pa s]
* \return resistance factor
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
    */
    double ksiOrifice ( double LengthOrifice, double visc );

/*!
* \brief function to determine safety factor against different criteria
*
* VDI Heat Atlas divides the range of mass velociy and pressure into 4 regions and recommends a method for each region\n
* Those methods are:\n
* Korneev for horizontal/inclined heated tubes\n
* Taitel-Dukler and Steiner for flow separation (only critical in heated tubes)\n
* Kon'kov for critical steam quality\n
* Doroshchuk for critical steam quality\n
* Katto-Ohno for critical steam quality\n
* Groeneveld for critical heat flux\n
* Kon'kov has following limitation: Dia < 33mm, MassVel > 200 kg/m2s; outside the limitation it will be disregarded\n
* Katto-Ohno is only valid for xIn <= 0\n
*
* the minimum safety factor is used
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*/
    void Criteria();

/*!
* \brief critical heat flux according table from Groeneveld
* 
* Source:  Groeneveld DC, Leung LKH, Kirillov PL, Bobkov VP, Smogalev IP, Vinogradov VN, Huang XC, Royer E (1996) \n
* The 1995 look-up table for criticalheat flux in tubes. Nucl Eng Des 163:1–23
*
* \param [in] pMPa Pressure [MPa]
* \param [in] rhoWSat Density of water at saturation [kg/m3]
* \param [in] rhoSSat Density of steam at saturation [kg/m3]
* \return critical heat flux [kW/m2]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
   */
    double chftable ( // return critical heat flux [kW/m2]
        double pMPa, // pressure [MPa]
        double rhoWSat, // Density of water at saturation [kg/m3]
        double rhoSSat ); // Density of steam at saturation [kg/m3]



/*!
* \brief determines the flow pattern in horizontal tubes
*
* Determines the flow pattern in horizontal tubes
* according the flow chart by Steiner
* 
* Source: VDI Heat Atlas, Second Edition, H3.1
*
* \param [in] pPa Pressure [Pa]
* \param [in] quality steam mass content, steam quality [-]
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*/
    double Steiner ( double pPa, double quality );
 
/*!
 * \brief determines the safety factor against wavy/stratified flow in horizontal tubes
 *
 * Determines the safety factor against wavy/stratified flow in horizontal tubes
 * according the criterion by Taitel- Dukler
 *
 * Source:  A Model for Predicting Flow Regime Transitions in Horizontal and Near Horizontal Gas-Liquid Flow, Taitel,Y,Dukler,A.E. , AIChE Journal (Vol.22, No 1)
 *
 * \param [in] pPa Pressure [Pa]
 * \param [in] quality steam mass content, steam quality [-]
 * \author Rainer_Jordan@<very, very warm>mail.com
 * #### Licence
 * Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
 * \date   January 2022
 */
    double Taitel_Dukler(double pPa, double quality);

/*!
* \brief reverses flow direction in tube
*
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*/
    void reverseDirection();

};

/*!
\var vector <_tube> Tubes;
\brief vector holding all tubes (of class _tube)
*/
EXTERN vector <_tube> Tubes;

/*!
* \class _branch 
* \brief class contains data related to one branch
*
*   Approximation for characteristic curve (pressure difference)
*    deltaP = dPLinear * g + dPConstant
*
* 
*/
class _branch {
public:
    size_t Number;           ///< index number of the branch within Branch-vector
    size_t NbPtIn;           ///< number of point at branch inlet
    size_t NbPtOut;          ///< number of point at branch outlet
    size_t NbNdIn;           ///< number of node at branch inlet
    size_t NbNdOut;          ///< number of node at branch outlet
    stVector NbTbInBr;       ///< tube numbers in branch
    size_t mTbInBr;          ///< max index# of tubes in branch (max number of tubes in a branch -1)
    KindOfBranch Kind;       ///< kind of branch
    bool isHorizontalHeated; ///< indicates if branch is horizontal and heated (needs special handling)
    double deltaH;           ///< height (elevation) difference between branch outlet and inlet[m]
    double enthIn;           ///< spec. enthalpy at inlet [kJ/kg]
    double enthOut;          ///< spec. enthalpy at outlet [kJ/kg]
    double qSum;             ///< sum of absorbed heat in kW
    double xIn;              ///< steam quality at inlet
    double g;                ///< Flow in kg/s
    double gNew;             ///< new flow
    double gPrev1;           ///< flow from last iteration step
    double gPrev2;           ///< flow from 2nd last iteration step
    double gPrev3;           ///< flow from 3rd last iteration step
    double gSteamIn;         ///< steam flow at inlet in kg/s
    double minArea;          ///< minimum cross section area  (to determine max. mass velocity)in Branch
    int neg;                 ///< counter how often negative flow was calculated
    int NoChanges;           ///< counter of flow direction changes
    bool isFlowSet2zero;     ///< determines whether the program has set the flow f.i. too many flow reversals -> 0.
    double dPdyn;            ///< total pressure difference in branch in Pa
    double dPstat;           ///< static pressure difference in branch in Pa
    double dPIn;             ///< pressure difference at branch inlet in Pa (by header or T-piece)
    double dPOut;            ///< pressure difference at branch outlet in Pa (by header or T-piece)
    double dPLinear;         ///< factor for pressure difference, characteristic curve (linear factor)
    double dPConstant;       ///< factor for pressure difference, characteristic curve (constant part)

    ///constructor
    _branch() {
        Number = MINUS1; // index number of the branch within Branch-vector
        NbPtIn = MINUS1; // number of point at branch inlet
        NbPtOut = MINUS1; // number of point at branch outlet
        NbNdIn = MINUS1; // number of node at branch inlet
        NbNdOut = MINUS1; // number of node at branch outlet
        isHorizontalHeated = false; // indicates if branch is horizontal and heated (needs special handling)
        deltaH = 0.; //height difference between branch inlet and outlet [m]
        NbTbInBr.clear(); // tube numbers in branch
        mTbInBr = MINUS1; // max index# of tubes in branch (max number of tubes in a branch -1)
        Kind = KindOfBranch::Undefined; // kind of branch
        enthIn = 0.; // spec. enthalpy at inlet [kJ/kg]
        enthOut = 0.; //spec. enthalpy at ourlet [kJ/kg]
        qSum = 0.; // sum of absorbed heat in kW
        xIn = 0.; // steam quality at inlet
        g = 0.; // Flow in kg/s
        gNew = 0.; // new flow
        gPrev1 = 0.; // flow from last iteration step
        gPrev2 = 0.; // flow from last iteration step
        gPrev3 = 0. ; // flow from last iteration step
        gSteamIn = 0.; // steam flow at inlet in kg/s
        neg = 0; // counter how often negative flow was calculated
        NoChanges = 0; // counter of flow direction changes
        isFlowSet2zero = false;
        minArea = 1e30;
        dPdyn = 0.; // total pressure difference in branch in Pa
        dPstat = 0.; // static pressure difference in branch in Pa
        dPIn = 0.; // pressure difference at branch inlet in Pa (by header or T-piece)
        dPOut = 0.; // pressure difference at branch outlet in Pa (by header or T-piece)
        dPLinear = 0.; // factor for dyn. pressure difference (linear factor)
        dPConstant = 0.; // factor for dyn. pressure difference (constant part)
        //			IsDirectionSet = false;
    }
    // functions
/*! 
* \brief calculation of pressure difference in branch
*
*     Calculation of pressure difference in branch
*     and coefficients of approximate characteristic curve
*     dp = dPLinear * g + dPConstant
* 
* the characteristic curve (pressure difference over flow) is approximated in vicinity of the flow in this step by a secant\n
* the pressure difference (dPdyn and dPstat) is calculated 3 times to get 3 data points for the secant linear approximation  
* \return int error code
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*/
    int dpBranch();

/*!
* \file   KsiTee.cpp
* \fn double KsiTee ( TFlow TCase, double Flowz, double Flowa, double Areaz, double Areaa, double Fillet, double angle, double dhyd )
* \brief Calculation for the resistance factor of tube inlet/outlet connected to a Tee
*
* It is also the inlet or outlet of a branch
* \param [in] TCase: 10 possible cases
* \param [in] Flowz: main flow [kg/s]
* \param [in] Flowa: flow in branch-off [kg/s]
* \param [in] Areaz: cross section area of main flow branch [m2]
* \param [in] Areaa: cross section area of branch-off [m2]
* \param [in] Fillet: fillet radius at branch-off tube inlet [m]
* \param [in] angle: angle between straight and branch-off [-]
* \param [in] dhyd: branch-off diameter [m]
* \return resistance factor
*
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*/
    double KsiTee(TFlow TCase, double Flowz, double Flowa, double Areaz, double Areaa,
       double Fillet, double angle, double dhyd);

/*!
* \brief reversing all parameters depending on flow direction in branch
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*
* \return int error code
*/
    int reverseDirection();
};

/*!
\var vector <_branch> Branches;
\brief vector holding all branches (of class _branch)
*/
EXTERN vector <_branch> Branches;

/*!
* \class _node
* \brief data related to one node
*
*/
class _node {
public:
    size_t Number;           ///< index number in Node-vector
    size_t NbPt;             ///< point number of this node
    stVector NbBr;           ///< number of branches in this node
    size_t mBrInNd;          ///< max index of branches in this node
    stVector NbBrArrive;     ///< number of branches ending in this node
    size_t mBrArrive;        ///< max index of arriving branches in this node
    stVector NbBrLeave;      ///< number of branches starting in this node
    size_t mBrLeave;         ///< max index of leaving branches per node
    bool isVisited;          ///< on the way to drum/downcomer this node is already visited
    bool IsT;                ///< Whether Node is a Tee piece
	 size_t NbTbTStraight[2]; ///< Tube numbers that are straight (inline) at a Tee piece
	 size_t NbTbTOff;         ///< tube number that branches off at a Tee piece (only one branch)
	 size_t NbBrTStraight[2]; ///< Branch numbers that are straight (inline) at a Tee piece
    size_t NbBrTOff;         ///< Branch number that branches off at a Tee piece (only one branch)
	 TeeOrientation TOrientation;    ///< it gives the orientation of Tee piece
            //NoTee,
            // Undetermined,           no distinct orientation
            // StraightHorOffHor,      Straight horizontal, Off horizontal
            // StraightHorOffVerUp,    Straight horizontal, Off vertical up
            // StraightHorOffVerDown,  Straight horizontal, Off vertical down
            // StraightVerOffHor,      Straight vertical, Off horizontal
            // OffHor,                 Straight every else, Off horizontal
            // OffVerUp,               Straight every else, Off vertical up
            // OffVerDown              Straight every else, Off vertical down
            // horizontal: +- 10 deg to horizontal, vertical: above 45 deg to horizontal
    double Elev;             ///< elevation of point in m (drum is set to 100m)
    double pNode;            ///< nodal pressure in Pa; to avoid problems with floating point precision
    //              it is the difference to drum pressure, i.e. pNode[0 = drum] = 0.
    //              For flow calculation only the differences between 2 Nodes are relevant
    double pPrev;            ///< nodal pressure from previous iteration step [Pa]
    double gSteam;           ///< steam flow in node in kg/s
    double gSum;             ///< sum of flow in node in kg/s (should be 0),
    // in mesh set-up it is used as sum of arriving or leaving flow for the first guess of flow in branches
    double enth;             ///< enthalpy of all arriving flows in node [kJ/kg]
    double gSumArrive;       ///< Sum of arriving flows [kg/s] used in enth and mesh
    double gSumLeave;        ///< Sum of leaving flows [kg/s] used in enth and mesh

    /// constructor
    _node() {
        Number = MINUS1;
        NbPt = MINUS1; // point number of this node
        mBrInNd = MINUS1; // max index of branches in this node
        mBrArrive = MINUS1; // max index of arriving branches in this node
        mBrLeave = MINUS1; // max index of leaving branches per node
        isVisited = false; //
        IsT = false; //Whether Node is a T-piece
        NbBrTStraight[0] = MINUS1; // Branch numbers that are straight (inline) at a T-piece
        NbBrTStraight[1] = MINUS1; // Branch numbers that are straight (inline) at a T-piece
        NbBrTOff = MINUS1; // Branch number that branches off at a T-piece (only one branch)
		  NbTbTStraight[0] = MINUS1; // Branch numbers that are straight (inline) at a T-piece
		  NbTbTStraight[1] = MINUS1; // Branch numbers that are straight (inline) at a T-piece
		  NbTbTOff = MINUS1; // Branch number that branches off at a T-piece (only one branch)
		  TOrientation = TeeOrientation::NoTee;
        gSum = 0.;
        gSumArrive = 0.;
        gSumLeave = 0.;
        gSteam = 0.;
        Elev = 0.; // elevation of point in m (drum is set to 100m)
        NbBr.clear(); // number of branches in this node
        NbBrArrive.clear(); // number of branches ending in this node
        NbBrLeave.clear(); // number of branches starting in this node
        pNode = 0.; // nodal pressure in Pa; to avoid problems with floating point precision
        //              it is the difference to drum pressure, i.e. pNode[0 = drum] = 0.
        //              For flow calculation only the differences between 2 Nodes are relevant
        pPrev = 0.; // nodal pressure from previous iteration step [Pa]
        enth = 0.; // enthalpy of all arriving flows in node [kJ/kg]
    }
 
/**
* \fn _node::isTee(_tube& firstTube, _tube& secondTube, _tube& thirdTube)
* \brief checks if this node is a Tee and sets the orientation of the Tee
*
* All 3 tubes are straight (no bends)\n
* 2 tubes should be in line (straight) and have the same diameter\n 
* the remaining tube should have equal or smaller diameter, the angle to the 2 other tubes should be between 45 deg and 135 deg  
* \param [in] firstTube reference to first tube
* \param [in] secondTube reference to second tube
* \param [in] thirdTube reference to third tube
* \return bool 
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*/
    bool isTee(_tube& firstTube, _tube& secondTube, _tube& thirdTube);
};
/*!
\var vector <_node> Nodes;
\brief vector holding all nodes (of class _node)
*/
EXTERN vector <_node> Nodes;

///functions
int main(int argc, char** argv);

int CalcPressureNodes(SparseMatrix<double>& S, /* System sparse matrix */
   VectorXd& B, /* right hand side vector............*/
   VectorXd& X /* solution vector ...................*/
   ); // error code as return value

extern int CalcEnthalpyNodes ( void );

extern int Mesh ( );

/*****************************************************************//**
* \file   StartResistance.cpp
* \brief calculates a resistance of branches to be used as "weight" for path finding
*
* for the determination of shortest (or split) path the branches have to have a "weight"\n
* this weight grows the more the branch is already flown\n
*
* in this case it is a pseudo resistance factor that is only needed in initFlow function
*
* it corresponds to the dynamic pressure drop at g = 1 kg/s
*
* the branch class member dPConstant is used for this purpose
*
* \return int error code
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/
extern int StartResistance();


extern bool LES(SparseMatrix<double>& S, /* System matrix as n*? sized vector.*/
   VectorXd& B, /* right hand side vector............*/
   VectorXd& X /* solution vector ...................*/
); // error code as return value

extern bool SingleEquation();

extern void initFlow ();

extern int readData ();

extern bool Step ( void );

/*****************************************************************//**
* \brief saves the flow in branches of this iteration step to file
*
* if the flow direction is opposite to the direction set in "initFlow" the flow in this branch is negative
* @return int error code
*
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/
extern int saveFlow ();

/*****************************************************************//**
* \brief calculates total energy dissipation
*
*Energy dissipation or energy loss in system try to minimize\n
*if there are 2 solutions the one with the lower energy dissipation should be used
*@return int error code
*
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/
extern double TotalDissipation();

/*****************************************************************//**
* \brief saves the calculation results as data file and different .dxf files
*
*  @return int error code
*
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
*********************************************************************/
extern int SaveResults();

extern int Print2dxf (ShowMode mode );

extern int DXFWrite ( ostream& outData,
                      ShowMode mode );
extern void processError(size_t iPt, const string& text);

#endif
