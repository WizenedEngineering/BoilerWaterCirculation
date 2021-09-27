/* 
 * File:   common.h
 * \author rainer_jordan-at-<very, very warm>mail.com
 * #### Licence
 * Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
 * \date Created on November 13, 2016, 8:59 PM
 */

#ifndef COMMON_H
#define COMMON_H

#include <cstdlib>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>      // std::invalid_argument

// No.. number of (count)
// Nb.. specific number of tube, branch etc.
/*!
\def MINUS1
\brief a -1 for size_t variable 
  gives 0 when incremented
*/
#define MINUS1 0xffffffffffffffffULL
/*!
\def M_PI
\brief pi
*/
#define M_PI 3.14159265

static const size_t DRUM = 0;

/*!
* \class _tube 
* \brief class handling one tube
*
*/
class _tube {
public: //just to avoid using structure 
    size_t PointIn;        ///< number of point at tube start (inlet)
    size_t PointOut;       ///< number of point at tube end (outlet)
    double dOut;           ///< outside diameter of tube [mm]
    double Thickness;      ///< tube thickness [mm]
    double RadiusBend;     ///< radius of bend [mm]
    double q;              ///< heat absorption of all parallel tubes [kW]
    double NoParallel;     ///< number of parallel tubes
    double ksiAdd;         ///< resistance factor for additional resistance
    double DiaOrificeIn;   ///< diameter of orifice at tube inlet [m], if <=0 not available
    double DiaOrificeOut;  ///< diameter of orifice at tube outlet [m], if <=0 not available
    double FactHeat;       ///< factor for position of heating
    double EnthGiven;      ///< given enthalpy at tube inlet (difference to drum water enthalpy) [kJ/kg] (to allow for heated downcomers)
    double OCSStartAngle;  ///< start angle of arc (bend, elbow) [deg]
    double OCSEndAngle;    ///< end angle of arc (bend, elbow) [deg]
    double OCSCenterX;     ///< x-coordinate center point (OCS) of arc (bend, elbow) [mm]
    double OCSCenterY;     ///< y-coordinate center point (OCS) of arc (bend, elbow) [mm]
    double OCSCenterZ;     ///< z-coordinate center point (OCS) of arc (bend, elbow) [mm]
    double Nx;             ///< x-extrusion of arc (bend, elbow) [-]
    double Ny;             ///< y-extrusion of arc (bend, elbow) [-]
    double Nz;             ///< z-extrusion of arc (bend, elbow) [-]
    std::string Name;      ///< name of the tube
    std::string LayerInfo; ///< string containing the information of the layer
   // constructor

   _tube() {
      PointIn = MINUS1;   // number of point at tube start 
      PointOut = MINUS1;  // number of point at tube end [1]
      dOut = 0.;          // outside diameter in mm
      Thickness = 0.;     // thickness in mm
      RadiusBend = 0.;    // radius of bend in m 
      q = 0.;             // heat absorption of all parallel tubes in kW
      NoParallel = 1.;    // number of parallel tubes
      ksiAdd = 0.;        // resistance factor for additional resistance
      DiaOrificeIn = 0.;  // diameter of orifice at tube inlet [m], if <=0 not available
      DiaOrificeOut = 0.; // diameter of orifice at tube outlet [m], if <=0 not available
      FactHeat = 1.;      // factor for position of heating
      EnthGiven = 0.;     // given enthalpy
      OCSStartAngle = 0.; // start angle of arc (bend, elbow) [deg]
      OCSEndAngle = 0.;   // end angle of arc (bend, elbow) [deg]
      OCSCenterX = 0.;    // x-coordinate center point (OCS) of arc (bend, elbow) [mm]
      OCSCenterY = 0.;    // y-coordinate center point (OCS) of arc (bend, elbow) [mm]
      OCSCenterZ = 0.;    // z-coordinate center point (OCS) of arc (bend, elbow) [mm]
      Nx = 0.;            // x-extrusion of arc (bend, elbow) [-]
      Ny = 0.;            // y-extrusion of arc (bend, elbow) [-]
      Nz = 1.;            // z-extrusion of arc (bend, elbow) [-]
      Name = " ";         // name of the tube 
      LayerInfo = " ";
   }
};

/*!
* \class _point 
* \brief data related to one point
*
*/
class _point {
public:
    double xCoord; ///< x Coordinate of point in mm
    double yCoord; ///< y Coordinate of point in mm
    double zCoord; ///< z Coordinate of point in mm
    int NoTb;      ///< number of tubes connected to the point
   //constructor

   _point() {
      xCoord = 0.; // x Coordinate of point in mm
      yCoord = 0.; // y Coordinate of point in mm
      zCoord = 0.; // z Coordinate of point in mm 
      NoTb = 0; // number of tubes connected to the point
   }
};
/**
 * \brief converts OCS data of arcs to WCS data.
 *
 * \param OCSCenterX x-coordinate of center point in OCS
 * \param OCSCenterY y-coordinate of center point in OCS
 * \param OCSCenterZ z-coordinate of center point in OCS
 * \param Nx x-component of extrusion vector
 * \param Ny y-component of extrusion vector
 * \param Nz z-component of extrusion vector
 * \param OCSRadius radius of arc
 * \param OCSAngle angle of point in OCS
 * \return
 */
extern _point OCS2WCS(double OCSCenterX, double OCSCenterY, double OCSCenterZ,
   double Nx, double Ny, double Nz, double OCSRadius, double OCSAngle);

#endif /* COMMON_H */

