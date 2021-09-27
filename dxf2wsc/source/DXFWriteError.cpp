/*!
 * \file DXFWriteError.cpp
 *
 * \author rainer_jordan-at-<very, very warm>mail.com
 * #### Licence
 * Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
 * \date   February 2021
 ******************************************************************** */

#include "common.h"

using namespace std;

extern int DXFWriteLine(ostream& outData, unsigned long long& handle, const string& layer, unsigned short color,
        double xStart, double yStart, double zStart, double xEnd, double yEnd, double zEnd);

/**
 * \fn DXFWriteError (ostream& outData, vector <size_t> &OrphanedPoints, vector <_point> &Points, vector<_point> &TubeErrors)
 * \brief writing lines between steam drum center and error point to .dxf file
 *
 * \param &outData outstream of .dxf file
 * \param &OrphanedPoints vector of points numbers with no or single tube connected
 * \param &Points vector of points
 * \param &TubeErrors vector of points indicating tubes with wrong tube data
 * \return int error code
 ******************************************************************** */

int DXFWriteError(ostream& outData,
        vector <size_t> &OrphanedPoints,
        vector <_point> &Points,
        vector<_point> &TubeErrors,
   vector<_point> &Intersect
) {
   int error = 0;
   //using new handles 
   unsigned long long handle = 100uLL;
   unsigned short color = 1;
   //! writing dxf minimal header 
   outData << "999\nDXF2WC\n  0\nSECTION\n  2\nCLASSES\n  0\nENDSEC\n  0\nSECTION\n  2\nENTITIES\n";
   //! writing lines to points that have one tube connected 
   if (!OrphanedPoints.empty()) {
      for (auto& iPt : OrphanedPoints) {
         error = DXFWriteLine(outData, handle, "0", color,
                 100000., 100000., 100000.,
                 Points[iPt].xCoord,
                 Points[iPt].yCoord,
                 Points[iPt].zCoord);
      }
   }
   //! writing lines to tubes with errors like too short or too small inside diameter
   color =6;
   if (!TubeErrors.empty()) {
      for (auto& Pt : TubeErrors) {
         error = DXFWriteLine(outData, handle, "0", color,
                 100000., 100000., 100000.,
                 Pt.xCoord,
                 Pt.yCoord,
                 Pt.zCoord);
      }
   }
   //! writing lines to intersections 
   color = 221;
   if (!Intersect.empty()) {
      for (auto& Pt : Intersect) {
         error = DXFWriteLine(outData, handle, "0", color,
            100000., 100000., 100000.,
            Pt.xCoord,
            Pt.yCoord,
            Pt.zCoord);
      }
   }
   //! writing dxf footer
   outData << "  0\nENDSEC\n  0\nSECTION\n  2\nOBJECTS\n  0\nENDSEC\n  0\nEOF\n";
   return error;
}
