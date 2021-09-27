/*****************************************************************//**
 * \file   cross.cpp
 * \brief  checks for intersection between start and end point of lines/tubes
 * 
 * \author rainer_jordan-at-<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  *********************************************************************/
#include "common.h"

extern size_t mTb; // max index of tubes (max number of tubes - 1)
extern size_t mPt; // max index of Points (number of points - 1)  

using namespace std;

/**
 * \fn check_intersection(vector <_tube>& Tubes, vector <_point>& Points, vector<_point>& Intersect )
 * \param Tubes vector of all tube(s)
 * \param Points vector of all point(s)
 * \param Intersect vector of intersection points, for error indication in dxf file
 * \return 
 */
int check_intersection(vector <_tube>& Tubes, vector <_point>& Points, vector<_point>& Intersect ) {
   _point Intersection;
   double t,u;

//   size_t maxTubes = Tubes.size();
   for (size_t i = 0; i <=  mTb - 1; i++) {
      for (size_t j = i + 1; j <= mTb; j++) {
         if (Tubes[i].PointIn == Tubes[j].PointIn ||
            Tubes[i].PointOut == Tubes[j].PointIn ||
            Tubes[i].PointIn == Tubes[j].PointOut ||
            Tubes[i].PointOut == Tubes[j].PointOut) {
            continue;
         }
         double x1 = Points[Tubes[i].PointIn].xCoord;
         double x2 = Points[Tubes[i].PointOut].xCoord;
         double x3 = Points[Tubes[j].PointIn].xCoord;
         double x4 = Points[Tubes[j].PointOut].xCoord;
         double y1 = Points[Tubes[i].PointIn].yCoord;
         double y2 = Points[Tubes[i].PointOut].yCoord;
         double y3 = Points[Tubes[j].PointIn].yCoord;
         double y4 = Points[Tubes[j].PointOut].yCoord;
         double z1 = Points[Tubes[i].PointIn].zCoord;
         double z2 = Points[Tubes[i].PointOut].zCoord;
         double z3 = Points[Tubes[j].PointIn].zCoord;
         double z4 = Points[Tubes[j].PointOut].zCoord;

         double diff21x = x2 - x1;
         double diff21y = y2 - y1;
         double diff21z = z2 - z1;

         double diff43x = x4 - x3;
         double diff43y = y4 - y3;
         double diff43z = z4 - z3;

         double diff31x = x3 - x1;
         double diff31y = y3 - y1;
         double diff31z = z3 - z1;

         double cx = diff21y * diff43z - diff21z * diff43y;
         double cy = diff21z * diff43x - diff21x * diff43z;
         double cz = diff21x * diff43y - diff21y * diff43x;
         /// distance between 2 tubes is a fraction expression. It can only be 0 if numerator is 0 and denominator is not 0
         double numerator = diff31x * cx + diff31y * cy + diff31z * cz;
         double denominator = sqrt(cx * cx + cy * cy + cz * cz);
         if (fabs(denominator) < 1e-9) {
            continue;
         }
         else {
            double dist = numerator / denominator;
            if (fabs(dist) < 1e-3) { /// distance between straight lines is 0 -> intersection or co-linear
               bool isCx = (fabs(cx) > fabs(cy) && fabs(cx) > fabs(cz));
               bool isCy = (fabs(cy) > fabs(cx) && fabs(cy) > fabs(cz));
               bool isCz = (fabs(cz) > fabs(cx) && fabs(cz) > fabs(cy)); 
               if (isCy) {
                  t = (diff31x * diff43z - diff43x * diff31z) / -cy;
                  u = (diff31x * diff21z - diff21x * diff31z) / -cy;
                  if (t >= 0. && t <= 1. && u >= 0. && u <= 1.) {
                     Intersection.xCoord = x1 + t * diff21x;
                     Intersection.yCoord = y1 + t * diff21y;
                     Intersection.zCoord = z1 + t * diff21z;
                     Intersect.push_back(Intersection);
                  }
               }
               else if (isCx) {
                  t = (diff31z * diff43y - diff43z * diff31y) / -cx;
                  u = (diff31z * diff21y - diff21z * diff31y) / -cx;
                  if (t >= 0. && t <= 1. && u >= 0. && u <= 1.) {
                     Intersection.xCoord = x1 + t * diff21x;
                     Intersection.yCoord = y1 + t * diff21y;
                     Intersection.zCoord = z1 + t * diff21z;
                     Intersect.push_back(Intersection);
                  }
               }
               else if(isCz) {
                  t = (diff31x * diff43y - diff43x * diff31y) / cz;
                  u = (diff31x * diff21y - diff21x * diff31y) / cz;
                  if (t >= 0. && t <= 1. && u >= 0. && u <= 1.) {
                     Intersection.xCoord = x1 + t * diff21x;
                     Intersection.yCoord = y1 + t * diff21y;
                     Intersection.zCoord = z1 + t * diff21z;
                     Intersect.push_back(Intersection);
                  }
               }
            }
         }
      }
   }
   return 0;
}

