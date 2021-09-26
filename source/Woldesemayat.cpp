//#include "stdafx.h"
#undef MAINFUNCTION
// It is assumed that the tube sections are small enough that the properties
// are mean value between inlet and outlet
#include "CommonHeader.h"

double _tube::Density_Woldesemayat(double x, double rhoW, double rhoS, double SurfTens, double PMPa, double& VoidFraction) {
   /* Local variables */
   double c0, vr, sinPhi, cosPhi, VelSG, VelSL, exp;

   sinPhi = HeightSection / LengthSection;
   cosPhi = sqrt(1. - sinPhi * sinPhi);
   VelSG = x * MassVel / rhoS;
   VelSL = (1. - x) * MassVel / rhoW;

   exp = pow(rhoS / rhoW, 0.1);
   c0 = VelSG / (VelSG + VelSL) * (1. + pow(VelSL / VelSG, exp));
   vr = 2.9 * pow(1.22 + 1.22 * sinPhi, 0.101325 / PMPa) * sqrt(sqrt(Dia * SurfTens * 9.80665 * (1. + cosPhi)*(rhoW - rhoS)) / rhoW);
   VoidFraction = VelSG / (c0 * (VelSG + VelSL) + vr);
   double VoidFractionHomogeneous = rhoW * x / (rhoW * x + rhoS * (1. - x)); //starting at homogeneous

   if (sinPhi > 0.) {
      /**
       * Upward flow: the range for void fraction is between homogeneous void fraction and x\n 
       * in some cases like low mass velocity or high x content the result can lay outside this range \n
       * that means we are outside the validity of this formula\n
       * for high x content a homogeneous solution is sensible (the flow pattern can change to mist flow and mist flow is more like homogeneous ) \n 
       * for low mass velocity the steam velocity can be significantly higher than water velocity. The void fraction approaching x\n 
       */
      VoidFraction = fmax(x, VoidFraction);
      VoidFraction = fmin(VoidFractionHomogeneous, VoidFraction);
   } else {
      /**
       * Downward flow: the range for void fraction is between homogeneous void fraction and 1.\n 
       * in some cases like low mass velocity or high x content the result can lay outside this range \n
       * that means we are outside the validity of this formula\n
       * for high x content a homogeneous solution is sensible (the flow pattern can change to mist flow and mist flow is more like homogeneous ) \n 
       * for low mass velocity there is a high chance for finely dispersed bubbles that again lead to no velocity difference of the phases -> homogeneous\n 
       */
      if (VoidFraction < VoidFractionHomogeneous || VoidFraction > 0.9999999) {
         VoidFraction = VoidFractionHomogeneous;
      }
   }
   return rhoW * (1. - VoidFraction) + rhoS * VoidFraction;
} /* Density_Woldesemayat */

