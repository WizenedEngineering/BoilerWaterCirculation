#ifndef _PROPH2O_H
#define _PROPH2O_H

#include <math.h>
#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#define FABS(x) ((x) >= 0. ? (x) : -(x))
//#define PI 3.14159265
#define WATER 2
#define STEAM 1
#define FLUID 3

using namespace std;

/**
* @brief calculates properties water and steam
*
 * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
*/
class H2O {
public:
   /**
   * @brief calculates specific enthalpy for given temperature and pressure
   *
   * @param tK temperature [K]
   * @param pMPa absolute pressure [MPa]
   * @param medium indication of state 1:steam, 2:water, 3:supercritical fluid
   * @return double enthalpy [kJ/kg]
 * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
   */
   static double enth(double tK, double pMPa, int medium);
   /**
   * @brief calculates specific entropy for given temperature and pressure
   *
   * @param tK temperature [K]
   * @param pMPa absolute pressure [MPa]
   * @param medium indication of state 1:steam, 2:water, 3:supercritical fluid
   * @return double
  * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  */
	static double entropy(double tK, double pMPa, int medium);

   /**
   * @brief calculates specific heat capacity at constant pressure
   *
   * @param tK temperature [K]
   * @param pMPa absolute pressure [MPa]
   * @param medium indication of state 1:steam, 2:water, 3:supercritical fluid
   * @return double specific heat capacity [kJ/(kg K)]
  * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  */
	static double cp(double tK, double pMPa, int medium);

   /**
   * @brief calculates specific heat capacity at constant volume
   *
   * @param tK temperature [K]
   * @param pMPa absolute pressure [MPa]
   * @param medium indication of state 1:steam, 2:water, 3:supercritical fluid
   * @return double specific heat capacity [kJ/(kg K)]
  * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  */
	static double cv(double tK, double pMPa, int medium);

   /**
   * @brief calculates specific volume for given temperature and pressure
   *
   * @param tK temperature [K]
   * @param pMPa absolute pressure [MPa]
   * @param medium indication of state 1:steam, 2:water, 3:supercritical fluid
   * @return double specific volume [m3/kg]
  * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  */
	static double specVol(double tK, double pMPa, int medium);

   /**
   * @brief calculates saturation pressure for given temperature
   *
   * @param tSatK temperature [K]
   * @return double saturation pressure [MPa]
  * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  */
	static double satPress(double tSatK);

   /**
   * @brief calculates saturation temperature for given pressure
   *
   * @param pMPa absolute pressure [MPa]
   * @return double saturation temperature [K]
  * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  */
	static double satTemp(double pMPa);

   /**
   * @brief calculates dynamic viscosity
   *
   * @param tK temperature [K]
   * @param specVol specific volume [m3/kg]
   * @return double dynamic viscosity [Pa s]
  * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  */
	static double dynVisc(double tK, double specVol);

   /**
   * @brief calculates conductivity
   *
   * @param tK temperature [K]
   * @param pMPa absolute pressure [MPa]
   * @param medium indication of state 1:steam, 2:water, 3:supercritical fluid
   * @return double Conductivity [W/mK]
  * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  */
	static double conduct(double tK, double pMPa, int medium);

   /**
   * @brief calculation of temperature from given spec. enthalpy
   *
   * @param ha specific enthalpy [kJ/kg]
   * @param pMPa absolute pressure [MPa]
   * @return double temperature [K]
 * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
   */
	static double temp(double ha, double pMPa);

   /**
   * @brief calculates speed of sound
   *
   * @param tK temperature [K]
   * @param pMPa absolute pressure [MPa]
   * @param medium indication of state 1:steam, 2:water, 3:supercritical fluid
   * @return double speed of sound [m/s]
  * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  */
	static double SoS(double tK, double pMPa, int medium);

   /**
   * @brief calculates surface tension at saturation
   *
   * @param tSatK saturation temperature [K]
   * @return double surface tension [mN/m]
  * #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
  */
	static double sigma(double tSatK);

private:
	static double gamma1(double pNorm, double tNorm); /* gamma function in region 1(WATER T<350degC)*/
	static double gamma1Tau(double pNorm, double tNorm); /* derivative of gamma with respect to tau region 1 */
	static double gamma1Pi(double pNorm, double tNorm); /* derivative of gamma with respect to pi region 1 */
	static double gamma1TauTau(double pNorm, double tNorm); /*second  derivative of gamma with respect to  tau region 1 */
	static double gamma1PiPi(double pNorm, double tNorm); /*derivative of gamma with respect to pi region 1 */
	static double gamma1PiTau(double pNorm, double tNorm); /*derivative of gamma with respect to pi and tau region 1 */

	static double gamma2r(double pNorm, double tNorm); /* gamma function in region 2(STEAM)*/
	static double gamma2rTau(double pNorm, double tNorm); /* derivative of gamma with respect to tau region 2 */
	static double gamma2rPi(double pNorm, double tNorm); /* derivative of gamma with respect to pi region 2 */
	static double gamma2rTauTau(double pNorm, double tNorm); /*second  derivative of gamma with respect to tau region 2 */
	static double gamma2rPiTau(double pNorm, double tNorm); /* derivative of gamma with respect to pi and tau region 2 */
	static double gamma2rPiPi(double pNorm, double tNorm); /*second derivative of gamma with respect to pi region 2 */
	static double gamma2(double pNorm, double tNorm); /* gamma function in region 2(STEAM)*/
	static double gamma2Tau(double pNorm, double tNorm); /* derivative of gamma with respect to tau region 2 */
	static double gamma2Pi(double pNorm, double tNorm); /* derivative of gamma with respect to pi region 2 */
	static double gamma2TauTau(double pNorm, double tNorm); /*second derivative of gamma with respect to tau region 2 */
	static double gamma2PiTau(double pNorm, double tNorm); /*derivative of gamma with respect to pi and tau region 2 */
	static double gamma2PiPi(double pNorm, double tNorm); /*second derivative of gamma with respect to pi region 2 */
	static double phi3(double rhoNorm, double tNorm);
	static double phi3Delta(double rhoNorm, double tNorm);
	static double phi3Tau(double rhoNorm, double tNorm);
	static double phi3TauTau(double rhoNorm, double tNorm);
	static double phi3DeltaTau(double rhoNorm, double tNorm);
	static double phi3DeltaDelta(double rhoNorm, double tNorm);
	static double p3_H2O(double tC, double rho);
	static double approx3_WATER(double pBar, double TempC);
	static double approx3_STEAM(double pBar, double TempC);
};

//static H2O *h2o = new H2O();
//auto h2o = std::make_unique(H2O);
#endif
