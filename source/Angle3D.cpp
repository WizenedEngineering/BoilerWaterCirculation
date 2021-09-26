//#include "stdafx.h"
#include "CommonHeader.h"
/* ----------------------------------------------------------------------- */
/*     function determines spatial angle between *this tube and other tube */
/* ----------------------------------------------------------------------- */

double _tube::Angle3d(const _tube& other) {
	double arg;
	arg = (Points[this->PointOut].xCoord - Points[this->PointIn].xCoord) / this->Length / 1e3 *
		((Points[other.PointOut].xCoord - Points[other.PointIn].xCoord) / other.Length / 1e3) +
		(Points[this->PointOut].yCoord - Points[this->PointIn].yCoord) / this->Length / 1e3 *
		((Points[other.PointOut].yCoord - Points[other.PointIn].yCoord) / other.Length / 1e3) +
		(Points[this->PointOut].zCoord - Points[this->PointIn].zCoord)	/ this->Length / 1e3 *
		((Points[other.PointOut].zCoord - Points[other.PointIn].zCoord) / other.Length / 1e3);
	if (fabs(arg) > 1.01) {
		cout << "\n error calculating angle between tube " << this->Number << " and " << other.Number;
		exit(1);
	}
	if (arg > 1.) arg = 1.;
	else if (arg < -1.) arg = -1.;

	return acos(arg);
} /* spatial angle in radians */

