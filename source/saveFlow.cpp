/*****************************************************************//**
 * \file   saveFlow.cpp
 * \brief saves the flow in branches of this iteration step to file
*
* \return int error code
 *
* \author Rainer_Jordan@<very, very warm>mail.com
* #### Licence
* Licensed under the European Union Public Licence (EUPL), Version 1.2 or later
* \date   September 2021
 *********************************************************************/


//#include "stdafx.h"
#undef MAINFUNCTION
#include "CommonHeader.h"

int saveFlow() {
	/* Local variables */
	string writefile;
	ofstream outData;

	writefile = PathFile + "Flow" + to_string(Base.iterg) + ".dat";

	outData.open(writefile.c_str());
	if (!outData.good()) {
		cout << "\n cannot open file " << writefile << endl;
		exit(1);
	}
	else {
		for (const auto& iBranch : Branches) {
			if (iBranch.NoChanges % 2) {
				outData << iBranch.Number << "    " << -iBranch.g << " dynLin " << iBranch.dPLinear << " dPConstant "
					<< iBranch.dPConstant << " gNew " << iBranch.gNew << " dpdyn " << iBranch.dPdyn << " dpstat " << iBranch.dPstat << "\n";
			}
			else {
				outData << iBranch.Number << "    " << iBranch.g << " dynLin " << iBranch.dPLinear << " dPConstant "
					<< iBranch.dPConstant << " gNew " << iBranch.gNew << " dpdyn " << iBranch.dPdyn << " dpstat " << iBranch.dPstat << "\n";
			}
		}
		outData.close();
	}//file opened correctly
	return (0);
}

