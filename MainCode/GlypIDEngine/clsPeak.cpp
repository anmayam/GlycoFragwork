// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#include ".\clspeak.h"

namespace GlypID
{
	namespace Peaks
	{
		clsPeak::clsPeak(void)
		{
		}

		clsPeak::~clsPeak(void)
		{
		}

		void clsPeak::Set(Engine::PeakProcessing::Peak &pk)
		{
			this->mdbl_mz = pk.mdbl_mz ; 
			this->mdbl_FWHM = pk.mdbl_FWHM ; 
			this->mdbl_intensity = pk.mdbl_intensity ; 
			this->mdbl_SN = pk.mdbl_SN ; 
			this->mint_data_index = pk.mint_data_index ;
			this->mint_peak_index = pk.mint_peak_index ; 
		}
	}
}
