// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)



#pragma once
#using <mscorlib.dll>
#include "PeakProcessor/Peak.h"

namespace GlypID
{
	namespace Peaks
	{
		public __gc class clsPeak
		{
			public:
			//! mz of the peak.
			double mdbl_mz ; 
			//!  intensity of peak.
			double mdbl_intensity ;
			//! Signal to noise ratio
			double mdbl_SN ; 
			//! index in PeakData::mvect_peak_tops std::vector. 
			int mint_peak_index ;
			//! index in mzs, intensity vectors that were used to create the peaks in PeakProcessor::DiscoverPeaks.
			int mint_data_index ;
			//! Full width at half maximum for peak.
			double mdbl_FWHM ;
			clsPeak(void);
			~clsPeak(void);
			void Set(Engine::PeakProcessing::Peak &pk) ;
		};
	}
}
