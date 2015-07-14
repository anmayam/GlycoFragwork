// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)



#pragma once
#using <mscorlib.dll>
#include "clsPeak.h"
#include "PeakProcessor/PeakProcessor.h"
#include "clsPeakProcessorParameters.h" 

namespace GlypID
{
	public __value enum enmProfileType {CENTROIDED = 0, PROFILE}; 

	namespace Peaks
	{
		public __gc class clsPeakProcessor
		{
			enmProfileType menmProfileType ; 
			clsPeakProcessorParameters *mobj_parameters ; 
			Engine::PeakProcessing::PeakProcessor __nogc *mobj_peak_processor ;
		public:
			clsPeakProcessor(void);
			~clsPeakProcessor(void);
			double GetBackgroundIntensity(float (&intensities) __gc []) ; 
			void SetPeakIntensityThreshold(double threshold);
			void DiscoverPeaks(float (&mzs) __gc [], float (&intensities) __gc [], GlypID::Peaks::clsPeak* (&peaks) __gc [], 
				float start_mz, float stop_mz, bool use_threshold) ;
			void SetOptions(clsPeakProcessorParameters *parameters) ;
			void InitializeUnprocessedData() ; 
			void GetClosestPeakMz(GlypID::Peaks::clsPeak* peak, float peak_mz) ; 


			__property enmProfileType get_ProfileType()
			{
				return  menmProfileType ; 
			}

			__property void set_ProfileType(enmProfileType type)
			{
				menmProfileType = type ; 
			}
		};
	}
}
