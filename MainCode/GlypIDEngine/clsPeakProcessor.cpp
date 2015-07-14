// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#include "clspeakprocessor.h"
#include "GlypIDEngineUtils.h"
#using <mscorlib.dll>
#include <vector>

namespace GlypID
{
	namespace Peaks
	{
		clsPeakProcessor::clsPeakProcessor(void)
		{
			mobj_peak_processor = new Engine::PeakProcessing::PeakProcessor () ; 
			mobj_parameters = new clsPeakProcessorParameters() ; 
			mobj_peak_processor->SetOptions(mobj_parameters->get_SignalToNoiseThreshold(), 0, false, (Engine::PeakProcessing::PEAK_FIT_TYPE) mobj_parameters->get_PeakFitType()) ; 
			menmProfileType = PROFILE ; 
		}

		clsPeakProcessor::~clsPeakProcessor(void)
		{
			if (mobj_peak_processor != NULL)
			{
				delete mobj_peak_processor ; 
				mobj_peak_processor = NULL ; 
			}
		}

		double clsPeakProcessor::GetBackgroundIntensity(float (&intensities) __gc [])
		{
			double thres = GlypID::Utils::GetAverage(intensities, FLT_MAX) ; 
			thres = GlypID::Utils::GetAverage(intensities, (float)(5*thres)) ;
			return thres ; 
		}		

		void clsPeakProcessor::SetPeakIntensityThreshold(double threshold)
		{
			mobj_peak_processor->SetPeakIntensityThreshold(threshold) ; 
		}

		void clsPeakProcessor::InitializeUnprocessedData()
		{
			mobj_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;
		}

		void clsPeakProcessor::GetClosestPeakMz(GlypID::Peaks::clsPeak *peak, float peak_mz)
		{
			Engine::PeakProcessing::Peak closest_peak ; 
			double mz = mobj_peak_processor->GetClosestPeakMz((double) peak_mz, closest_peak) ; 
			if (mz >0)
				peak->Set(closest_peak) ; 
		}

	
		void clsPeakProcessor::DiscoverPeaks(float (&mzs) __gc [], float (&intensities) __gc [], 
			GlypID::Peaks::clsPeak* (&peaks) __gc [], float start_mz, float stop_mz, bool use_threshold) 
		{
			std::vector<double> vectMzs ;
			std::vector<double> vectIntensities ;
			int numPoints = mzs->Length ; 
			for (int ptNum = 0 ; ptNum < numPoints ; ptNum++)
			{
				vectMzs.push_back((double)mzs[ptNum]) ; 
				vectIntensities.push_back((double)intensities[ptNum]) ; 
			}
			
			if (use_threshold)
			{
				double backgroundIntensity = GetBackgroundIntensity(intensities) ; 
				mobj_peak_processor->SetPeakIntensityThreshold(backgroundIntensity*mobj_parameters->get_PeakBackgroundRatio()) ; 
			}
			
		
			if (menmProfileType == PROFILE)
				mobj_peak_processor->SetPeaksProfileType(true) ; 
			else
				mobj_peak_processor->SetPeaksProfileType(false) ; 

			int numPeaks = this->mobj_peak_processor->DiscoverPeaks(&vectMzs, &vectIntensities, start_mz, stop_mz) ; 

			

			peaks = new clsPeak* __gc [numPeaks] ; 
			for (int pkNum = 0 ; pkNum < numPeaks ; pkNum++)
			{
				peaks[pkNum] = new clsPeak() ; 
				Engine::PeakProcessing::Peak pk ; 
				this->mobj_peak_processor->mobj_peak_data->GetPeak(pkNum, pk) ; 
				peaks[pkNum]->Set(pk) ; 
			}
		}

		void clsPeakProcessor::SetOptions(clsPeakProcessorParameters *parameters)
		{
			mobj_parameters = parameters ; 
			// the minimum intensity is not set till the actual data is available in DiscoverPeaks
			mobj_peak_processor->SetOptions(mobj_parameters->get_SignalToNoiseThreshold(), 0, mobj_parameters->get_ThresholdedData(), (Engine::PeakProcessing::PEAK_FIT_TYPE)mobj_parameters->get_PeakFitType()) ; 
		}
	}
}
