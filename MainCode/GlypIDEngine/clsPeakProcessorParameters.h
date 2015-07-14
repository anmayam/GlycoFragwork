// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#pragma once
#using <mscorlib.dll>
#using <System.Xml.dll>
using namespace System::Xml ; 

namespace GlypID
{
	namespace Peaks
	{
		//! enumeration for type of fit. 
		public __value enum PEAK_FIT_TYPE { APEX = 0, QUADRATIC, LORENTZIAN } ; 

		public __gc class clsPeakProcessorParameters: public System::ICloneable
		{
			double mdbl_SNThreshold ; 
			double mdbl_PeakBackgroundRatio ; 
			bool mbln_thresholded_data ; 
			PEAK_FIT_TYPE menm_FitType ; 
		public:
			clsPeakProcessorParameters(void);
			clsPeakProcessorParameters(double sn, double peak_bg_ratio, bool thresholded_data, PEAK_FIT_TYPE fit_type);
			~clsPeakProcessorParameters(void);
			void LoadPeakParameters(XmlReader *rdr) ; 
			void SavePeakParameters(System::Xml::XmlTextWriter *xwriter) ; 

			virtual Object* Clone()
			{
				clsPeakProcessorParameters *new_params = new clsPeakProcessorParameters(mdbl_SNThreshold,
					mdbl_PeakBackgroundRatio, mbln_thresholded_data, menm_FitType) ; 
				return new_params ; 
			}

			__property bool get_ThresholdedData()
			{
				return mbln_thresholded_data ; 
			}
			__property void set_ThresholdedData(bool value)
			{
				mbln_thresholded_data = value ; 
			}

			__property double get_PeakBackgroundRatio()
			{
				return mdbl_PeakBackgroundRatio ; 
			}
			__property void set_PeakBackgroundRatio(double value)
			{
				mdbl_PeakBackgroundRatio = value ; 
			}

			__property double get_SignalToNoiseThreshold()
			{
				return mdbl_SNThreshold ; 
			}
			__property void set_SignalToNoiseThreshold(double value)
			{
				mdbl_SNThreshold = value ; 
			}

			__property PEAK_FIT_TYPE get_PeakFitType()
			{
				return menm_FitType ; 
			}
			__property void set_PeakFitType(PEAK_FIT_TYPE value)
			{
				menm_FitType = value ; 
			}

		};
	}
}
