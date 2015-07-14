// Written by Anoop Mayampurath

// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

#pragma once
#include "MS2HCDScore/HCDInformationRecord.h" 
namespace GlypID
{
	namespace HCDScoring
	{
		public __gc class clsHCDScoringScanResults
		{
		public:
			///! scan number of MSn_scan
			int mint_msn_scan_num ; 
			//! msNLevel of MSn_scan
			int mint_msn_scan_level ; 
			//! parent scam
			int mint_parent_scan_num ; 
			//! msNLevel of parent scan
			int mint_parent_scan_level ; 
			//! m/z value of parent
			double mdbl_parent_mz ;
			//! mono m/z value of parent
			double mdbl_mono_mz ; 
			//! charge state 
			short mshort_cs ;			
			//! fit value .
			double mdbl_fit ;
			//! monoisotopic mw of feature.
			double mdbl_mono_mw ;			
			//! intensity of monoisotopic peak observed.
			int mint_mono_intensity  ;
			//! intensity of parent peak observed.
			int mint_parent_intensity  ;
			//! hcd  score
			double mdbl_hcd_score ; 
			//! identified peaks
			int marr_peak_indices __gc [] ; 
			//! number of identified peals
			int mint_num_peaks_observed ; 
			//! glycan type
			Engine::MS2HCDScoring::GLYCAN_TYPE menm_glycan_type ; 



		public:
			clsHCDScoringScanResults(void);
			~clsHCDScoringScanResults(void);
			void Set(Engine::MS2HCDScoring::HCDInformationRecord &hcdRecord)
			{

				mint_msn_scan_num = hcdRecord.mint_msn_scan_num ; 
				mint_msn_scan_level = hcdRecord.mint_msn_scan_level ; 
				mint_parent_scan_num = hcdRecord.mint_parent_scan_num ; 
				mint_parent_scan_level = hcdRecord.mint_parent_scan_level ; 
				mdbl_parent_mz = hcdRecord.mdbl_parent_mz ; 
				mdbl_mono_mz = hcdRecord.mdbl_mono_mz ; 
				mshort_cs = hcdRecord.mshort_cs ; 
				mdbl_fit = hcdRecord.mdbl_fit ; 
				mdbl_mono_mw = hcdRecord.mdbl_mono_mw ; 
				mint_mono_intensity = hcdRecord.mint_mono_intensity ; 
				mint_parent_intensity = hcdRecord.mint_parent_intensity ; 
				mdbl_hcd_score = hcdRecord.mdbl_hcd_score ; 
				mint_num_peaks_observed = hcdRecord.mint_num_peaks_identified ; 
				menm_glycan_type = hcdRecord.menm_glycan_type; 
				marr_peak_indices = new int __gc [hcdRecord.mint_num_peaks_identified] ; 
				for (int i = 0 ; i < hcdRecord.mint_num_peaks_identified ; i++)
				{
					marr_peak_indices[i] = hcdRecord.mvect_peak_indices[i] ; 
				}
			}
		};
	}
}
