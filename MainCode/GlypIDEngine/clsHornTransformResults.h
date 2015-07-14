// Written by Navdeep Jaitly for the Department of Energy (PNNL, Richland, WA)
// Copyright 2006, Battelle Memorial Institute
// E-mail: navdeep.jaitly@pnl.gov
// Website: http://ncrr.pnl.gov/software
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

#pragma once
#include "HornTransform/IsotopeFitRecord.h" 
namespace GlypID
{
	namespace HornTransform
	{
		public __gc class clsHornTransformResults
		{
		public:
			//! peak index of the peak.
			int mint_peak_index ; 
			//! scan number of peak
			int mint_scan_num ; 
			//! charge state 
			short mshort_cs ;
			//! intensity of feature.
			int mint_abundance ;
			//! m/z value of most abundant peak in the feature.
			double mdbl_mz ;
			//! fit value .
			double mdbl_fit ;
			//! average mw for the feature. 
			double mdbl_average_mw ;
			//! monoisotopic mw of feature.
			double mdbl_mono_mw ;
			//! mw at the most intense isotope.
			double mdbl_most_intense_mw ; 
			//! full width at half maximum of the peak.
			double mdbl_fwhm ; 
			//! signal to noise for the most intense isotopic peak.
			double mdbl_sn ; 
			//! intensity of monoisotopic peak observed.
			int mint_mono_intensity  ;
			//! intensity of the third isotopic peak observed. Used by other software for processing of O16/O18  data.
			int mint_iplus2_intensity  ;
			//! difference between obsered m/z and m/z from theoretical distribution of composition from Averagine
			double mdbl_delta_mz ; 
			//! need multiple isotopes to determine charge
			bool mbln_need_multiple_isotopes ; 
			//! number of isotope peaks
			int mint_num_isotopes_observed ; 
			//! array of indices of peak tops
			int marr_isotope_peak_indices __gc [] ; 

		public:
			clsHornTransformResults(void);
			~clsHornTransformResults(void);
			void Set(Engine::HornTransform::IsotopeFitRecord &fitRecord)
			{
				mint_peak_index = fitRecord.mint_peak_index ; 
				mint_scan_num = fitRecord.mint_scan_num ; 
				mshort_cs = fitRecord.mshort_cs ;
				mint_abundance = fitRecord.mint_abundance ;
				mdbl_mz = fitRecord.mdbl_mz ;
				mdbl_fit = fitRecord.mdbl_fit ;
				mdbl_average_mw = fitRecord.mdbl_average_mw ;
				mdbl_mono_mw = fitRecord.mdbl_mono_mw ;
				mdbl_most_intense_mw = fitRecord.mdbl_most_intense_mw ; 
				mdbl_fwhm = fitRecord.mdbl_fwhm ; 
				mdbl_sn = fitRecord.mdbl_sn ; 
				mint_mono_intensity = fitRecord.mint_mono_intensity ;
				mint_iplus2_intensity = fitRecord.mint_iplus2_intensity ;
				mdbl_delta_mz = fitRecord.mdbl_delta_mz ; 
				mint_num_isotopes_observed = fitRecord.mint_num_isotopes_observed ; 
				marr_isotope_peak_indices = new int __gc [mint_num_isotopes_observed] ; 
				for (int i = 0 ; i < mint_num_isotopes_observed ; i++)
				{
					marr_isotope_peak_indices[i] = fitRecord.marr_isotope_peak_indices[i] ; 
				}
			}
		};
	}
}
