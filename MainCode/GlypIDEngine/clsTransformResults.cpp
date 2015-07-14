// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#include "clstransformresults.h"
#using <mscorlib.dll>
#include <algorithm>


namespace GlypID
{
	namespace Results
	{
		clsTransformResults::clsTransformResults(void)
		{
			mobj_lcms_results = NULL ; 
			menmFileType = GlypID::Readers::FileType::UNDEFINED ; 
		}

		clsTransformResults::~clsTransformResults(void)
		{
			if (mobj_lcms_results != NULL)
			{
				delete mobj_lcms_results ; 
				mobj_lcms_results = NULL ; 
			}
		}
		
		void clsTransformResults::SetLCMSTransformResults(Engine::Writers::LCMSTransformResults *results)
		{
			if (mobj_lcms_results != NULL)
			{
				delete mobj_lcms_results ; 
				mobj_lcms_results = NULL ; 
			}
			mobj_lcms_results = results ; 
		}

		int clsTransformResults::GetMinScan()
		{
			if(mobj_lcms_results == NULL)
				return -1 ; 
			return mobj_lcms_results->GetMinScan() ;
		}

		int clsTransformResults::GetMaxScan()
		{
			if(mobj_lcms_results == NULL)
				return -1 ; 
			return mobj_lcms_results->GetMaxScan() ;
		}

		bool clsTransformResults::IsDeisotoped()
		{
			if (mobj_lcms_results == NULL)
				throw new System::Exception(S"No results stored.") ;
			return mobj_lcms_results->IsDeisotoped() ; 
		}

		
		// Used by the stl algorithm sort to sort std::vector of peaks in descending order of mdbl_intensity.
		bool RawPeaksIntensityComparison(Engine::Writers::LCMSPeak<Engine::Writers::PeakMinInfo> &pk1, Engine::Writers::LCMSPeak<Engine::Writers::PeakMinInfo> &pk2)
		{
			if (pk1.mobj_peak.mflt_intensity < pk2.mobj_peak.mflt_intensity)
				return true ; 
			if (pk1.mobj_peak.mflt_intensity > pk2.mobj_peak.mflt_intensity)
				return false ; 
			return pk1.mobj_peak.mflt_mz < pk2.mobj_peak.mflt_mz ; 
		}

		void clsTransformResults::GetRawDataSortedInIntensity(GlypID::Writers::clsLCMSPeak* (&lcms_peaks) __gc [])
		{
			if (mobj_lcms_results == NULL)
			{
				lcms_peaks = NULL ; 
				return ; 
			}

			std::vector<Engine::Writers::LCMSPeak<Engine::Writers::PeakMinInfo> > vectPeaks ; 
			int num_peaks = mobj_lcms_results->GetNumPeaks() ; 
			vectPeaks.reserve(num_peaks) ; 
			mobj_lcms_results->GetAllPeaks(vectPeaks) ; 
			std::sort(vectPeaks.begin(), vectPeaks.end(), &RawPeaksIntensityComparison) ; 

			lcms_peaks = new GlypID::Writers::clsLCMSPeak* __gc [num_peaks] ; 
			int min_scan = mobj_lcms_results->GetMinScan() ; 
			int max_scan = mobj_lcms_results->GetMaxScan() ; 
			Engine::Writers::LCMSPeak<Engine::Writers::PeakMinInfo> pk ; 

			for (int pk_num = 0 ; pk_num < num_peaks ; pk_num++)
			{
				mint_percent_done = (pk_num * 100)/num_peaks ; 
				pk = vectPeaks[pk_num] ; 
				lcms_peaks[pk_num] = new GlypID::Writers::clsLCMSPeak(pk.mint_scan_num, pk.mobj_peak.mflt_mz, pk.mobj_peak.mflt_intensity) ; 
			}
			mint_percent_done = 100 ; 
		}

		void clsTransformResults::GetRawData(GlypID::Writers::clsLCMSPeak* (&lcms_peaks) __gc [])
		{
			if (mobj_lcms_results == NULL)
			{
				lcms_peaks = NULL ; 
				return ; 
			}

			int num_peaks = mobj_lcms_results->GetNumPeaks() ; 
			lcms_peaks = new GlypID::Writers::clsLCMSPeak* __gc [num_peaks] ; 
			int min_scan = mobj_lcms_results->GetMinScan() ; 
			int max_scan = mobj_lcms_results->GetMaxScan() ; 
			Engine::Writers::LCMSPeak<Engine::Writers::PeakMinInfo> pk ; 

			for (int pk_num = 0 ; pk_num < num_peaks ; pk_num++)
			{
				mint_percent_done = (pk_num * 100)/num_peaks ; 
				pk = mobj_lcms_results->GetPeak(pk_num) ; 
				lcms_peaks[pk_num] = new GlypID::Writers::clsLCMSPeak(pk.mint_scan_num, pk.mobj_peak.mflt_mz, pk.mobj_peak.mflt_intensity) ; 
			}
			mint_percent_done = 100 ; 
		}
		
		
		void clsTransformResults::GetSIC(int min_scan, int max_scan, float mz, float mz_tolerance, float (&peak_intensities) __gc [])
		{
			std::vector<float> vect_intensities ; 
			mobj_lcms_results->GetSIC(min_scan, max_scan, mz-mz_tolerance, mz+mz_tolerance, vect_intensities) ; 
			int num_scans = max_scan - min_scan + 1 ; 
			peak_intensities = new float __gc [num_scans] ; 

			for (int scan_num = min_scan ; scan_num <= max_scan ; scan_num++)
			{
				peak_intensities[scan_num-min_scan] = vect_intensities[scan_num-min_scan] ; 
			}
		}

		void clsTransformResults::GetScanPeaks(int scan_num, float (&peak_mzs) __gc [], float (&peak_intensities) __gc [])
		{
			std::vector<float> vect_mzs ; 
			std::vector<float> vect_intensities ; 

			mobj_lcms_results->GetScanPeaks(scan_num, vect_mzs, vect_intensities) ; 

			int num_pts = vect_intensities.size() ; 
			peak_mzs = new float __gc [num_pts] ; 
			peak_intensities = new float __gc [num_pts] ; 

			for (int pt_num = 0 ; pt_num < num_pts ; pt_num++)
			{
				peak_mzs[pt_num] = vect_mzs[pt_num] ; 
				peak_intensities[pt_num] = vect_intensities[pt_num] ; 
			}
		}

		int clsTransformResults::GetNumPeaks()
		{
			return mobj_lcms_results->GetNumPeaks() ; 
		}

		void clsTransformResults::ReadResults(System::String *fileName)
		{
			if (mobj_lcms_results == NULL)
			{
				mobj_lcms_results = new Engine::Writers::LCMSTransformResults() ; 
			}

			char fileNameCh[512] ; 
			GlypID::Utils::GetStr(fileName, fileNameCh) ; 
			mobj_lcms_results->LoadResults(fileNameCh) ; 
		}

		void clsTransformResults::WriteResults(System::String *fileName, bool save_signal_range)
		{
			char fileNameCh[512] ; 
			try
			{
				GlypID::Utils::GetStr(fileName, fileNameCh) ; 
				mobj_lcms_results->SaveResults(fileNameCh, save_signal_range) ; 
			}
			catch (char *mesg)
			{
				System::String *exception_msg = new System::String(mesg) ; 
				throw new System::Exception(exception_msg) ; 
			}
		}

		void clsTransformResults::WriteScanResults(System::String *fileName)
		{
			char fileNameCh[512] ; 
			try
			{
				GlypID::Utils::GetStr(fileName, fileNameCh) ; 				
				mobj_lcms_results->SaveScanResults(fileNameCh, true) ; 
			}
			catch (char *mesg)
			{
				System::String *exception_msg = new System::String(mesg) ; 
				throw new System::Exception(exception_msg) ; 
			}
		}


	}
}