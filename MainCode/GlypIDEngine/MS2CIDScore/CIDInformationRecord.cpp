// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#include "CIDInformationRecord.h"

namespace Engine
{
	namespace MS2CIDScoring
	{
		CIDInformationRecord::CIDInformationRecord()
		{
			mint_msn_scan_num = 0 ; 
			mint_parent_scan_num = 0 ; 

			mint_msn_scan_level = 0 ; 
			mint_parent_scan_level = 0 ; 

			mdbl_parent_mz = 0 ;
			mdbl_mono_mz = 0 ; 
			mdbl_y1_mz = 0 ; 


			mint_parent_intensity = 0 ; 
			mint_mono_intensity = 0 ; 

			mshort_cs = -1 ;
			mdbl_fit = 1 ;			
			mdbl_mono_mw = 0 ;
			mdbl_average_mw = 0 ; 
			

			mdbl_cid_score = 0 ; 
			mdbl_cid_score_p_value = 1 ; 
			mbln_contains_oxonium_ion = false ; 

			mstr_pep_seq_name = "" ;
			mstr_pro_seq_name = "" ;
			mstr_glyco_site = "" ; 
			mdbl_seq_mass = 0 ; 
			mdbl_glycan_mass = 0 ; 
			mstr_glycan_composition = " " ; 
			mdbl_mass_error = 0; 
			mbln_mass_error_in_ppm = true ; 
			mdbl_parent_scan_time = 0.0 ; 
			mdbl_cid_scan_time = 0.0 ; 

			mstr_glyco_site = ""; 

			
		}

		CIDInformationRecord::~CIDInformationRecord()
		{

		}

		void CIDInformationRecord::Clear()
		{
			mint_msn_scan_num = 0 ; 
			mint_parent_scan_num = 0 ; 

			mint_msn_scan_level = 0 ; 
			mint_parent_scan_level = 0 ; 

			mdbl_parent_mz = 0 ;
			mdbl_mono_mz = 0 ; 
			mdbl_y1_mz = 0 ; 


			mint_parent_intensity = 0 ; 
			mint_mono_intensity = 0 ; 

			mshort_cs = -1 ;
			mdbl_fit = 1 ;			
			mdbl_mono_mw = 0 ;
			mdbl_average_mw = 0 ; 
			mdbl_cid_score = 0 ; 
			mdbl_cid_score_p_value = 1 ; 
			mbln_contains_oxonium_ion = false ; 

			mstr_pep_seq_name = "" ;
			mstr_pro_seq_name = "" ;
			mstr_glyco_site = "" ; 
			mstr_nglyco_site = "";
			mdbl_seq_mass = 0 ; 
			mdbl_glycan_mass = 0 ; 
			mstr_glycan_composition = " " ; 
			mdbl_mass_error = 0; 
			mbln_mass_error_in_ppm = true ; 
			mdbl_parent_scan_time = 0.0 ; 
			mdbl_cid_scan_time = 0.0 ; 


		}

		void CIDInformationRecord::AddInfoToCIDRecord(int msnScan, int parentScan, int msnLevel, int parentLevel, double parentMz, 
			double monoMz, int parentIntensity, int monoIntensity, short charge, double fit, double monoMass, double avgMass, double cidScore, 
			double pvalue, bool contains_oxonium_ion, double y1_ion, double parent_time, double time)
		{
			mint_msn_scan_num  = msnScan ; 
			mint_parent_scan_num = parentScan ; 
			mint_msn_scan_level = msnLevel ; 
			mint_parent_scan_level = parentLevel ; 
			mdbl_parent_mz = parentMz ; 
			mdbl_mono_mz = monoMz ; 
			mint_parent_intensity = parentIntensity ; 
			mint_mono_intensity = monoIntensity ; 
			mshort_cs = charge ; 
			mdbl_fit = fit ; 
			mdbl_mono_mw = monoMass ; 
			mdbl_average_mw = avgMass ; 
			mdbl_cid_score = cidScore ;
			mdbl_cid_score_p_value = pvalue ; 
			mbln_contains_oxonium_ion = contains_oxonium_ion ; 
			mdbl_y1_mz = y1_ion ; 
			mdbl_parent_scan_time = parent_time ; 
			mdbl_cid_scan_time = time ; 
		}

		void CIDInformationRecord::AddSearchInfoToCIDRecord(std::string proSeq, double pepMass, double glycanMass, std::string glycanComp, std::string pepSeq, std::string sitepositions, bool constains_site, double massError, bool massErrorInPPM) 
		{
			mstr_pep_seq_name = pepSeq ; 
			mstr_pro_seq_name = proSeq ;
			mstr_nglyco_site = sitepositions; 
			mdbl_seq_mass = pepMass ; 
			mdbl_glycan_mass = glycanMass ; 
			mstr_glycan_composition = glycanComp ; 
			mdbl_mass_error = massError ; 
			mbln_mass_error_in_ppm = massErrorInPPM ; 
			if(constains_site)
				mstr_glyco_site = "yes" ; 
			else
				mstr_glyco_site = "no" ; 

		}

		void CIDInformationRecord::AddPeaks(Engine::PeakProcessing::PeakData &peak_data)
		{
			m_peak_data.Clear() ; 
			Engine::PeakProcessing::Peak peak ; 
			for (int i = 0 ; i < peak_data.GetNumPeaks() ; i++)
			{
				peak_data.GetPeak(i, peak) ; 
				m_peak_data.AddPeak(peak) ; 
			}
		}

		
	
	}

}