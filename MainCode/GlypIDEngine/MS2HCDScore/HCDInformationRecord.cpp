// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#include "HCDInformationRecord.h"

namespace Engine
{
	namespace MS2HCDScoring
	{
		HCDInformationRecord::HCDInformationRecord()
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
			mdbl_hcd_score = 1 ; 

			menm_glycan_type = Engine::MS2HCDScoring::GLYCAN_TYPE::NA ; 

			mstr_pep_seq_name = "" ;
			mstr_pro_seq_name = "" ;
			mstr_glyco_site = "" ; 
			mdbl_seq_mass = 0 ; 
			mdbl_glycan_mass = 0 ; 
			mstr_glycan_composition = " " ; 
			mdbl_mass_error = 0; 
			mbln_mass_error_in_ppm = true ; 
			mdbl_parent_scan_time = 0.0 ; 
		}

		HCDInformationRecord::~HCDInformationRecord()
		{
		}

		void HCDInformationRecord::Clear()
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
			mdbl_hcd_score = 1 ; 

			menm_glycan_type = Engine::MS2HCDScoring::GLYCAN_TYPE::NA ; 

			mstr_pep_seq_name = "" ;
			mstr_pro_seq_name = "" ;
			mstr_glyco_site = "" ; 
			mdbl_seq_mass = 0 ; 
			mdbl_glycan_mass = 0 ; 
			mstr_glycan_composition = " " ; 
			mdbl_mass_error = 0; 
			mbln_mass_error_in_ppm = true ; 
			mdbl_parent_scan_time = 0.0 ; 
		}

		void HCDInformationRecord::AddInfoToHCDRecord(int msnScan, int parentScan, int msnLevel, int parentLevel, double parentMz, 
			double monoMz, int parentIntensity, int monoIntensity, short charge, double fit, double monoMass, double hcdScore, std::vector <int> &peakIndices,
			Engine::MS2HCDScoring::GLYCAN_TYPE glycan_type)
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
			mdbl_hcd_score = hcdScore ; 	
			mint_num_peaks_identified = peakIndices.size() ; 
			menm_glycan_type = glycan_type ; 
			
			for (int i=0 ; i < mint_num_peaks_identified ; i++)
			{
				mvect_peak_indices.push_back(peakIndices[i]) ; 
			}

		}

		void HCDInformationRecord::AddSearchInfoToHCDRecord(std::string proSeq, double pepMass, double glycanMass, std::string &glycanComp, std::string pepSeq, bool constains_site, double massError, bool massErrorInPPM) 
		{
			mstr_pep_seq_name = pepSeq ; 
			mstr_pro_seq_name = proSeq ;
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

		
	
	}

}