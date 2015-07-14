// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#include "CIDHCDCombinedInformationRecord.h"

namespace Engine
{
	namespace MS2CIDHCDCombinedScoring
	{
		CIDHCDCombinedInformationRecord::CIDHCDCombinedInformationRecord()
		{
			mint_parent_scan_num = 0 ; 
			mint_cid_scan_num = 0 ; 
			mint_hcd_scan_num = 0 ; 

			mdbl_parent_mz = 0 ;
			mdbl_mono_mz = 0 ; 
			mdbl_y1_mz = 0 ; 

			mint_parent_intensity = 0 ; 
			mint_mono_intensity = 0 ; 

			mshort_cs = -1 ;
			mdbl_fit = 1 ;			
			mdbl_mono_mw = 0 ;

			mdbl_cid_score = 0 ; 
			mdbl_cid_score_p_value = 1; 
			mdbl_hcd_score = 1 ; 
			mint_num_hcd_peaks_identified = 0 ; 

			menm_glycan_type = Engine::MS2HCDScoring::GLYCAN_TYPE::NA ; 
			mstr_pep_seq_name = "" ;
			mstr_pro_seq_name = "" ;
			mdbl_seq_mass = 0 ; 
			mdbl_glycan_mass = 0 ; 
			mstr_glycan_composition = " " ; 
			mstr_nglyco_site ="" ; 
			mdbl_ppm_error = 0; 

			mdbl_parent_scan_time = 0.0 ; 
			mdbl_cid_scan_time = 0.0 ; 
			mdbl_hcd_scan_time = 0.0 ; 
		}

		CIDHCDCombinedInformationRecord::~CIDHCDCombinedInformationRecord()
		{
		}

		void CIDHCDCombinedInformationRecord::AddInfoToCombinedRecord(int parentScan, int cidScan, int hcdScan, double parentMz, 
			double monoMz, int parentIntensity, int monoIntensity, short charge, double fit, double monoMass, double cidScore, double pvalue, double hcdScore,  
			bool contains_oxonium, std::vector <int> &peakIndices, Engine::MS2HCDScoring::GLYCAN_TYPE glycan_type, double y1MZ, double parent_time, double cid_time, double hcd_time)
		{
			mint_parent_scan_num = parentScan ; 
			mint_cid_scan_num = cidScan ; 
			mint_hcd_scan_num = hcdScan ; 
			
			mdbl_parent_mz = parentMz ; 
			mdbl_mono_mz = monoMz ; 
			mdbl_y1_mz = y1MZ ; 
			mdbl_parent_scan_time = parent_time ; 
			mint_parent_intensity = parentIntensity ; 
			mint_mono_intensity = monoIntensity ; 
			mshort_cs = charge ; 
			mdbl_fit = fit ; 
			mdbl_mono_mw = monoMass ; 
			mdbl_cid_score = cidScore ; 
			mdbl_cid_score_p_value = pvalue ;
			mdbl_cid_scan_time = cid_time ; 
			mdbl_hcd_score = hcdScore ; 
			mbln_contains_oxonium_ions = contains_oxonium  ; 
			mdbl_hcd_scan_time = hcd_time ; 

			mint_num_hcd_peaks_identified = peakIndices.size() ; 
			menm_glycan_type = glycan_type ; 
			
			for (int i=0 ; i < mint_num_hcd_peaks_identified ; i++)
			{
				mvect_peak_indices.push_back(peakIndices[i]) ; 
			}
		}
		void CIDHCDCombinedInformationRecord::AddSearchInfoToCombinedRecord(std::string proSeq, double pepMass, double glycanMass, std::string &glycanComp, std::string pepSeq, std::string site, double ppmError, bool containsSite) 
		{
			mstr_pep_seq_name = pepSeq ; 
			mstr_pro_seq_name = proSeq ;
			mstr_nglyco_site = site;
			mdbl_seq_mass = pepMass ; 
			mdbl_glycan_mass = glycanMass ; 
			mstr_glycan_composition = glycanComp ; 
			mdbl_ppm_error = ppmError ; 
			if(containsSite)
				mstr_glyco_site = "yes";
			else
				mstr_glyco_site = "no" ; 			

		}

	}

}