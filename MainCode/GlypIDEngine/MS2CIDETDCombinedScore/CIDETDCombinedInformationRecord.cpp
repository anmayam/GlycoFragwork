// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#include "CIDETDCombinedInformationRecord.h"

namespace Engine
{
	namespace MS2CIDETDCombinedScoring
	{
		CIDETDCombinedInformationRecord::CIDETDCombinedInformationRecord()
		{
			mint_parent_scan_num = 0 ; 
			mint_cid_scan_num = 0 ; 
			mint_etd_scan_num = 0 ; 

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
			mdbl_etd_score = 0 ; 
			mdbl_etd_score_fdr = 100 ; 
			

			mstr_pep_seq_name = "" ;
			mstr_pro_seq_name = "" ;
			mdbl_seq_mass = 0 ; 
			mdbl_glycan_mass = 0 ; 
			mstr_glycan_composition = " " ; 
			mdbl_ppm_error = 0; 

			mdbl_parent_scan_time = 0.0 ; 
		}

		CIDETDCombinedInformationRecord::~CIDETDCombinedInformationRecord()
		{
		}
		void CIDETDCombinedInformationRecord::Clear()
		{
			mint_parent_scan_num = 0 ; 
			mint_cid_scan_num = 0 ; 
			mint_etd_scan_num = 0 ; 

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
			mdbl_etd_score = 0 ; 
			mdbl_etd_score_fdr = 100 ; 
			

			mstr_pep_seq_name = "" ;
			mstr_pro_seq_name = "" ;
			mdbl_seq_mass = 0 ; 
			mdbl_glycan_mass = 0 ; 
			mstr_glycan_composition = " " ; 
			mdbl_ppm_error = 0; 

			mdbl_parent_scan_time = 0.0 ; 

		}

		void CIDETDCombinedInformationRecord::AddInfoToCombinedRecord(int parentScan, int cidScan, int etdScan, double parentMz, 
			double monoMz, int parentIntensity, int monoIntensity, short charge, double fit, double monoMass, double cidScore, double pvalue, double etdScore,  
			double fdr, bool contains_oxonium,  double y1MZ, double parent_time)
		{
			mint_parent_scan_num = parentScan ; 
			mint_cid_scan_num = cidScan ; 
			mint_etd_scan_num = etdScan ; 
			
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
			mdbl_etd_score = etdScore ; 
			mdbl_etd_score_fdr = fdr ; 
			mbln_contains_oxonium_ions = contains_oxonium  ; 
			mdbl_y1_mz = y1MZ ; 
			mdbl_parent_scan_time = parent_time ; 
			

		}
		void CIDETDCombinedInformationRecord::AddSearchInfoToCombinedRecord(std::string proSeq, double pepMass, double glycanMass, std::string &glycanComp, std::string pepSeq, double ppmError, bool containsSite) 
		{
			mstr_pep_seq_name = pepSeq ; 
			mstr_pro_seq_name = proSeq ;
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