// Written by Anoop Mayampurath

// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

#pragma once
#include "MS2ETDScore/ETDInformationRecord.h" 
#include "MS2HCDScore/HCDInformationRecord.h" 
#include "clsCIDScoringScanResults.h" 
#include "clsHCDScoringScanResults.h" 
#include "clsHornTransformResults.h" 
namespace GlypID
{
	namespace ETDScoring
	{
		public __gc class clsETDScoringScanResults
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
			//! average mw of feature
			double mdbl_average_mw ; 
			//! intensity of monoisotopic peak observed.
			int mint_mono_intensity  ;
			//! intensity of parent peak observed.
			int mint_parent_intensity  ;
			//! etd  score
			double mdbl_etd_score ; 
			//! etd fdr
			double mdbl_etd_fdr_score ; 
			//! protein name
			System::String *mstr_pro_seq_name ;
			//! peptide match
			System::String  *mstr_pep_seq;
			//! peptide sequnce mass
			double mdbl_seq_mass ; 		
			//! peptide has glycosylation sote
			System::String  *mstr_glyco_site ; 
			//! position of glycosylation site
			System::String *mstr_nglyco_site ; 
			//! position of oglycosylation site
			System::String *mstr_oglyco_site;
			//! glycan mass
			double mdbl_glycan_mass ; 
			//! glycan composition
			System::String  *mstr_glycan_composition ; 
			//! mass error (either ppm or Da)
			double mdbl_mass_error ;	
			//! to set if the mass error was in ppm or DA
			bool mbln_mass_error_in_ppm ; 			
			//! parent elution time, not tech the peak of the elution profile, but still
			double mdbl_parent_scan_time ; 		
			//! msn_elution time
			double mdbl_etd_scan_time ; 
			//! true hit or false hit
			bool mbln_true_hit ; 
			// ! y1 mz 
			double mdbl_y1_mz ;
			// hcd score
			double mdbl_hcd_score ; 
			//  glycan type				
			Engine::MS2HCDScoring::GLYCAN_TYPE menm_glycan_type ; 
			// cid score
			double mdbl_cid_score ; 
			// cid sequenceing score
			double mflt_cid_sequencing_score ; 
			// glycan sequencng
			System::String *mstr_glycan_sequence; 
			


		public:
			clsETDScoringScanResults(void);
			~clsETDScoringScanResults(void);


			

			
	

			void Set(Engine::MS2ETDScoring::ETDInformationRecord &etdRecord)
			{

				mint_msn_scan_num = etdRecord.mint_msn_scan_num ; 
				mint_msn_scan_level = etdRecord.mint_msn_scan_level ; 
				mint_parent_scan_num = etdRecord.mint_parent_scan_num ; 
				mint_parent_scan_level = etdRecord.mint_parent_scan_level ; 
				mdbl_parent_mz = etdRecord.mdbl_parent_mz ; 
				mdbl_mono_mz = etdRecord.mdbl_mono_mz ; 
				mshort_cs = etdRecord.mshort_cs ; 
				mdbl_fit = etdRecord.mdbl_fit ; 
				mdbl_mono_mw = etdRecord.mdbl_mono_mw ; 
				mint_mono_intensity = etdRecord.mint_mono_intensity ; 
				mint_parent_intensity = etdRecord.mint_parent_intensity ; 
				mdbl_etd_score = etdRecord.mdbl_etd_score ; 
				mdbl_average_mw = etdRecord.mdbl_average_mw ; 
				
			}

			void Set(GlypID::CIDScoring::clsCIDScoringScanResults &cidResult)
			{

				mint_msn_scan_num = cidResult.mint_msn_scan_num ; 
				mint_msn_scan_level = cidResult.mint_msn_scan_level ; 
				mint_parent_scan_num = cidResult.mint_parent_scan_num ; 
				mint_parent_scan_level = cidResult.mint_parent_scan_level ; 
				mdbl_parent_mz = cidResult.mdbl_parent_mz ; 
				mdbl_mono_mz = cidResult.mdbl_mono_mz ; 
				mshort_cs = cidResult.mshort_cs ; 
				mdbl_fit = cidResult.mdbl_fit ; 
				mdbl_mono_mw = cidResult.mdbl_mono_mw ; 
				mint_mono_intensity = cidResult.mint_mono_intensity ; 
				mint_parent_intensity = cidResult.mint_parent_intensity ; 
			}

			void Set(Engine::MS2CIDScoring::CIDInformationRecord &cidRecord)
			{				
				mint_msn_scan_num = cidRecord.mint_msn_scan_num ; 
				mint_msn_scan_level = cidRecord.mint_msn_scan_level ; 
				mint_parent_scan_num = cidRecord.mint_parent_scan_num ; 
				mint_parent_scan_level = cidRecord.mint_parent_scan_level ; 
				mdbl_parent_mz = cidRecord.mdbl_parent_mz ; 
				mdbl_mono_mz = cidRecord.mdbl_mono_mz ; 
				mshort_cs = cidRecord.mshort_cs ; 
				mdbl_fit = cidRecord.mdbl_fit ; 
				mdbl_mono_mw = cidRecord.mdbl_mono_mw ; 
				mdbl_average_mw = cidRecord.mdbl_average_mw ; 
				mint_mono_intensity = cidRecord.mint_mono_intensity ; 
				mint_parent_intensity = cidRecord.mint_parent_intensity ; 
				mstr_pro_seq_name = cidRecord.mstr_pro_seq_name.c_str(); 
				mstr_pep_seq = cidRecord.mstr_pep_seq_name.c_str() ; 
				mstr_glycan_composition= cidRecord.mstr_glycan_composition.c_str() ; 
				mdbl_glycan_mass = cidRecord.mdbl_glycan_mass ; 
				mbln_mass_error_in_ppm = cidRecord.mbln_mass_error_in_ppm ;
				mdbl_mass_error = cidRecord.mdbl_mass_error;
				mdbl_seq_mass = cidRecord.mdbl_seq_mass ;	
				mdbl_y1_mz = cidRecord.mdbl_y1_mz ; 
				mstr_glyco_site = cidRecord.mstr_glyco_site.c_str() ; //yes or no value
				mstr_nglyco_site = cidRecord.mstr_nglyco_site.c_str() ; // Actual site in protein
			}
		};
	}
}
