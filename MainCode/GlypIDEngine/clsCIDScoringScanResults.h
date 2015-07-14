// Written by Anoop Mayampurath

// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

#pragma once
#include "MS2CIDScore/CIDInformationRecord.h" 
namespace GlypID
{
	namespace CIDScoring
	{
		public __gc class clsCIDScoringScanResults
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
			//! cid  score
			double mdbl_cid_score ; 
			//! putative peptide
			System::String *mstr_pro_seq_name ;
			//! putative protein
			System::String  *mstr_pep_seq_name;


		public:
			clsCIDScoringScanResults(void);
			~clsCIDScoringScanResults(void);
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
				mdbl_cid_score = cidRecord.mdbl_cid_score ; 
			}			
		};
	}
}
