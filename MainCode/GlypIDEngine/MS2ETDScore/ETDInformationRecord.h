// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#pragma once
#include <vector>

namespace Engine
{
	namespace MS2ETDScoring
	{
		//! class to store information that logs into the log file of DeconMSn
		class  ETDInformationRecord
		{
		public:			
			//! scan number of MSn_scan
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
			//! average mw
			double mdbl_average_mw ; 
			//! intensity of monoisotopic peak observed.
			int mint_mono_intensity  ;
			//! intensity of parent peak observed.
			int mint_parent_intensity  ;
			//! cid  score
			double mdbl_etd_score ; 
			//! etd score p value
			double mdbl_etd_fdr_score;			
			//! protein name
			std::string mstr_pro_seq_name ;
			//! peptide match
			std::string mstr_pep_seq_name;
			//! peptide sequnce mass
			double mdbl_seq_mass ; 
			//! peptide has glycosylation sote
			std::string mstr_glyco_site ; 
			//! glycan mass
			double mdbl_glycan_mass ; 
			//! glycan composition
			std::string mstr_glycan_composition ; 
			//! mass error (either ppm or Da)
			double mdbl_mass_error ;	
			//! to set if the mass error was in ppm or DA
			bool mbln_mass_error_in_ppm ; 			
			//! parent elution time, not tech the peak of the elution profile, but still
			double mdbl_parent_scan_time ; 			

			//! default constructor
			ETDInformationRecord() ; 
			//! destructor.
			~ETDInformationRecord() ;
			//! record details
			void AddInfoToETDRecord(int msnScan, int parentScan, int msnLevel, int parentLevel, double parentMz, 
			double monoMz, int parentIntensity, int monoIntensity, short charge, double fit, double monoMass, double avgMass,double parent_time) ; 
			//! clear up record
			void Clear() ; 
			//! adding peptide search information
			void AddSearchInfoToETDRecord(std::string proSeq, double pepMass, double glycanMass, std::string &glycanComp, std::string pepSeq,   double massError, bool massErrorInPPM, double etdScore, double fdr) ; 

		

		};
	}
}