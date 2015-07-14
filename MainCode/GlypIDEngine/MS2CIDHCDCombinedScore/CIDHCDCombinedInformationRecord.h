// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#pragma once
#include "../MS2HCDScore/HCDInformationRecord.h"
namespace Engine
{
	namespace MS2CIDHCDCombinedScoring
	{
		//! class to store information 
		class  CIDHCDCombinedInformationRecord
		{
		public:			
			//! scan number of Parent Scan
			int mint_parent_scan_num ; 
			//! scan number of CID scan
			int mint_cid_scan_num ; 
			//! scan number of HCD scan
			int mint_hcd_scan_num ; 
			//! parent m/z
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
			//! cid  score
			double mdbl_cid_score ; 
			//!cid p_value score
			double mdbl_cid_score_p_value ;
			//! hcd score
			double mdbl_hcd_score ;
			//! cid contains oxonium ions
			bool mbln_contains_oxonium_ions ; 
			//! peak indices
			std::vector <int> mvect_peak_indices ; 
			//! number of identified hcd peak
			int mint_num_hcd_peaks_identified ; 
			//! glycan type
			Engine::MS2HCDScoring::GLYCAN_TYPE menm_glycan_type ; 
			//! protein name
			std::string mstr_pro_seq_name ;
			//! peptide match
			std::string mstr_pep_seq_name;
			//! peptide sequence mass
			double mdbl_seq_mass ; 
			//! contains glycosylation site
			std::string mstr_glyco_site ; 
			//! glycan mass
			double mdbl_glycan_mass ; 
			//! glycan composition
			std::string mstr_glycan_composition ; 
			//! Site of glycosylation
			std::string mstr_nglyco_site; 
			//! ppm error
			double mdbl_ppm_error ;
			//! y1 ion
			double mdbl_y1_mz ; 
			//! parent elution time, not tech the peak of the elution profile, but still
			double mdbl_parent_scan_time ; 
			//! cid elution time
			double mdbl_cid_scan_time; 
			//! hcd elution time
			double mdbl_hcd_scan_time ; 
			//! default constructor
			CIDHCDCombinedInformationRecord() ; 
			//! destructor.
			~CIDHCDCombinedInformationRecord() ;
			//! record details
			void AddInfoToCombinedRecord(int parentScan, int cidScan, int hcdScan, double parentMz, double monoMz, int parentIntensity, int monoIntensity, short charge, double fit, double monoMass, double cidScore, double pvalue, double hcdScore, 
				bool contains_oxonium, std::vector <int> &peakIndices, Engine::MS2HCDScoring::GLYCAN_TYPE glycan_type, double y1MZ, double parent_time, double cid_time, double hcd_time) ; 
			void AddSearchInfoToCombinedRecord(std::string proSeq, double pepMass, double glycanMass, std::string &glycanComp, std::string pepSeq, std::string site, double ppmError, bool containsSite) ; 
		

		};
	}
}