// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University

#include "HCDInformationRecord.h"
#include "../PeakProcessor/PeakProcessor.h"
#include "../PeakProcessor/PeakData.h"
#include "../MS2CIDScore/CIDInformationRecord.h"
#include "../Utilities/util.h"

#include <vector>


namespace Engine
{
	namespace MS2HCDScoring
	{		

		class HCDScoring
		{
			private:
				// minimum number of peaks 
				int mint_min_number_peaks ; 
				// min mz 
				double mdbl_min_mz ; 
				// max mz 
				double mdbl_max_mz ; 
				// array of mzs
				double marr_theoretical_peaks[7] ; 
				// charge carrier mass
				double mdbl_cc_mass; 
				// theoretical intensitiy distribution for each class
				float marr_theoretical_hm[7] ; 
				float marr_theoretical_ca[7] ; 
				float marr_theoretical_cs[7] ; 
				float marr_theoretical_hybrid[7] ; 
				// identified peaks
				std::vector <int> mint_peak_indices ; 				


			public:
				
				HCDScoring(void); 
				~HCDScoring(void) ; 

				//double CalculateHCDScore(Engine::PeakProcessing::PeakData &pk_data, Engine::PeakProcessing::PeakData &detected_pk_data) ; 


				double CalculateHCDScore(Engine::PeakProcessing::PeakData &pk_data, std::vector <int> &identified_peak_indices) ; 
				double CalculateHCDScorePValue(Engine::PeakProcessing::PeakData &pk_data, std::vector <int> &identified_peak_indices) ; 
				double DetermineGlycanType(Engine::PeakProcessing::PeakData &pk_data, std::vector <int> &identified_peak_indices, Engine::MS2HCDScoring::GLYCAN_TYPE &glycan_type) ; 								
				void SetOptions(int mint_num_peaks, double min_mz, double max_mz) ; 
				double DetermineY1IonThroughCorrelation(Engine::PeakProcessing::PeakData &hcd_peaks, Engine::PeakProcessing::PeakData &cid_peaks, double min_mz, double max_mz) ;				
				double DetermineY1IonThroughPeptideSearching(Engine::PeakProcessing::PeakData &hcd_data,  Engine::MS2CIDScoring::CIDInformationRecord &cid_record) ; 
				//double DetermineY1IonThroughPeptideSearching(Engine::PeakProcessing::PeakData &hcd_data, std::vector<Engine::MS2CIDScoring::CIDInformationRecord> &vect_cid_records) ; 
				void GetOptions(int &mint_num_peaks, double &min_mz, double &max_mz) ; 

		}; 
	}
}




