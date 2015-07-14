// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
#pragma once
#include <vector>
#include <map> 
#include "../MS2HCDScore/HCDInformationRecord.h"

#include <fstream> 


namespace Engine
{
	namespace Writers
	{
		//! This class stores information from HCD scoring an LC-MS/MS dataset.
		class HCDScoringResults
		{
			
			//!keeping hcd records as vectors for the moment but might shift to deque structures
			std::vector<Engine::MS2HCDScoring::HCDInformationRecord> mvect_hcd_records ; 			


			int mint_num_records ; 
			
			

			
		public:			

			
			HCDScoringResults() ; 
			~HCDScoringResults() ; 
			
			void SaveResultsHCD(char *hcdFile) ; 			
			void LoadResultsHCD(char *inpFile) ; 
			void ClearHCDRecords() ; 
			
			void AddMultipleHCDRecords(std::vector<Engine::MS2HCDScoring::HCDInformationRecord> &vect_hcd_records) ; 			
			void AddHCDRecord(Engine::MS2HCDScoring::HCDInformationRecord hcd_record) ; 
			
		
		};
	}
}