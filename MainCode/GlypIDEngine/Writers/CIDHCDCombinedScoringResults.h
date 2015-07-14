// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
#pragma once
#include <vector>
#include <map> 
#include "../MS2CIDHCDCombinedScore/CIDHCDCombinedInformationRecord.h"

#include <fstream> 


namespace Engine
{
	namespace Writers
	{
		//! This class stores information from HCD/CID scoring along with peptide search results for an LC-MS/MS dataset.
		class CIDHCDCombinedScoringResults
		{
			
			//!keeping cid records as vectors for the moment but might shift to deque structures
			std::vector<Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord> mvect_score_records	 ; 			


			int mint_num_records ; 			
		public:					
			CIDHCDCombinedScoringResults() ; 
			~CIDHCDCombinedScoringResults() ; 
			
			void SaveResultsCombinedScoring(char *hcdFile) ; 			
			void LoadResultsCombinedScoring(char *inpFile) ; 
			void ClearCombinedScoringRecords() ; 
			
			void AddMultipleCombinedRecords(std::vector<Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord> &vect_combined_records) ; 			
			void AddCombinedScoringRecord(Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord score_record) ; 		
		};
	}
}