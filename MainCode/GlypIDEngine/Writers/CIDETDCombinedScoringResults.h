// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
#pragma once
#include <vector>
#include <map> 
#include "../MS2CIDETDCombinedScore/CIDETDCombinedInformationRecord.h"

#include <fstream> 


namespace Engine
{
	namespace Writers
	{
		//! This class stores information from HCD/CID scoring along with peptide search results for an LC-MS/MS dataset.
		class CIDETDCombinedScoringResults
		{
			
			//!keeping cid records as vectors for the moment but might shift to deque structures
			std::vector<Engine::MS2CIDETDCombinedScoring::CIDETDCombinedInformationRecord> mvect_score_records	 ; 			


			int mint_num_records ; 			
		public:					
			CIDETDCombinedScoringResults() ; 
			~CIDETDCombinedScoringResults() ; 
			
			void SaveResultsCombinedCIDETDScoring(char *etdFile) ; 			
			void LoadResultsCombinedCIDETDScoring(char *inpFile) ; 
			void ClearCombinedScoringCIDETDRecords() ; 
			
			void AddMultipleCombinedCIDETDRecords(std::vector<Engine::MS2CIDETDCombinedScoring::CIDETDCombinedInformationRecord> &vect_combined_records) ; 			
			void AddCombinedScoringCIDETDRecord(Engine::MS2CIDETDCombinedScoring::CIDETDCombinedInformationRecord score_record) ; 		
		};
	}
}