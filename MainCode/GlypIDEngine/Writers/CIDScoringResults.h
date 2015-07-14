// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
#pragma once
#include <vector>
#include <map> 
#include "../MS2CIDScore/CIDInformationRecord.h"

#include <fstream> 


namespace Engine
{
	namespace Writers
	{
		//! This class stores information from HCD scoring an LC-MS/MS dataset.
		class CIDScoringResults
		{
			
			//!keeping cid records as vectors for the moment but might shift to deque structures
			std::vector<Engine::MS2CIDScoring::CIDInformationRecord> mvect_cid_records	 ; 			


			int mint_num_records ; 			
		public:					
			CIDScoringResults() ; 
			~CIDScoringResults() ; 
			
			void SaveResultsCID(char *hcdFile) ; 			
			void LoadResultsCID(char *inpFile) ; 
			void ClearCIDRecords() ; 
			
			void AddMultipleCIDRecords(std::vector<Engine::MS2CIDScoring::CIDInformationRecord> &vect_cid_records) ; 			
			void AddCIDRecord(Engine::MS2CIDScoring::CIDInformationRecord cid_record) ; 		
		};
	}
}