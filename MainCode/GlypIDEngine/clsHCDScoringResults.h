// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#pragma once
#include "Writers/HCDScoringResults.h"
#include "GlypIDEngineUtils.h" 
#include "clsRawData.h"

// NEEDS work

namespace GlypID
{
	namespace Results
	{
		public __gc class clsHCDScoringResults
		{
			int mint_percent_done ; 
			GlypID::Readers::FileType menmFileType ; 
		public:
			Engine::Writers::HCDScoringResults __nogc *mobj_hcd_scoring_results ;
			void SetHCDScoringResults(Engine::Writers::HCDScoringResults *results) ; 			
			void WriteResults(System::String *fileName) ; 
			void ReadResults(System::String *fileName) ; 
			clsHCDScoringResults(void);
			~clsHCDScoringResults(void);		

			__property GlypID::Readers::FileType get_FileType()
			{
				return menmFileType ; 
			}

			__property void set_FileType(GlypID::Readers::FileType fileType)
			{
				menmFileType = fileType ; 
			}
		};
			
	}
}