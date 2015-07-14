// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#pragma once
#include "Writers/CIDScoringResults.h"
#include "GlypIDEngineUtils.h" 
#include "clsRawData.h"

// NEEDS work

namespace GlypID
{
	namespace Results
	{
		public __gc class clsCIDScoringResults
		{
			int mint_percent_done ; 
			GlypID::Readers::FileType menmFileType ; 
		public:
			Engine::Writers::CIDScoringResults __nogc *mobj_cid_scoring_results ;
			clsCIDScoringResults(void);
			~clsCIDScoringResults(void);		

			
			void SetCIDScoringResults(Engine::Writers::CIDScoringResults *results) ; 			
			void WriteResults(System::String *fileName) ; 
			void ReadResults(System::String *fileName) ; 
			
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