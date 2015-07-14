// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#include "clsCIDETDScoringResults.h"
#using <mscorlib.dll>
#include <algorithm>
namespace GlypID
{
	namespace Results
	{
		clsCIDETDCombinedScoringResults::clsCIDETDCombinedScoringResults(void)
		{
			mobj_cid_etd_combined_scoring_results = NULL;
			menmFileType = GlypID::Readers::FileType::UNDEFINED ; 
		}

		clsCIDETDCombinedScoringResults::~clsCIDETDCombinedScoringResults(void)
		{
			if (mobj_cid_etd_combined_scoring_results != NULL)
			{
				delete mobj_cid_etd_combined_scoring_results ; 
				mobj_cid_etd_combined_scoring_results = NULL ; 
			}
		}

		void clsCIDETDCombinedScoringResults::SetCombinedScoringResults(Engine::Writers::CIDETDCombinedScoringResults *results)
		{
			if (mobj_cid_etd_combined_scoring_results != NULL)
			{
				delete mobj_cid_etd_combined_scoring_results ; 
				mobj_cid_etd_combined_scoring_results = NULL ; 
			}
			mobj_cid_etd_combined_scoring_results = results ; 
		}
		
		void clsCIDETDCombinedScoringResults::ReadResults(System::String *fileName)
		{
			
		}

		void clsCIDETDCombinedScoringResults::WriteResults(System::String *fileName)
		{
			char fileNameCh[512] ;
			try
			{
				GlypID::Utils::GetStr(fileName, fileNameCh) ; 
				mobj_cid_etd_combined_scoring_results->SaveResultsCombinedCIDETDScoring(fileNameCh) ; 
			}
			catch (char *mesg)
			{
				System::String *exception_msg = new System::String(mesg) ; 
				throw new System::Exception(exception_msg) ; 
			}

			
		}		

	}
}