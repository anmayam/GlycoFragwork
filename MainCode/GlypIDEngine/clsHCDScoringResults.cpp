// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#include "clsHCDScoringResults.h"
#using <mscorlib.dll>
#include <algorithm>
namespace GlypID
{
	namespace Results
	{
		clsHCDScoringResults::clsHCDScoringResults(void)			
		{
			mobj_hcd_scoring_results = NULL ; 			
			menmFileType = GlypID::Readers::FileType::UNDEFINED ; 
		}

		clsHCDScoringResults::~clsHCDScoringResults(void)
		{
			if (mobj_hcd_scoring_results != NULL)
			{
				delete mobj_hcd_scoring_results ; 
				mobj_hcd_scoring_results = NULL ; 
			}
		}

		void clsHCDScoringResults::SetHCDScoringResults(Engine::Writers::HCDScoringResults *results)
		{
			if (mobj_hcd_scoring_results != NULL)
			{
				delete mobj_hcd_scoring_results ; 
				mobj_hcd_scoring_results = NULL ; 
			}
			mobj_hcd_scoring_results = results ; 
		}
		
		void clsHCDScoringResults::ReadResults(System::String *fileName)
		{
			
		}

		void clsHCDScoringResults::WriteResults(System::String *fileName)
		{
			char fileNameCh[512] ;
			try
			{
				GlypID::Utils::GetStr(fileName, fileNameCh) ; 
				mobj_hcd_scoring_results->SaveResultsHCD(fileNameCh) ; 
			}
			catch (char *mesg)
			{
				System::String *exception_msg = new System::String(mesg) ; 
				throw new System::Exception(exception_msg) ; 
			}

			
		}		

	}
}