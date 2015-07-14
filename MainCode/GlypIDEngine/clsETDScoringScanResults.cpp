// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University

#using <mscorlib.dll>
#include "clsETDScoringScanResults.h"
#include <exception>

namespace GlypID
{
	namespace ETDScoring
	{
		clsETDScoringScanResults::clsETDScoringScanResults(void)
		{
			mbln_true_hit = false ; 	
			menm_glycan_type = Engine::MS2HCDScoring::GLYCAN_TYPE::NA ; 
			mdbl_hcd_score = 1; 
		}

		clsETDScoringScanResults::~clsETDScoringScanResults(void)
		{
			
		}
	}
}