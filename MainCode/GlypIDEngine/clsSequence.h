// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)



#pragma once
#using <mscorlib.dll>
#include "SequenceManager/Sequence.h"

namespace GlypID
{
	namespace Sequence
	{
		public __gc class clsSequence
		{
			public:
			//! sequence string
			System::String *sequence;
			// ! Protein name
			System::String *proteinName;
			// ! Site name
			System::String *nGlycoSite ; 
			// ! O glyco sites
			System::String *oGlycoSites ; 
			//! sequence mass
			double mass ; 
			//! decoy sequence
			bool is_decoy ; 
			
			clsSequence(void);
			~clsSequence(void);
			void Set(Engine::SequenceManager::Sequence &sk) ;
			double CalculateSequenceMass(bool return_mono); 	
			
			
		};
	}
}
