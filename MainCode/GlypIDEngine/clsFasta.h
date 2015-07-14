// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#pragma once
#include ".\clsSequence.h"
#include "Readers\FastaFile.h"




#using <mscorlib.dll>

namespace GlypID
{
 namespace Readers
	{
		
		public __gc class clsFasta
		{

		public:
			 clsFasta();
			~clsFasta(void);

			
			void LoadFile(System::String *file_name, GlypID::Sequence::clsSequence* (&sequences) __gc []) ; 
			void Close() ; 
		};
	}
}