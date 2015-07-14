// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#using <mscorlib.dll>
#include "clsHornTransformResults.h"
#include <exception>

namespace GlypID
{
	namespace HornTransform
	{
		clsHornTransformResults::clsHornTransformResults(void)
		{
		}

		clsHornTransformResults::~clsHornTransformResults(void)
		{
			delete marr_isotope_peak_indices ; 
		}




	}
}