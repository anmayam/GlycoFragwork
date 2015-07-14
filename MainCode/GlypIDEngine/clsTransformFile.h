// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#pragma once
#include "clsHornTransformParameters.h"
#include "clsPeakProcessorParameters.h"
#include "clsTransformResults.h"
#include <fstream> 

namespace GlypID
{
	public __gc class clsTransformFile
	{
		void WritePeakParameters(std::fstream &fout, Peaks::clsPeakProcessorParameters *peak_parameters) ; 
		void WriteTransformParameters(std::fstream &fout, HornTransform::clsHornTransformParameters *transform_parameters) ;
		void WritePeaks(std::fstream &fout, GlypID::Results::clsTransformResults *results) ; 
	public:
		clsTransformFile(void);
		~clsTransformFile(void);
		void ReadFile(System::String *file_name, GlypID::Results::clsTransformResults *results, 
			Peaks::clsPeakProcessorParameters *peak_parameters, HornTransform::clsHornTransformParameters *transform_parameters) ; 
		void WriteFile(System::String *file_name, GlypID::Results::clsTransformResults *results, 
			Peaks::clsPeakProcessorParameters *peak_parameters, HornTransform::clsHornTransformParameters *transform_parameters) ; 
		void WriteFile(char *file_name, GlypID::Results::clsTransformResults *results, 
			Peaks::clsPeakProcessorParameters *peak_parameters, HornTransform::clsHornTransformParameters *transform_parameters) ; 
	};
}