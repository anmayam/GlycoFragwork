// Written by Anoop Mayampurath, Indiana University

#pragma once
#using <mscorlib.dll>
#include "MS2HCDScore/HCDInformationRecord.h"

namespace GlypID
{
	namespace Glycan
	{
		public __gc class clsGlycan
		{
		public:
			clsGlycan(void) ; 
			~clsGlycan(void) ; 

			double accurate_mass ; 
			double average_mass ; 
			System::String *composition ; 
			bool is_decoy ; 

			bool GlycanCompositionHasNeuAC(System::String *comp) ; 
			double CalculateGlycanMass(bool accurate); 
			bool GlycanCompositionHasDeHex(System::String *comp) ; 
			void SetMonosaccharideCompostion(System::String *comp) ; 
			int numHex ; 
			int numHexNAc ; 
			int numDeHex ; 
			int numNeuAc ; 


			

		};
	}
}

