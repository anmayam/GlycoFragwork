// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University

#include ".\clsglycan.h"
#include "GlycanTheoretical\Nglycan.h"
#include "GlypIDEngineUtils.h" 


namespace GlypID
{
	namespace Glycan
	{
		clsGlycan::clsGlycan(void)
		{		
			accurate_mass = 0 ; 
			average_mass = 0 ; 
			is_decoy = false ; 
			numHex = 0 ; 
			numDeHex = 0 ; 
			numHexNAc = 0 ; 
			numNeuAc = 0 ; 
			
		}

		clsGlycan::~clsGlycan(void)
		{			
		}



		void clsGlycan::SetMonosaccharideCompostion(System::String *comp) 
		{
			if (comp->Length == 0)
			{
				numHex = -1 ; 
				numDeHex = -1 ; 
				numHexNAc = -1 ; 
				numNeuAc = -1 ; 
			}	
			Engine::GlycanTheoretical::GlycanComposition glyc; 
			char gly_ch[512] ; 
			gly_ch[0] = '\0' ; 
			GlypID::Utils::GetStr(comp, gly_ch); 
			glyc.initNGlycanConstants() ; 
			glyc.parseCompositionString(gly_ch); 

			numHex = glyc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_HEX] ; 
			numDeHex = glyc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_DEHEX] ; 
			numHexNAc = glyc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_HEXNAC] ; 
			numNeuAc = glyc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_NEUAC] ; 
		}

		double clsGlycan::CalculateGlycanMass(bool accurate)
		{
			if (composition->Length == 0)
			{
				return false ; 
			}	

			Engine::GlycanTheoretical::GlycanComposition glyc; 
			char gly_ch[512] ; 
			gly_ch[0] = '\0' ; 
			GlypID::Utils::GetStr(composition, gly_ch); 
			glyc.initNGlycanConstants() ; 
			glyc.parseCompositionString(gly_ch); 
			if (accurate)
				return ((double) glyc.accurateMass()) ; 
			else
				return ((double) glyc.averageMass()) ; 

		}

		bool clsGlycan::GlycanCompositionHasNeuAC(System::String *comp) 
		{
			if (comp->Length == 0)
			{
				return false ; 
			}	

			Engine::GlycanTheoretical::GlycanComposition glyc; 
			char gly_ch[512] ; 
			gly_ch[0] = '\0' ; 
			GlypID::Utils::GetStr(comp, gly_ch); 
			glyc.initNGlycanConstants() ; 
			glyc.parseCompositionString(gly_ch); 
			return (glyc.hasNeuAc()) ; 
		}

		bool clsGlycan::GlycanCompositionHasDeHex(System::String *comp)
		{
			if (comp->Length ==0)
			{
				return false ; 
			}

			Engine::GlycanTheoretical::GlycanComposition glyc ; 
			char gly_ch[512] ; 
			gly_ch[0] = '\0' ; 
			GlypID::Utils::GetStr(comp, gly_ch); 
			glyc.initNGlycanConstants() ; 
			glyc.parseCompositionString(gly_ch); 
			return (glyc.hasDeHex()) ; 
			
		}


	}
}
