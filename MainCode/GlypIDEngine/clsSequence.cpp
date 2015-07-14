// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University

#include ".\clssequence.h"
#include "GlypIDEngineUtils.h" 
#include "./SequenceManager/AminoAcid.h"
#include "./Utilities/system.h"

using namespace System;
using namespace System::Collections;

namespace GlypID
{
	namespace Sequence
	{
		clsSequence::clsSequence(void)
		{
		}

		clsSequence::~clsSequence(void)
		{
		}

		double clsSequence::CalculateSequenceMass(bool return_mono)
		{
			Engine::SequenceManager::Sequence sq; 
			char seq_string[512] ; 
			char seq_name[512] ; 
			seq_string[0] = '\0' ; 
			seq_name[0] = '\0' ;
			GlypID::Utils::GetStr(sequence, seq_string); 
			GlypID::Utils::GetStr(proteinName, seq_name) ; 
			sq.setSeqString(seq_string);
			sq.setSeqName(seq_name) ; 
			
			//add modification when calculating peptide mass
			Engine::SequenceManager::TAminoacids tas;
			tas['C'].resetMass(tas['C'].monoMass() + IODOACEDAMIDE_MONO, tas['C'].averageMass() + IODOACEDAMIDE_AVG); 
			
			string str = sq.getSeqString();
			float pep_mass;
			if (return_mono)
				pep_mass = (float) tas.monoPeptideMass(str);
			else
				pep_mass = (float) tas.averagePeptideMass(str) ; 
			
			return (double) pep_mass;
		}

		

		
		void clsSequence::Set(Engine::SequenceManager::Sequence &sq)
		{
			/*peptideSeq = __nogc new Engine::SequenceManager::Sequence() ; 
			peptideSeq->
			this->proteinName = sq.getName() ; */
			
			
		}
	}
}
