// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University

#include "clsFasta.h"
#include "GlypIDEngineUtils.h"

#using <mscorlib.dll>
#include <vector>

namespace GlypID
{
	namespace Readers
	{
		clsFasta::clsFasta(void)
		{	
		}

		clsFasta::~clsFasta(void)
		{
		}

		void clsFasta::LoadFile(System::String __gc *file_name, GlypID::Sequence::clsSequence *(__gc &sequences) __gc[])
		{
			/*std::vector<Engine::SequenceManager::Sequence> vect_candidate_proteins ;
			char fasta_name_ch[256] ; 
			GlypID::Utils::GetStr(file_name, fasta_name_ch) ;
			bool readfasta ; 
			readfasta = Engine::Readers::TFastaFile::ReadFastaFile(fasta_name_ch, vect_candidate_proteins );  
			int numSeq = vect_candidate_proteins.size() ; 


			sequences = new GlypID::Sequence::clsSequence* __gc [numSeq] ; 
			for (int sqNum = 0 ; sqNum < numSeq ; sqNum++)
			{
				sequences[sqNum] = new GlypID::Sequence::clsSequence() ; 
				Engine::SequenceManager::Sequence &seq = vect_candidate_proteins.at[sqNum]
				sequences[sqNum]->Set(& ) ; 				
			}*/
		}

		void clsFasta::Close()
		{
		}

		
	}
}
