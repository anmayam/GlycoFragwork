// Written by Anoop Mayampurath, Indiana University
#pragma once
#using <mscorlib.dll>
#define ARR_SIZE 2 ; 
using namespace System ;
#include <vector>
#include "SequenceManager/Sequence.h"
#include "GlycanCompositionManager/GlycanProcessor.h"
#include "clsSequence.h"
#include "clsGlycan.h"
#include "clsScoringParameters.h"
#include "clsETDScoringScanResults.h"
#include "clsPeak.h" 
#include "clsHCDScoring.h"


namespace GlypID
{
	namespace Glycopeptide
	{
		public __gc class clsGlycopeptide
		{
		
			
			Engine::GlycanCompositionManager::GlycanProcessor __nogc *mobj_glycan_processor ; 			
			Engine::Readers::TFastaFile __nogc *mobj_fasta_file ; 
			
		

			public:		
			clsGlycopeptide(void);
			~clsGlycopeptide(void);			
			void SetSearchParameters(GlypID::Scoring::clsScoringParameters *scoring_parameters) ; 
			void LoadProteinsFromFasta(System::String *fFileName, GlypID::Sequence::clsSequence* (&proteins) __gc[]); 
			void LoadNGlycopeptidesFromFasta(System::String *fFileName, GlypID::Sequence::clsSequence* (&nglyco_peptides) __gc[], bool add_decoy_peptide);
			void LoadGlycansFromList(System::String *gFileName, GlypID::Glycan::clsGlycan* (&glycans) __gc[], bool add_decoy_glycan) ; 
			//void SetGlycopeptides(GlypID::Sequence::clsSequence* (&proteins) __gc[]) ; 	
			
			void SearchForGlycopeptides(GlypID::ETDScoring::clsETDScoringScanResults* (&etdResults) __gc[], GlypID::Glycan::clsGlycan* (&glycans) __gc[], GlypID::Sequence::clsSequence* (&nglyco_peps) __gc[], GlypID::enmGlycanType type);			
			void GetNGlycopeptideSequences(GlypID::Sequence::clsSequence* (&proteins) __gc[], GlypID::Sequence::clsSequence* (&nglyco) __gc[]);
			double CalculateSequenceMass(System::String* sequence) ; 

			
			




		};
	}
}