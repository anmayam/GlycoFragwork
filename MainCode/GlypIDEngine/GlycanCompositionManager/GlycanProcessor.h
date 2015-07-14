// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath and Yin Wu, Indiana University

#include "../Readers/FastaFile.h" 
#include "../Readers/GlycanIo.h" 
#include "../MS2CIDHCDCombinedScore/CIDHCDCombinedInformationRecord.h"
#include "../MS2CIDScore/CIDInformationRecord.h" 
#include "../GlycanTheoretical/Elemental.h"
#include "../SequenceManager/Sequence.h" 
#include "../SequenceManager/ProteaseManager.h"
#include "../SequenceManager/Digest.h"
#include "../SequenceManager/Protease.h"
#include "../SequenceManager/RangeLocationFilter.h"
#include "../SequenceManager/RangeLocation.h"
#include "../SequenceManager/AminoAcid.h"

namespace Engine
{
	namespace GlycanCompositionManager
	{
		class GlycanProcessor
		{
			float mflt_mass_error ; 			
			bool mbln_use_ppm ; 
			std::vector<Engine::GlycanTheoretical::GlycanComposition> mvect_glycan_compositions ; 
			std::vector<Engine::SequenceManager::Sequence> mvect_nglyco_peps;

			double mdbl_GlcNac ; 

			double mdbl_H2O; 
			double mdbl_proton_mass ;  
			float mflt_ft_mz_error ; 
			std::vector <float> mflt_seq_masses ; 

		public:
			//! default constructor
			GlycanProcessor() ; 
			//! destructor.
			~GlycanProcessor() ;
			//! search records for peptide ids and append the information (for cid/hcd combined)
			/*void SearchAndFillPeptideInformation(std::vector<Engine::GlycanTheoretical::GlycanComposition> &glycan_composition, 
				std::vector<Engine::SequenceManager::Sequence> &nglyco_sequences, 
				std::vector<Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord> &score_records, bool filter_list) ;*/
			//! search recoreds for peptide ids and append information (for CID)
			void GlycanProcessor::SearchAndFillPeptideInformation(std::vector<Engine::GlycanTheoretical::GlycanComposition> &glycan_composition, 
				std::vector<Engine::SequenceManager::Sequence> &nglyco_sequences, 
				std::vector<Engine::MS2CIDScoring::CIDInformationRecord> &score_records, Engine::MS2HCDScoring::GLYCAN_TYPE &glycanType, 
				bool filter_glycan_list, bool include_only_best_hit) ;
			//! Search one single record for peptide ids and append information (for CID)
			void SearchAndFillPeptideInformationSingleRecord(std::vector<Engine::GlycanTheoretical::GlycanComposition> &glycan_composition, 
			std::vector<Engine::SequenceManager::Sequence> &nglyco_sequences, Engine::MS2CIDScoring::CIDInformationRecord &cid_record, Engine::MS2HCDScoring::GLYCAN_TYPE &glycanType, bool filter_glycan_list) ; //, bool include_best_hit) ;
			//! Search one single record for peptide ids aand append information (for HCD) 
			void SearchAndFillPeptideInformationSingleRecord(std::vector<Engine::GlycanTheoretical::GlycanComposition> &glycan_composition, 
			std::vector<Engine::SequenceManager::Sequence> &nglyco_sequences, Engine::MS2HCDScoring::HCDInformationRecord &hcd_record, Engine::MS2HCDScoring::GLYCAN_TYPE &glycanType, bool filter_list) ; 

			//! set ppm error
			void SetPPMError(double error) ; 
			//! set Da error
			void SetDAError(double error) ; 
			//! Process protein list and get sequence masses
			void GetSequenceMassesFromProteins(std::vector<Engine::SequenceManager::Sequence> &protein_sequences, std::vector <double> & sequence_masses);			
			//! Get Nglyocpeptides from protein list
			void GetNGlycopeptideSeqList(std::vector <Engine::SequenceManager::Sequence> &glyco_peps, std::vector <Engine::SequenceManager::Sequence> &protein_list) ; 
			//! Get Nglycopeptide from single protein
			void GetNGlycopeptides(Engine::SequenceManager::Sequence& protein_seq, Engine::SequenceManager::Digest& digest, std::vector<Engine::SequenceManager::Sequence>& nglyco_peps) ; 
			// ! Creates a set of reverse peptides given a set of forward peptides
			void CreateReverseGlycopeptides(std::vector <Engine::SequenceManager::Sequence> &in_glyco_peps, std::vector <Engine::SequenceManager::Sequence> &out_glyco_peps) ; 		
			//! Get composition given a mass
			Engine::GlycanTheoretical::GlycanComposition* MassToGlycanComposition(float mass, float err_win, std::vector<Engine::GlycanTheoretical::GlycanComposition> &gc_list) ; 
			//! String version of above
			string MassToGlycanCompositionString(float mass, float err_win, std::vector<Engine::GlycanTheoretical::GlycanComposition> &gc_list) ; 
			//! Get MZ of glycopeptide 
			double GetGlycopeptideMz(double peptide_mass, double glycan_mass, int charge ); 
			//! Identify if peptide has a glycosylation site
			bool GlycanProcessor::IdentifyIfContainsGlycoSite(std::string peptide_sequence) ; 
			//! Reverses a sequence but keeps the motif at position intact
			string ReverseSequence(string inp_sequence, int position);
			//! returns the first motif position after search_start_position
			int GetMotifPositionWithinSequence(string sequence, int search_start_position);
		};
	}
}