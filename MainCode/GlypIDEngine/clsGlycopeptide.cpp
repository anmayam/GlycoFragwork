#include "clsGlycopeptide.h"
#include "GlypIDEngineUtils.h" 
#include "MS2CIDScore/CIDInformationRecord.h"
#include <string>

 

namespace GlypID
{
	namespace Glycopeptide
	{
		clsGlycopeptide::clsGlycopeptide(void)
		{
			try
			{
				mobj_glycan_processor = new Engine::GlycanCompositionManager::GlycanProcessor() ; 
				mobj_fasta_file = new Engine::Readers::TFastaFile() ; 
			}
			catch (char *mesg)
			{
				System::String *exception_msg = new System::String(mesg) ; 
				throw new System::Exception(exception_msg) ; 
			}

		}

		clsGlycopeptide::~clsGlycopeptide(void)
		{
			mobj_fasta_file = NULL ; 
			mobj_fasta_file = NULL ; 
		}
		
		void clsGlycopeptide::LoadProteinsFromFasta(System::String *fFileName, GlypID::Sequence::clsSequence* (&proteins) __gc[])
		{
			std::vector <Engine::SequenceManager::Sequence> vect_candidate_proteins ; 
			bool readfasta ; 
			char fFileNameCh[512] ;

			GlypID::Utils::GetStr(fFileName, fFileNameCh) ; 
			readfasta = Engine::Readers::TFastaFile::ReadFastaFile(fFileNameCh, vect_candidate_proteins) ; 			 

			int numProteins = vect_candidate_proteins.size();
			proteins = new GlypID::Sequence::clsSequence* __gc[numProteins] ; 
			for (int i = 0; i< numProteins ; i++)
			{
				Engine::SequenceManager::Sequence sq = vect_candidate_proteins[i]; 
				proteins[i] = new GlypID::Sequence::clsSequence() ; 
				proteins[i]->sequence	= sq.getSeqString().c_str(); 
				proteins[i]->proteinName = sq.getName().c_str() ; 
			}
		}

		

		
		void clsGlycopeptide::LoadNGlycopeptidesFromFasta(System::String *fFileName, GlypID::Sequence::clsSequence* (&nglyco) __gc[], bool add_decoy_peptide)
		{
			std::vector <Engine::SequenceManager::Sequence> vect_candidate_proteins ; 
			std::vector <Engine::SequenceManager::Sequence> vect_nglyco_peps ; 
			std::vector <Engine::SequenceManager::Sequence> vect_rev_nglyco_peps ; 

			bool readfasta ; 
			char fFileNameCh[512] ;

			GlypID::Utils::GetStr(fFileName, fFileNameCh) ; 
			readfasta = Engine::Readers::TFastaFile::ReadFastaFile(fFileNameCh, vect_candidate_proteins) ; 	

			if (readfasta)
			{
				mobj_glycan_processor->GetNGlycopeptideSeqList(vect_nglyco_peps, vect_candidate_proteins) ; 

				if (add_decoy_peptide)
				{
					mobj_glycan_processor->CreateReverseGlycopeptides(vect_nglyco_peps, vect_rev_nglyco_peps) ; 
				}

				int fwd_size = vect_nglyco_peps.size() ; 
				int rev_size = vect_rev_nglyco_peps.size() ; 
				if (add_decoy_peptide)
					nglyco = new GlypID::Sequence::clsSequence* __gc[vect_nglyco_peps.size() + vect_rev_nglyco_peps.size()]; 
				else
					nglyco = new GlypID::Sequence::clsSequence* __gc[vect_nglyco_peps.size()]; 

				for (int i= 0 ; i < vect_nglyco_peps.size(); i++)
				{
					nglyco[i] = new GlypID::Sequence::clsSequence() ; 
					nglyco[i]->sequence = vect_nglyco_peps[i].getSeqString().c_str();
					nglyco[i]->proteinName = vect_nglyco_peps[i].getName().c_str() ; 
					nglyco[i]->nGlycoSite = vect_nglyco_peps[i].getSitePosition().c_str() ; 
					nglyco[i]->mass =  vect_nglyco_peps[i].seqMass ;
					nglyco[i]->is_decoy= false ; 

				}
				
				if (add_decoy_peptide)
				{				
					int index = 0 ; 
					for (int i= fwd_size ; i < fwd_size+rev_size; i++)
					{
						nglyco[i] = new GlypID::Sequence::clsSequence() ; 
						nglyco[i]->sequence = vect_rev_nglyco_peps[index].getSeqString().c_str();
						nglyco[i]->proteinName = vect_rev_nglyco_peps[index].getName().c_str() ; 
						nglyco[i]->nGlycoSite = vect_rev_nglyco_peps[index].getSitePosition().c_str() ; 					
						nglyco[i]->mass = vect_rev_nglyco_peps[index].seqMass ; 
						nglyco[i]->is_decoy = true ; 
						index++; 
					}
				}
			}
			
		}

		

		void clsGlycopeptide::GetNGlycopeptideSequences(GlypID::Sequence::clsSequence* (&proteins) __gc[], GlypID::Sequence::clsSequence* (&nglyco) __gc[])
		{
			

			std::vector <Engine::SequenceManager::Sequence> vect_candidate_proteins ;
			std::vector <Engine::SequenceManager::Sequence> vect_nglyco_peps ; 
			for (int i=0 ; i < proteins->Length; i++)
			{
				Engine::SequenceManager::Sequence sq; 
				char seq_string[1000] ; 
				char seq_name[512] ; 
				seq_string[0] = '\0' ; 
				seq_name[0] = '\0' ;
				GlypID::Utils::GetStr(proteins[i]->sequence, seq_string); 
				GlypID::Utils::GetStr(proteins[i]->proteinName, seq_name) ; 
				sq.setSeqString(seq_string);
				sq.setSeqName(seq_name) ; 
				vect_candidate_proteins.push_back(sq) ; 
			}

			mobj_glycan_processor->GetNGlycopeptideSeqList(vect_nglyco_peps, vect_candidate_proteins) ; 
			nglyco = new GlypID::Sequence::clsSequence* __gc[vect_nglyco_peps.size()]; 
			for (int i = 0 ; i < vect_nglyco_peps.size() ; i++)
			{
				nglyco[i] = new GlypID::Sequence::clsSequence() ; 
				nglyco[i]->sequence = vect_nglyco_peps[i].getSeqString().c_str();
				nglyco[i]->proteinName = vect_nglyco_peps[i].getName().c_str() ; 
				nglyco[i]->mass = vect_nglyco_peps[i].seqMass ; 			
				
			}
		}

		void clsGlycopeptide::LoadGlycansFromList(System::String *gFileName, GlypID::Glycan::clsGlycan* (&glycans) __gc[], bool add_decoy_glycan)
		{
			std::vector<Engine::GlycanTheoretical::GlycanComposition> vect_glycans;
			char gFileNameCh[512]; 

			GlypID::Utils::GetStr(gFileName, gFileNameCh) ; 
			Engine::Readers::GlycanIo *glycan_data = new Engine::Readers::GlycanIo(); 
			glycan_data->LoadGlycanListFromFile(gFileNameCh, vect_glycans) ; 

			int numGlycans = 0 ; 			
			
			if (add_decoy_glycan)
				numGlycans = 2 * (int)vect_glycans.size() ; 
			else
				numGlycans = (int) vect_glycans.size()   ; 

			
			glycans = new GlypID::Glycan::clsGlycan* __gc[numGlycans] ; 
			int index = 0 ;
			int i = 0 ; 

			

	
			while (i < (int)vect_glycans.size())
			{
				Engine::GlycanTheoretical::GlycanComposition glyc = vect_glycans[i] ; 

				glycans[index] = new GlypID::Glycan::clsGlycan() ; 
				glycans[index]->composition = glyc.getCompositionString().c_str() ; 
				glycans[index]->accurate_mass = glyc.accurateMass() ; 
				glycans[index]->average_mass = glyc.averageMass() ; 
				glycans[index]->is_decoy = false; 
				glycans[index]->SetMonosaccharideCompostion(glyc.getCompositionString().c_str()) ; 

				
				
				if(add_decoy_glycan)
				{					
					glycans[index+1] = new GlypID::Glycan::clsGlycan() ; 					
					char decoy[512] ; 
					char comp[512] ; 
					//strcpy(decoy, "Decoy_") ; 
					GlypID::Utils::GetStr(glycans[i]->composition , comp);
					glycans[index+1]->composition = comp ; 
					glycans[index+1]->accurate_mass = glyc.accurateMass() + 40  ; 
					glycans[index+1]->average_mass = glyc.averageMass() + 40 ;
					glycans[index+1]->is_decoy = true ; 
					glycans[index]->SetMonosaccharideCompostion(glyc.getCompositionString().c_str()) ; 
					index += 2; 
				}
				else
					index++ ;
				i++ ;

			}
		}

		void clsGlycopeptide::SetSearchParameters(GlypID::Scoring::clsScoringParameters *scoring_parameters)
		{
			if(scoring_parameters->get_UsePPM())
				mobj_glycan_processor->SetPPMError(scoring_parameters->get_PPMTolerance()) ; 
			else
				mobj_glycan_processor->SetDAError(scoring_parameters->get_DaTolerance()) ; 
		}

		/*void clsGlycopeptide::SetGlycopeptides(GlypID::Sequence::clsSequence* (&proteins) __gc[])
		{
			std::vector <Engine::SequenceManager::Sequence> vect_candidate_proteins ;
			for (int i=0 ; i < proteins->Length; i++)
			{
				Engine::SequenceManager::Sequence sq; 
				char seq_string[10000000] ; 
				char seq_name[512] ; 
				seq_string[0] = '\0' ; 
				seq_name[0] = '\0' ;
				GlypID::Utils::GetStr(proteins[i]->sequence, seq_string); 
				GlypID::Utils::GetStr(proteins[i]->proteinName, seq_name) ; 
				sq.setSeqString(seq_string);
				sq.setSeqName(seq_name) ; 
				vect_candidate_proteins.push_back(sq) ; 
			}
			
			std::vector<double> peptide_masses ; 
			mobj_glycan_processor->GetSequenceMassesFromProteins(vect_candidate_proteins, peptide_masses) ;
		}*/

		double clsGlycopeptide::CalculateSequenceMass(System::String *seq)
		{
			Engine::SequenceManager::Sequence sq; 
			char seq_string[512] ; 
			seq_string[0] = '\0' ; 
			GlypID::Utils::GetStr(seq, seq_string); 
			sq.setSeqString(seq_string);			
			float mass = sq.calculateMass(true) ; 
			return (double) mass ;
		}


		

		void clsGlycopeptide::SearchForGlycopeptides(GlypID::ETDScoring::clsETDScoringScanResults* (&etdResults) __gc[], GlypID::Glycan::clsGlycan* (&glycans) __gc[],
			GlypID::Sequence::clsSequence* (&nglyco_peps) __gc[], GlypID::enmGlycanType type)		
		{
			std::vector<Engine::GlycanTheoretical::GlycanComposition> vect_glycans ; 
			std::vector<Engine::SequenceManager::Sequence> vect_candidate_nglyco_peps ; 


			// -- Will have to copy from array to vectors since you can't declare vectors as part of a class in managed spce -- //
			for (int i=0 ; i < nglyco_peps->Length; i++)
			{
				
				Engine::SequenceManager::Sequence sq; 
				char seq_string[512] ; 
				char seq_name[512] ; 
				char nglyco_site[512] ; 
				seq_string[0] = '\0' ; 
				seq_name[0] = '\0' ;
				nglyco_site[0] = '\0';				
				GlypID::Utils::GetStr(nglyco_peps[i]->sequence, seq_string); 
				GlypID::Utils::GetStr(nglyco_peps[i]->proteinName, seq_name) ; 
				GlypID::Utils::GetStr(nglyco_peps[i]->nGlycoSite, nglyco_site) ; 
				sq.setSeqString(seq_string);
				sq.setSeqName(seq_name) ; 		
				sq.seqMass = nglyco_peps[i]->mass ; 
				sq.is_decoy = nglyco_peps[i]->is_decoy ; 
				sq.setNSitePosition(nglyco_site) ; 
				vect_candidate_nglyco_peps.push_back(sq) ; 

				
			}
		
			for (int i=0; i<glycans->Length ; i++)
			{
				Engine::GlycanTheoretical::GlycanComposition glyc ; 
				char gly_ch[512] ; 
				gly_ch[0] = '\0' ; 
				GlypID::Utils::GetStr(glycans[i]->composition, gly_ch); 
				glyc.initNGlycanConstants() ; 
				glyc.parseCompositionString(gly_ch); 
				glyc.glycan_accurate_mass =  glycans[i]->accurate_mass;
				glyc.glycan_avg_mass = glycans[i]->average_mass;
				glyc.is_decoy = glycans[i]->is_decoy ;
				vect_glycans.push_back(glyc) ; 
				
			
			}
			
			

			std::vector<Engine::MS2CIDScoring::CIDInformationRecord> vect_cid_record;
			Engine::MS2HCDScoring::GLYCAN_TYPE thisType = (Engine::MS2HCDScoring::GLYCAN_TYPE) type;
			for (int i = 0; i < etdResults->Length ; i++)
			{
				Engine::MS2CIDScoring::CIDInformationRecord cid_record ; 
				cid_record.AddInfoToCIDRecord(etdResults[i]->mint_msn_scan_num, etdResults[i]->mint_parent_scan_num, etdResults[i]->mint_msn_scan_level, 
				etdResults[i]->mint_parent_scan_level, etdResults[i]->mdbl_parent_mz, etdResults[i]->mdbl_mono_mz, etdResults[i]->mint_mono_intensity, etdResults[i]->mint_mono_intensity,
				etdResults[i]->mshort_cs, etdResults[i]->mdbl_fit, etdResults[i]->mdbl_mono_mw, etdResults[i]->mdbl_average_mw,
				0, 0, true, 0, etdResults[i]->mdbl_parent_scan_time,0) ; 
				
				//mobj_glycan_processor->SearchAndFillPeptideInformationSingleRecord(vect_glycans, vect_candidate_nglyco_peps, cid_record, thisType, true, false); 
				vect_cid_record.push_back(cid_record) ; 
			}		
			

			mobj_glycan_processor->SearchAndFillPeptideInformation(vect_glycans, vect_candidate_nglyco_peps, vect_cid_record, thisType, true, false); 

			etdResults->Clear();
			etdResults = new GlypID::ETDScoring::clsETDScoringScanResults* __gc[vect_cid_record.size()] ; 
			for (int i=0 ; i<vect_cid_record.size(); i++)
			{
				Engine::MS2CIDScoring::CIDInformationRecord cid_record  = vect_cid_record.at(i);
				GlypID::ETDScoring::clsETDScoringScanResults *etdResult = new GlypID::ETDScoring::clsETDScoringScanResults() ; 
				etdResult->Set(cid_record);
				etdResults[i]= new GlypID::ETDScoring::clsETDScoringScanResults() ; 
				etdResults[i] = etdResult ; 
			}
		}
	}
}