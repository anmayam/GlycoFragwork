// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath and Yin Wu, Indiana University,
#include "GlycanProcessor.h"
#include "../Utilities/system.h"
#include "../Utilities/util.h"

#include <vector>
#include <map>
namespace Engine
{
	namespace GlycanCompositionManager
	{
		GlycanProcessor::GlycanProcessor(void)
		{
			mflt_mass_error = (float) 20 ; 			
			mbln_use_ppm = true ; 
			mdbl_H2O = 18.01056; 
			mdbl_GlcNac =  203.0794; 
			mdbl_proton_mass = 1.00727638; // AM:this might not be enough precision
			mflt_ft_mz_error = (float)0.015 ; 
		}

		GlycanProcessor::~GlycanProcessor(void)
		{
			mvect_glycan_compositions.clear() ; 
			mvect_nglyco_peps.clear() ; 
		}

		void GlycanProcessor::SetPPMError(double error)
		{
			mbln_use_ppm = true ; 
			mflt_mass_error = (float) error ; 
		}

		void GlycanProcessor::SetDAError(double error)
		{
			mbln_use_ppm = false ; 
			mflt_mass_error = (float) error ; 
		}

		void GlycanProcessor::CreateReverseGlycopeptides(std::vector <Engine::SequenceManager::Sequence> &in_glyco_peps, std::vector <Engine::SequenceManager::Sequence> &out_glyco_peps)
		{
			for (int i = 0 ; i < in_glyco_peps.size(); i++)
			{
				Engine::SequenceManager::Sequence fwd_nglycopeptide = in_glyco_peps[i] ; 
				string fwd_sequence = fwd_nglycopeptide.getSeqString() ; 
				Engine::SequenceManager::Sequence rev_nglycopeptide; 

				int found = GetMotifPositionWithinSequence(fwd_sequence, 0) ; 
				while (found != -1)
				{
					string rev_sequence ; 
					rev_sequence = ReverseSequence(fwd_sequence, found) ; 
					rev_nglycopeptide.is_decoy  = true ; 
					rev_nglycopeptide.setSeqName("decoy_" + fwd_nglycopeptide.getName()) ; 
					rev_nglycopeptide.setSeqString(rev_sequence) ; 
					rev_nglycopeptide.setNSitePosition("N00"); 
					rev_nglycopeptide.seqMass = rev_nglycopeptide.calculateMass(true) ;
					out_glyco_peps.push_back(rev_nglycopeptide) ; 
					found = GetMotifPositionWithinSequence(fwd_sequence, found+1);
					if (found >0)
					{
						bool debug = true ; 
					}
				}
			}
		}

		string GlycanProcessor::ReverseSequence(string inp_sequence, int position)
		{
			int length = inp_sequence.length() ; 			
			string reverse_sequence ;
			
			char char_reverse_sequence[512] ; 
			int rindex = 1 ; 
			// AM: May 1, 2013. Switching K,R back to end. // Need to test			
			for (int i = length-2 ; i > 0; i--)
			{
				char_reverse_sequence[rindex] = inp_sequence[i] ; 
				rindex++; 
			}
			char_reverse_sequence[rindex] = '\0' ; 			

			for (int j = 0; j <3 ; j++)
			{
				int lindex = (length-position-1) - j;
				int rindex = position+j;
				if (lindex <0)
				{
					bool debug = true ; 
				}
				if (lindex >=0 && rindex <length)
					swap(char_reverse_sequence[lindex], char_reverse_sequence[rindex]) ; 
			}

			

			reverse_sequence = char_reverse_sequence ; 
			return reverse_sequence ; 

		}

		int GlycanProcessor::GetMotifPositionWithinSequence(string sequence, int search_start_position)
		{
			int pos = -1 ;		
			bool found = false ; 

			if (search_start_position >= sequence.length()-1) // To get rid of xxxNK after N has been picked up
			{
				return pos ; 
			}
			
			for (int i = search_start_position ; i < (int)sequence.length() - 2 ; i++)
			{
				if (sequence[i] == 'N')
				{
					pos = i ; 				
					if(sequence[i+2] == 'S' || sequence[i+2] == 'T')
					{
						if(sequence[i+1] != 'P')
						{
							pos = i ; 
							found = true ; 
							break ; 
						}
					}
				}
			}

			if (!found)
			{					
				// check  if it is xxxxNK.
				if ((sequence[sequence.length()-2] == 'N') && (sequence[sequence.length()-1] != 'P'))
				{
					pos = sequence.length()-2 ; 					
				}
				else
				{
					pos = -1 ; 
				}

				
			}

			return pos ; 
		}


		void GlycanProcessor::GetNGlycopeptideSeqList(std::vector <Engine::SequenceManager::Sequence> &glyco_peps, std::vector <Engine::SequenceManager::Sequence> &protein_list)
		{
			//build the peptide mass list
			//prepare for tryptic digestion
			if (!Engine::SequenceManager::ProteaseManager::registered("Trypsin"))
				Engine::SequenceManager::ProteaseManager::createProtease("RK", true, "P", "Trypsin");

			Engine::SequenceManager::Protease& trypsin = Engine::SequenceManager::ProteaseManager::getProteaseByName("Trypsin");			
			Engine::SequenceManager::Digest digest;
			
			digest.setMaxMissedCleavages(2);
			digest.setProtease(&trypsin);

			Engine::SequenceManager::IRangeLocationFilter* weightFilter = new Engine::SequenceManager::PeptideWeightFilter(500, 5000);
			Engine::SequenceManager::IRangeLocationFilter* nglycanFilter = new Engine::SequenceManager::NGlycanFilter();
			Engine::SequenceManager::AndRangeLocationFilter* filters = new Engine::SequenceManager::AndRangeLocationFilter();
			
			filters->addFilter(weightFilter);
			filters->addFilter(nglycanFilter);
			digest.setRangeLocationFilter(filters);
			
			std::vector<Engine::SequenceManager::Sequence> tmp_peps; 	
			
			for (int i=0; i< (int) protein_list.size(); i++)
			{
				string protein_name = protein_list[i].getName();	

				//generate the possible N-glyco tryptic peptides
				tmp_peps.clear();
				digest.setSequence(&(protein_list[i]));
				digest.addDigestFeatures();
				GetNGlycopeptides(protein_list[i], digest, tmp_peps);				
				glyco_peps.insert(glyco_peps.end(), tmp_peps.begin(), tmp_peps.end());
				
			}			
			
		}

		void GlycanProcessor::GetNGlycopeptides(Engine::SequenceManager::Sequence& protein_seq, Engine::SequenceManager::Digest& digest,
			std::vector<Engine::SequenceManager::Sequence>& nglyco_peps)
		{
			
			std::vector<string>& nglycopeptides = protein_seq.getAnnotation()[Engine::SequenceManager::Digest::PEPTIDE_FEATURE_TYPE];

			// below are Yin's comments
			//simulate one-time nonspecific protease cleavage. 
			//remove this code later!!!???
			//vector<string> nglycopeptides;
			//nonspecificProteaseCleavage(nglycopeptide, nglycopeptides);
			
			//filter out peptide duplications
			std::vector<string> distinct_peps;
			for ( int i=0; i < (int) nglycopeptides.size(); i++)
			{
				string& str1 = nglycopeptides[i];
				bool is_duplicated = false;
				for (int j=0; j < (int) distinct_peps.size(); j++)
				{
					string& str2 = distinct_peps[j];
					if ( str1.compare(str2) == 0 )
					{
						is_duplicated = true;
						break;
					}
				}

				if ( ! is_duplicated )
					distinct_peps.push_back(str1);
			}

			

			// Anoop April 2012 Setting glycosylation position			
			string protein = protein_seq.getSeqString() ; 
			
			for (int i=0; i < (int) distinct_peps.size() ; i++)
			{			
				string site_position ; 
				string &str3 = distinct_peps[i] ; 				
				size_t found = protein.find(str3); 
				while(found != string::npos)
				{
					bool nfound = false ; 
					for (int j= (int)found ; j < ((int)found + str3.length()-1) ; j++)
					{
						if (protein[j] == 'N')
						{
							if(protein[j+2] == 'S' || protein[j+2] == 'T')
							{
								if(protein[j+1] != 'P')
								{									
									if (!site_position.empty())
									{
										site_position = site_position + "/" ; 										
									}
									stringstream ss ;
									ss<<(j+1) ; 
									site_position = site_position +  "N"; 
									site_position = site_position + ss.str()  ; 
								}
							}
						}
					}	
					found = protein.find(str3, found+1); 
					
				}
				// Anoop July 2012 calculate mass here itself
				Engine::SequenceManager::Sequence this_seq ; 
				this_seq.setSeqName(protein_seq.getName()) ; 
				this_seq.setSeqString(distinct_peps[i]) ; 
				this_seq.setNSitePosition(site_position); 
				this_seq.seqMass = this_seq.calculateMass(true) ; 				
				nglyco_peps.push_back(this_seq);

			}


			/*for (int i=0; i< (int) distinct_peps.size(); i++)
			{
				nglyco_peps.push_back(Engine::SequenceManager::Sequence(protein_seq.getName(), distinct_peps[i]));
			}*/

			distinct_peps.clear();
		}

	

		Engine::GlycanTheoretical::GlycanComposition* GlycanProcessor::MassToGlycanComposition(float mass, float err_win,
			std::vector<Engine::GlycanTheoretical::GlycanComposition> &gc_list)
		{
			assert(err_win > 0);
			Engine::GlycanTheoretical::GlycanComposition *result = NULL;
			float min_err = err_win;
			
			for (int i=0; i< (int) gc_list.size(); i++) 
			{
				Engine::GlycanTheoretical::GlycanComposition &gc = gc_list[i];
				float err = fabs(gc.accurateMass() - mass);
				if ( err <= err_win && err < min_err )
				{ //match found
					result = &(gc_list[i]);
				}
			}
			return result;
		}

		std::string GlycanProcessor::MassToGlycanCompositionString(float mass, float err_win, std::vector<Engine::GlycanTheoretical::GlycanComposition> &gc_list)
		{
			assert(err_win > 0);
			string result = string("");
			for (int i=0; i< (int) gc_list.size(); i++) 
			{
				Engine::GlycanTheoretical::GlycanComposition &gc = gc_list[i];
				float err = fabs(gc.accurateMass() - mass);
				
				if ( err <= err_win ) 
				{ //match found
					if ( result.length() > 0 )
						result.append(string(" OR "));
					result.append(gc.getCompositionString());
				}
			}
			return result;
		}

		double GlycanProcessor::GetGlycopeptideMz(double peptide_mass, double glycan_mass, int charge )
		{
			assert( charge > 0 );
			return (peptide_mass + glycan_mass - mdbl_H2O + mdbl_proton_mass * charge)/charge;
		}

		void GlycanProcessor::GetSequenceMassesFromProteins(std::vector<Engine::SequenceManager::Sequence> &protein_sequences, std::vector <double> & sequence_masses)
		{
			// ---------- Function that cuts each protein and stores sequence masses in a vector ----------//
			std::vector<Engine::SequenceManager::Sequence> nglyco_peps;			
			GetNGlycopeptideSeqList(nglyco_peps, protein_sequences);
	
			for (int j=0; j< (int) nglyco_peps.size(); j++) 
			{
				const Engine::SequenceManager::Sequence &seq = nglyco_peps[j];
				float seq_mass = seq.seqMass ; //SequenceToMass(seq, true);
				sequence_masses.push_back((double) seq_mass) ; 				
			}
			
		}
		

		void GlycanProcessor::SearchAndFillPeptideInformationSingleRecord(std::vector<Engine::GlycanTheoretical::GlycanComposition> &glycan_composition, 
			std::vector<Engine::SequenceManager::Sequence> &nglyco_sequences, Engine::MS2HCDScoring::HCDInformationRecord &hcd_record, Engine::MS2HCDScoring::GLYCAN_TYPE &glycanType, bool filter_glycan_list)
		{
			//---------- Function that appends peptide informatin to one HCD record ---------- //
			//std::vector<float> glycan_mass_list;			
			bool found_hit = false ; 
			float seq_mass = 0 ; 
			bool check_sialylated = false ; 
			bool consider_glycan = false ; 
			Engine::GlycanTheoretical::Elemental ele = Engine::GlycanTheoretical::Elemental(string("C"));			

			if ((nglyco_sequences.size() ==0) || (glycan_composition.size() == 0))
				return ; 

			int num_MH_mass = 0;			
			num_MH_mass = glycan_composition.size();

			

			//-------- Filtering based on Glycan type--------//
			// Change July 9, 2012, Adding hybrid along with sialylated predictions to filter candidates		
			if(filter_glycan_list)
			{
				if ((glycanType == Engine::MS2HCDScoring::GLYCAN_TYPE::COMPLEX_SIALYLATED)  || (glycanType == Engine::MS2HCDScoring::GLYCAN_TYPE::HYBRID))
					check_sialylated = true ; 
				else
					check_sialylated = false ; 
			}
			
						
			// first take peptide
			for (int j = 0; j < (int) nglyco_sequences.size() ; j++)
			{
				const Engine::SequenceManager::Sequence &seq = nglyco_sequences[j];
				//float seq_mass = SequenceToMass(seq, true);	
				float seq_mass = seq.seqMass; 
				bool found_glyco_site = IdentifyIfContainsGlycoSite(seq.getSeqString()) ;

				// Add glycan combination
				for (int k = 0 ; k <num_MH_mass ; k++)
				{
					Engine::GlycanTheoretical::GlycanComposition *gc = &glycan_composition[k];

					consider_glycan = false ; 
					if (filter_glycan_list)
					{
						// Yes to filtering so see if predicted glycan type is same as this particular gc type
						if (check_sialylated == gc->hasNeuAc())
						{
							consider_glycan = true ; 
						}
						else
						{
							consider_glycan = false ; 
						}
					}
					else
					{	
						consider_glycan = true ; // No filtering, so consider them all
						
					}


					if(consider_glycan)
					{
						double glycan_mass = gc->glycan_accurate_mass ; // gc->accurateMass() ; 
							
						assert( gc != NULL );

						// make glyco peptide
					
						double glycopeptide_mass = seq_mass + glycan_mass - mdbl_H2O; 
										
						
						// Calculate error
						float err ;
						hcd_record.mbln_mass_error_in_ppm = mbln_use_ppm ; // INVARIANT : setting it here so that accurate header printing can be done
						if (mbln_use_ppm)
							err = float((glycopeptide_mass - hcd_record.mdbl_mono_mw)/(hcd_record.mdbl_mono_mw) * 1000000) ; 
						else
							err = float(glycopeptide_mass - hcd_record.mdbl_mono_mw) ; 


						if ( fabsf(err) <=  mflt_mass_error ) 
						{ 
							if (found_hit)
							{
								//already assigned a sequence, so have to compare ppm error now
								if ((float)hcd_record.mdbl_mass_error > fabsf(err))
								{								
									hcd_record.AddSearchInfoToHCDRecord(seq.getName(), (double) seq_mass, (double) glycan_mass, gc->getCompositionString(), 
										seq.getSeqString(), found_glyco_site, fabs(err), mbln_use_ppm) ; 
									found_hit = true ; 
								}
							}	
							else
							{
								//new sequence
								hcd_record.AddSearchInfoToHCDRecord(seq.getName(), (double) seq_mass, (double) glycan_mass, gc->getCompositionString(), 
										seq.getSeqString(), found_glyco_site, fabs(err), mbln_use_ppm) ; 
								found_hit = true ; 
							}
						}
					}
				}
				
			}
			
		}

		


		void GlycanProcessor::SearchAndFillPeptideInformationSingleRecord(std::vector<Engine::GlycanTheoretical::GlycanComposition> &glycan_composition, 
			std::vector<Engine::SequenceManager::Sequence> &nglyco_sequences, Engine::MS2CIDScoring::CIDInformationRecord &cid_record, Engine::MS2HCDScoring::GLYCAN_TYPE &glycanType, bool filter_glycan_list)
		{		

			//---------- Function that appends peptide informatin to one CID record, but uses HCD information ---------- //	
			
			bool found_hit = false ; 
			float seq_mass = 0 ; 
			bool check_sialylated = false ; 
			bool consider_glycan = false ; 
			Engine::GlycanTheoretical::Elemental ele = Engine::GlycanTheoretical::Elemental(string("C"));			

			//---------- GlycanList ----------//
			// Do the work of getGlycanMassList except just for external list 
			int num_MH_mass = 0;			
			num_MH_mass = glycan_composition.size();

			if ((nglyco_sequences.size() ==0) || (glycan_composition.size() == 0))
				return ; 

			//-------- Filtering based on Glycan type--------//
			// July 2012 Update : inclusing hybrid structures
			if (glycanType == Engine::MS2HCDScoring::GLYCAN_TYPE::NA)
				filter_glycan_list = false ; 
			if(filter_glycan_list)
			{
				if ((glycanType == Engine::MS2HCDScoring::GLYCAN_TYPE::COMPLEX_SIALYLATED) || (glycanType == Engine::MS2HCDScoring::GLYCAN_TYPE::HYBRID))
					check_sialylated = true ; 
				else
					check_sialylated = false ; 
			}
			
			
			// first take peptide
			for (int j = 0; j < (int) nglyco_sequences.size() ; j++)
			{
				const Engine::SequenceManager::Sequence &seq = 	nglyco_sequences[j];
				float seq_mass = seq.seqMass ; //SequenceToMass(seq, true);	
				bool found_glyco_site = IdentifyIfContainsGlycoSite(seq.getSeqString()) ;
				
				// Add glycan combination
				for (int k = 0 ; k <num_MH_mass ; k++)
				{
					Engine::GlycanTheoretical::GlycanComposition *gc = &glycan_composition[k];

					string str = seq.getSeqString() ; 							

					consider_glycan = false ; 
					if (filter_glycan_list)
					{
						// Yes to filtering so see if predicted glycan type is same as this particular gc type
						if (check_sialylated == gc->hasNeuAc())
						{
							consider_glycan = true ; 
						}
						else
						{
							consider_glycan = false ; 
						}
					}
					else
					{	
						consider_glycan = true ; // No filtering, so consider them all
						
					}


					if(consider_glycan)
					{
						double glycan_mass = gc->glycan_accurate_mass; //gc->accurateMass() ; 
							
						assert( gc != NULL );

					
						double glycopeptide_mass = seq_mass + glycan_mass - mdbl_H2O; 
										
						
						// Calculate error
						float err ;
						cid_record.mbln_mass_error_in_ppm = mbln_use_ppm ; // INVARIANT : setting it here so that accurate header printing can be done
						if (mbln_use_ppm)
							err = float((glycopeptide_mass - cid_record.mdbl_mono_mw)/(cid_record.mdbl_mono_mw) * 1000000) ; 
						else
							err = float(glycopeptide_mass - cid_record.mdbl_mono_mw) ; 


						if ( fabsf(err) <=  mflt_mass_error ) 
						{
							if (found_hit)
							{
								//already assigned a sequence, so have to compare ppm error now
								if ((float)cid_record.mdbl_mass_error > fabsf(err))
								{								
									cid_record.AddSearchInfoToCIDRecord(seq.getName(), (double) seq_mass, (double) glycan_mass, gc->getCompositionString(), 
										seq.getSeqString(), seq.getSitePosition() , found_glyco_site, fabs(err), mbln_use_ppm) ; 
									found_hit = true ; 
								}
							}	
							else
							{
								//new sequence
								cid_record.AddSearchInfoToCIDRecord(seq.getName(), (double) seq_mass, (double) glycan_mass, gc->getCompositionString(), 
									seq.getSeqString(), seq.getSitePosition(), found_glyco_site, fabs(err), mbln_use_ppm) ; 
								found_hit = true ; 
							}
						}			

						
					}
				}
				
			}
			
		}


		void GlycanProcessor::SearchAndFillPeptideInformation(std::vector<Engine::GlycanTheoretical::GlycanComposition> &glycan_composition, 
				std::vector<Engine::SequenceManager::Sequence> &nglyco_sequences, 
				std::vector<Engine::MS2CIDScoring::CIDInformationRecord> &score_records, Engine::MS2HCDScoring::GLYCAN_TYPE &glycanType, 
				bool filter_glycan_list, bool include_only_best_hit) 
		{
			// ---------- Function that appends protein information to each CID record ----------//	
		
			std::vector<Engine::MS2CIDScoring::CIDInformationRecord> cid_records ; 			
			bool found_hit = false ; 
			bool check_sialylated = false ; 
			bool check_highmannose = false ; 
			Engine::GlycanTheoretical::Elemental ele = Engine::GlycanTheoretical::Elemental(string("C"));
			Engine::MS2CIDScoring::CIDInformationRecord this_record ; 	

			if ((nglyco_sequences.size() ==0) || (glycan_composition.size() == 0))
				return ; 

			int num_MH_mass = glycan_composition.size();
					


			//-------- Filtering based on Glycan type--------//
			// Change July 9, 2012, Adding hybrid along with sialylated predictions to filter candidates
			if (glycanType == Engine::MS2HCDScoring::GLYCAN_TYPE::HYBRID)
			{
				bool debug = true ; 
			}
			if (glycanType == Engine::MS2HCDScoring::GLYCAN_TYPE::NA)
				filter_glycan_list = false ; 
			if(filter_glycan_list)
			{
				if ((glycanType == Engine::MS2HCDScoring::GLYCAN_TYPE::COMPLEX_SIALYLATED) ||(glycanType == Engine::MS2HCDScoring::GLYCAN_TYPE::HYBRID))
					check_sialylated = true ; 
				else
					check_sialylated = false ; 				
			}

			
			


			//---------- Start processing every record  by looking through all possible peptide + glycan combination ---------- //
			for (int i = 0 ; i < (int) score_records.size() ; i++)
			{
				this_record = score_records[i] ; 
				found_hit = false ; 			

				// first take peptide
				for (int j = 0; j < (int) nglyco_sequences.size() ; j++)
				{
					const Engine::SequenceManager::Sequence &seq = nglyco_sequences[j];
					float seq_mass = seq.seqMass ; //SequenceToMass(seq, true);	
					bool found_glyco_site =  IdentifyIfContainsGlycoSite(seq.getSeqString()) ;

					
					assert (seq_mass !=0) ; 

					if (seq.getSeqString() == "KNGTLILLHAFHEQGGVYRSIRDTE")
					{
						bool debug = true; 
					}
					if (seq.getSeqString() == "RPTGEVYDIEIDTLETTCHVLDPTPLANCSVR")
					{
						bool debug = true; 
					}



					// Add glycan combination
					for (int k = 0 ; k <num_MH_mass ; k++)
					{
						Engine::GlycanTheoretical::GlycanComposition *gc = &glycan_composition[k];
						double glycan_mass = 0 ; 
						
									
						glycan_mass = gc->glycan_accurate_mass ; //->accurateMass() ; 

						assert (glycan_mass != 0) ; 
						
						assert( gc != NULL );


						if (abs(glycan_mass - 2604.93042) < 0.001)
						{
							bool devug = true ; 
						}
						
				

						bool consider_glycan = false ; 
						if (filter_glycan_list && (!gc->is_decoy) && (!seq.is_decoy))
						{
							// Yes to filtering so see if predicted glycan type is same as this particular gc type							
							if (check_sialylated == gc->hasNeuAc())
							{
								consider_glycan = true ; 
							}
							else
							{
								consider_glycan = false ; 
							}
						}
						else
						{	
							consider_glycan = true ; // No filtering or is a decoy so consider them all
							
						}


						if(consider_glycan)
						{
							
							double glycopeptide_mass = seq_mass + glycan_mass - mdbl_H2O; 
							
							
							// Calculate error
							float err ;
							if (mbln_use_ppm)
								err = float((glycopeptide_mass - this_record.mdbl_mono_mw)/(this_record.mdbl_mono_mw) * 1000000) ; 
							else
								err = float(glycopeptide_mass - this_record.mdbl_mono_mw) ; 

							if ( fabsf(err) <=  mflt_mass_error ) 
							{ 
								
								std::string comp ; 
								if (gc->is_decoy)			// Make sure it is reported as decoy								
								{
									comp = "Decoy_" ; 
									comp.append(gc->getCompositionString());
								}
								else
									comp = gc->getCompositionString() ; 

								if (include_only_best_hit)
								{
									if (found_hit)
									{
										//already assigned a sequence, so have to compare ppm error now
										if ((float)this_record.mdbl_mass_error > fabsf(err))
										{
											cid_records.pop_back() ;
											this_record.AddSearchInfoToCIDRecord(seq.getName(), (double) seq_mass, (double) glycan_mass, comp, 
												seq.getSeqString(), seq.getSitePosition(), found_glyco_site, fabs(err), mbln_use_ppm) ; 							
											cid_records.push_back(this_record) ; 
											found_hit = true ; 
										}
									}	
									else
									{
										this_record.AddSearchInfoToCIDRecord(seq.getName(), (double) seq_mass, (double) glycan_mass, comp, 
											seq.getSeqString(), seq.getSitePosition(), found_glyco_site, fabs(err), mbln_use_ppm) ;  							
										cid_records.push_back(this_record) ; 
										found_hit = true ; 
									}
								}
								else
								{
									// Include all hits ; will use some other method ( e.g. ETD )to filter out resutls
									if (found_hit)
									{
										bool debug = true; 
									}
									
									this_record.AddSearchInfoToCIDRecord(seq.getName(), (double) seq_mass, (double) glycan_mass, comp, 
										seq.getSeqString(), seq.getSitePosition(),  found_glyco_site, fabs(err), mbln_use_ppm) ;  							
									cid_records.push_back(this_record) ; 
									found_hit = true ; 
									
								}
							}	
						}
					}
				}
				
				if (!found_hit)
					cid_records.push_back(this_record) ; 
			}

		

			//---------- Store and return ----------//
			score_records.clear() ; 
			score_records.insert(score_records.begin(), cid_records.begin(), cid_records.end()); 


		}

		bool GlycanProcessor::IdentifyIfContainsGlycoSite(std::string peptide_sequence)
		{
			//--------- Function that scans a peptide sequence and looks for N-glycosylation motif ---------//
			char *seq = new char[peptide_sequence.size() + 1]; 
			seq[peptide_sequence.size()] = 0 ; 
			memcpy(seq, peptide_sequence.c_str(), peptide_sequence.size()) ;
			bool found = false ; 
			for (int i = 0 ; i < peptide_sequence.size()-2 ; i++)
			{
				if (seq[i] == 'N')
				{
					if(seq[i+2] == 'S' || seq[i+2] == 'T')
					{
						if(seq[i+1] != 'P')
						{
							found = true ; 
							break ; 
						}
					}
				}
			}
			return found ; 		
		}

		

		/*void GlycanProcessor::SearchAndFillPeptideInformation(std::vector<Engine::GlycanTheoretical::GlycanComposition> &glycan_composition, 
				std::vector<Engine::SequenceManager::Sequence> &protein_sequences, 
				std::vector<Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord> &score_records, bool filter_list) 
		{
			// ---------- Function that appends protein information to each CID + HCD (combined) record ----------//

			std::vector<float> glycan_mass_list;
			std::vector<Engine::SequenceManager::Sequence> nglyco_peps;
			std::vector<Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord> all_records ; 			
			bool found_hit = false ; 
			Engine::GlycanTheoretical::Elemental ele = Engine::GlycanTheoretical::Elemental(string("C"));
			Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord this_record ; 	


			//---------- GlycanList ----------//
			// Do the work of getGlycanMassList except just for external list 
			int num_MH_mass = 0;			
			num_MH_mass = glycan_composition.size();
			
		
			//---------- PeptideList ------------//			
			GetNGlycopeptideSeqList(nglyco_peps, protein_sequences);


			//---------- Start processing every record  by looking through all possible peptide + glycan combination ---------- //
			for (int i = 0 ; i < (int) score_records.size() ; i++)
			{
				this_record = score_records[i] ; 

				
				found_hit = false ; 
				// first take peptide
				for (int j = 0; j < (int) nglyco_peps.size() ; j++)
				{
					const Engine::SequenceManager::Sequence &seq = nglyco_peps[j];
					float seq_mass = SequenceToMass(seq, true);	
					bool found_glyco_site = IdentifyIfContainsGlycoSite(seq.getSeqString()) ;

					// Add glycan combination
					for (int k = 0 ; k <num_MH_mass ; k++)
					{
						Engine::GlycanTheoretical::GlycanComposition *gc = &glycan_composition[k];
						double glycan_mass = gc->accurateMass() ; 
						
						assert( gc != NULL );

						// make glyco peptide						
						//string formula = ele.glycoPepToFormula(glycopep);	
						//float highest_isotope_idx = (float) ele.highestIsotopeOfFormula(1, formula); // not sure why these are here but still keeping them
						//float second_isotope_idx = (float) ele.highestIsotopeOfFormula(2, formula);
						//assert( highest_isotope_idx >= 0 && second_isotope_idx >= 0);
						//double glycopeptide_mass1 = ele.formulaToMass(formula) - mdbl_H2O ; 
						double glycopeptide_mass = seq_mass + glycan_mass - mdbl_H2O; 
						
						

						// Calculate error
						// INVARIANT : this cannot be DA, has to be be ppm. 
						float err ;
						err = float((glycopeptide_mass - this_record.mdbl_mono_mw)/(this_record.mdbl_mono_mw) * 1000000) ; 
						
						if ( fabsf(err) <=  mflt_mass_error ) 
						{ 
							if (found_hit)
							{
								//already assigned a sequence, so have to compare ppm error now
								if ((float)this_record.mdbl_ppm_error > Engine::Utilities::roundf(fabs(err)))
								{
									this_record.AddSearchInfoToCombinedRecord(seq.getName(), (double) seq_mass, (double) glycan_mass, gc->getCompositionString(), 
										seq.getSeqString(), seq.getSitePosition(), Engine::Utilities::roundf(fabs(err)), found_glyco_site) ; 
									all_records.push_back(this_record) ; 
									found_hit = true ; 
								}
							}	
							else
							{
								//new sequence
								this_record.AddSearchInfoToCombinedRecord(seq.getName(), (double) seq_mass, (double) glycan_mass, gc->getCompositionString(), 
										seq.getSeqString(),seq.getSitePosition(), Engine::Utilities::roundf(fabs(err)), found_glyco_site) ; 
								all_records.push_back(this_record) ; 
								found_hit = true ; 
							}
						}						
					}
				}
				
				if (!found_hit)
					all_records.push_back(this_record) ; 
			}			

			//---------- Store and return ----------//
			score_records.clear() ; 
			score_records.insert(score_records.begin(), all_records.begin(), all_records.end()); 
		}*/
	}
}

