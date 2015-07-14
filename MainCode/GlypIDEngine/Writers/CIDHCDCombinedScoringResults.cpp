// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
#include "CIDHCDCombinedScoringResults.h"
#include "../Utilities/String_Utils.h" 
#include <stdlib.h>
#include <fstream> 
#include <iostream> 
#include <time.h> 

namespace Engine
{
	namespace Writers
	{
		const int MAX_INFORMATION_BLOCK_SIZE_V1 = 512 ; 

		CIDHCDCombinedScoringResults::CIDHCDCombinedScoringResults()
		{					
			mint_num_records = 0 ; 
		}

		CIDHCDCombinedScoringResults::~CIDHCDCombinedScoringResults()
		{
			ClearCombinedScoringRecords() ; 
		}

		void CIDHCDCombinedScoringResults::ClearCombinedScoringRecords()
		{
			mvect_score_records.clear() ; 
			mint_num_records = 0 ; 
		}

		void CIDHCDCombinedScoringResults::AddCombinedScoringRecord(Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord score_record)
		{
			mvect_score_records.push_back(score_record) ; 
		}

		void CIDHCDCombinedScoringResults::AddMultipleCombinedRecords(std::vector<Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord> &vect_score_records)
		{
			if (vect_score_records.size() == 0)
				return ; 
				
			mvect_score_records.insert(mvect_score_records.end(), vect_score_records.begin(), vect_score_records.end()) ; 		
			mint_num_records += vect_score_records.size() ; 
		}	
		
		void CIDHCDCombinedScoringResults::SaveResultsCombinedScoring(char *outFile) 
		{
			char out_file_name[512] ; 
			strcpy(out_file_name, outFile) ; 	
			strcat(out_file_name, "_cid_hcd.csv") ;
			std::ofstream fout(out_file_name);

			//------------ Printing -----------//			
			fout<<"CID_Scan"<<","<<"HCD_Scan"<<","<<"Parent_Scan"<<","<<"Parent_Mz" ; 
			fout<<","<<"Mono_Mz"<<","<<"Charge_State"<<","<<"Monoisotopic_Mass"<<","<<"Isotopic_Fit"; 
			fout<<","<<"Mono_Intensity"<<","<<"HCD_Score"<<","<<"CID_score"<<","<<"CID_P_Value"<<","<<"Contains_Oxonium"<<","<<"Glycan_Type"<<","; 
			fout << "Protein_name"<<","<<"Peptide_mass"<<","<<"Sequence"<<","<<"Site"<<","<<"Glycan_Mass"<<",";
			fout<<"Glycan_Composition"<<","<<"PPM_Error"<<",";
			fout<<"Y1_ion"<<","<<"Parent_scan_time"<<","<<"CID_scan_time"<<","<<"HCD_scan_time"; 
			fout<<"\n" ; 
			fout.precision(4) ; 
			fout.setf(std::ios::fixed, std::ios::floatfield);
			Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord score_record ; 
			char type[20] ; 
			 
			Engine::MS2HCDScoring::GLYCAN_TYPE gType ; 
			for (int i = 0 ; i < (int)mvect_score_records.size() ; i++)
			{
				score_record = mvect_score_records[i] ; 
				gType = (Engine::MS2HCDScoring::GLYCAN_TYPE) score_record.menm_glycan_type ; 
			    switch (gType)
				{
						case Engine::MS2HCDScoring::GLYCAN_TYPE::HIGH_MANNOSE:
                            strcpy(type,"HIGH_MANNOSE");
                            break;
						case Engine::MS2HCDScoring::GLYCAN_TYPE::COMPLEX_ASIALYLATED:                        
                            strcpy(type, "C_ASIALYLATED");
                            break; 
						case Engine::MS2HCDScoring::GLYCAN_TYPE::COMPLEX_SIALYLATED:                                                
                            strcpy(type,"C_SIALYLATED");
                            break; 
						case Engine::MS2HCDScoring::GLYCAN_TYPE::HYBRID:
                            strcpy(type,"HYBRID");
                            break; 
						case Engine::MS2HCDScoring::GLYCAN_TYPE::NA:
                            strcpy(type,"NA");
                            break; 
                        default:
                            break;
				}	

				const char *protein_seq ; 
				const char *peptide_seq ;
				const char *glycan_comp ; 
				const char *contains_site ; 
				const char *contains_oxonium ; 
				const char *site_position ; 

				if (score_record.mbln_contains_oxonium_ions)
					contains_oxonium = "yes" ; 
				else
					contains_oxonium = "no" ;
				
				protein_seq = score_record.mstr_pro_seq_name.c_str();
				peptide_seq = score_record.mstr_pep_seq_name.c_str() ; 
				glycan_comp = score_record.mstr_glycan_composition.c_str() ; 	
				contains_site = score_record.mstr_glyco_site.c_str() ; 
				site_position = score_record.mstr_nglyco_site.c_str() ; 

				fout<<score_record.mint_cid_scan_num<<","<<score_record.mint_hcd_scan_num<<","<<score_record.mint_parent_scan_num<<","<<score_record.mdbl_parent_mz ;
				fout<<","<<score_record.mdbl_mono_mz<<","<<score_record.mshort_cs; 
				fout<<","<<score_record.mdbl_mono_mw<<","<<score_record.mdbl_fit ; 
				fout<<","<<score_record.mint_mono_intensity<<","; 
				char buffer[100] ; 
				sprintf(buffer, "%E",score_record.mdbl_hcd_score);  
				fout<<buffer<<"," ; 
				fout<<score_record.mdbl_cid_score<<",";
				fout<<score_record.mdbl_cid_score_p_value<<","<<contains_oxonium<<","<<type<<",";
				fout<< protein_seq<< "," << score_record.mdbl_seq_mass<<","<<peptide_seq<<","<<site_position<<","<<score_record.mdbl_glycan_mass<<","; 
				fout<< glycan_comp<< "," <<score_record.mdbl_ppm_error <<","; 
				fout<<score_record.mdbl_y1_mz<<","<<score_record.mdbl_parent_scan_time<<","<<score_record.mdbl_cid_scan_time<<","<<score_record.mdbl_hcd_scan_time ; 
				fout<<"\n" ; 
			}
			fout.close() ; 
		}		
		
		
		void CIDHCDCombinedScoringResults::LoadResultsCombinedScoring(char *inpFile)
		{
			std::ifstream fin(inpFile, std::ios::binary) ; 
			mvect_score_records.clear() ; 
			Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord score_record ; 

			// TO DO
		}

		
	}
}