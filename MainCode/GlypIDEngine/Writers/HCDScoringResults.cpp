// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
#include "HCDScoringResults.h"
#include <fstream> 
#include <iostream> 
#include <time.h> 

namespace Engine
{
	namespace Writers
	{
		const int MAX_INFORMATION_BLOCK_SIZE_V1 = 512 ; 

		HCDScoringResults::HCDScoringResults()
		{		
			
			mint_num_records = 0 ; 
		}

		HCDScoringResults::~HCDScoringResults()
		{
		
			ClearHCDRecords() ; 
			
		}

		void HCDScoringResults::ClearHCDRecords()
		{
			mvect_hcd_records.clear() ; 
			mint_num_records = 0 ; 
		}

		void HCDScoringResults::AddHCDRecord(Engine::MS2HCDScoring::HCDInformationRecord hcd_record)
		{
			mvect_hcd_records.push_back(hcd_record) ; 
		}

		void HCDScoringResults::AddMultipleHCDRecords(std::vector<Engine::MS2HCDScoring::HCDInformationRecord> &vect_hcd_records)
		{
			
				if (vect_hcd_records.size() == 0)
					return ; 

				mvect_hcd_records.insert(mvect_hcd_records.end(), vect_hcd_records.begin(), vect_hcd_records.end()) ; 
				mint_num_records += vect_hcd_records.size() ; 
		}

	

		
		void HCDScoringResults::SaveResultsHCD(char *hcdFile) 
		{
			char hcd_file_name[512] ; 
			strcpy(hcd_file_name, hcdFile) ; 
			strcat(hcd_file_name, "_hcd.csv") ;
			std::ofstream fout(hcd_file_name);

			fout<<"MSn_Scan"<<","<<"MSn_Level"<<","<<"Parent_Scan"<<","<<"Parent_Scan_Level"<<","<<"Parent_Mz" ; 
			fout<<","<<"Mono_Mz"<<","<<"Charge_State"<<","<<"Monoisotopic_Mass"<<","<<"Isotopic_Fit"; 
			fout<<","<<"Mono_Intensity"<<","<<"HCD_Score"<<","<<"Glycan_Type"; 
			fout << "Protein_name"<<","<<"Peptide_mass"<<","<<"Sequence"<<","<<"Contains_Site"<<","<<"Glycan_Mass"<<",";
			fout<<"Glycan_Composition"<<"," ; 
			if (mvect_hcd_records[0].mbln_mass_error_in_ppm) 
				fout<<"Error(ppm)"<<",";
			else
				fout<<"Error(Da)"<<"," ; 
			fout<<"Y1_ion"<<","<<"Parent_scan_time"; 
			fout<<"\n" ; 
			
		
			Engine::MS2HCDScoring::HCDInformationRecord hcd_record ; 
			char type[20] ; 
			Engine::MS2HCDScoring::GLYCAN_TYPE gType ; 

			for (int i = 0 ; i < mvect_hcd_records.size() ; i++)
			{
				hcd_record = mvect_hcd_records[i] ; 
				const char *protein_seq ; 
				const char *peptide_seq ;
				const char *glycan_comp ; 
				const char *contains_site ;
				protein_seq = hcd_record.mstr_pro_seq_name.c_str();
				peptide_seq = hcd_record.mstr_pep_seq_name.c_str() ; 
				glycan_comp = hcd_record.mstr_glycan_composition.c_str() ; 
				contains_site = hcd_record.mstr_glyco_site.c_str() ; 

				gType = (Engine::MS2HCDScoring::GLYCAN_TYPE) hcd_record.menm_glycan_type ; 
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


				
				fout.precision(4) ; 
				fout.setf(std::ios::fixed, std::ios::floatfield);
				fout<<hcd_record.mint_msn_scan_num<<","<<hcd_record.mint_msn_scan_level<<","<<hcd_record.mint_parent_scan_num ; 
				fout<<","<<hcd_record.mint_parent_scan_level<<","<<hcd_record.mdbl_parent_mz ; 
				fout<<","<<hcd_record.mdbl_mono_mz<<","<<hcd_record.mshort_cs; 
				fout<<","<<hcd_record.mdbl_mono_mw<<","<<hcd_record.mdbl_fit ; 
				fout<<","<<hcd_record.mint_mono_intensity<<","; 
				//fout.setf(std::ios::floatfield) ; 
				char buffer[100] ; 
				sprintf(buffer, "%E",hcd_record.mdbl_hcd_score);  
				fout<<buffer<<","<<type ; 
				fout<<","<<protein_seq<< "," << hcd_record.mdbl_seq_mass<<","<<peptide_seq<<","<<contains_site<<","<<hcd_record.mdbl_glycan_mass<<","; 
				fout<<glycan_comp<< "," <<hcd_record.mdbl_mass_error <<","; 	
				fout<<hcd_record.mdbl_y1_mz<<","<<hcd_record.mdbl_parent_scan_time ; 
				fout<<"\n" ; 
			}
			fout.close() ; 


		}		
		
		
		void HCDScoringResults::LoadResultsHCD(char *inpFile)
		{
			std::ifstream fin(inpFile, std::ios::binary) ; 
			mvect_hcd_records.clear() ; 
			Engine::MS2HCDScoring::HCDInformationRecord hcd_record ; 

			// TO DO
		}

		
	}
}