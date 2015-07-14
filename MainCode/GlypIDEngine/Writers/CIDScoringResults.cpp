// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
#include "CIDScoringResults.h"
#include <fstream> 
#include <iostream> 
#include <time.h> 

namespace Engine
{
	namespace Writers
	{		

		CIDScoringResults::CIDScoringResults()
		{		
			
			mint_num_records = 0 ; 
		}

		CIDScoringResults::~CIDScoringResults()
		{
		
			ClearCIDRecords() ; 
			
		}

		void CIDScoringResults::ClearCIDRecords()
		{
			
			mvect_cid_records.clear() ; 
			mint_num_records = 0 ; 
		}

		void CIDScoringResults::AddCIDRecord(Engine::MS2CIDScoring::CIDInformationRecord cid_record)
		{
			mvect_cid_records.push_back(cid_record) ; 
		}

		void CIDScoringResults::AddMultipleCIDRecords(std::vector<Engine::MS2CIDScoring::CIDInformationRecord> &vect_cid_records)
		{			
				if (vect_cid_records.size() == 0)
					return ; 

				mvect_cid_records.insert(mvect_cid_records.end(), vect_cid_records.begin(), vect_cid_records.end()) ; 
				mint_num_records += vect_cid_records.size() ; 
		}

	
		
		void CIDScoringResults::SaveResultsCID(char *cidFile) 
		{
			char cid_file_name[512] ; 
			strcpy(cid_file_name, cidFile) ; 
			strcat(cid_file_name, "_cid.csv") ;
			std::ofstream fout(cid_file_name);

			fout<<"MSn_Scan"<<","<<"MSn_Level"<<","<<"Parent_Scan"<<","<<"Parent_Scan_Level"<<","<<"Parent_Mz" ; 
			fout<<","<<"Mono_Mz"<<","<<"Charge_State"<<","<<"Monoisotopic_Mass"<<","<<"Isotopic_Fit"; 
			fout<<","<<"Mono_Intensity"<<","<<"CID_Score"<<","<<"CID_P_Value"<<","<<"Contains_Oxonium"<<","; 
			fout << "Protein_name"<<","<<"Peptide_mass"<<","<<"Sequence"<<","<<"Contains_Site"<<","<<"Glycan_Mass"<<",";
			fout<<"Glycan_Composition"<<"," ; 
			if (mvect_cid_records[0].mbln_mass_error_in_ppm) 
				fout<<"Error(ppm)"<<",";
			else
				fout<<"Error(Da)"<<"," ; 
			fout<<"Y1_ion"<<","<<"Parent_time"<<","<<"CID_time"; 
			fout<<"\n" ; 

			fout.precision(4) ; 
			fout.setf(std::ios::fixed, std::ios::floatfield);

			Engine::MS2CIDScoring::CIDInformationRecord cid_record ; 

			for (int i = 0 ; i < mvect_cid_records.size() ; i++)
			{
				cid_record = mvect_cid_records[i] ; 

				const char *protein_seq ; 
				const char *peptide_seq ;
				const char *glycan_comp ; 
				const char *contains_site ;
				const char *contains_oxonium ; 

				if (cid_record.mbln_contains_oxonium_ion)
					contains_oxonium = "yes" ; 
				else
					contains_oxonium = "no" ; 
				
					

				protein_seq = cid_record.mstr_pro_seq_name.c_str();
				peptide_seq = cid_record.mstr_pep_seq_name.c_str() ; 
				glycan_comp = cid_record.mstr_glycan_composition.c_str() ; 
				contains_site = cid_record.mstr_glyco_site.c_str() ; 

				fout<<cid_record.mint_msn_scan_num<<","<<cid_record.mint_msn_scan_level<<","<<cid_record.mint_parent_scan_num ; 
				fout<<","<<cid_record.mint_parent_scan_level<<","<<cid_record.mdbl_parent_mz ; 
				fout<<","<<cid_record.mdbl_mono_mz<<","<<cid_record.mshort_cs; 
				fout<<","<<cid_record.mdbl_mono_mw<<","<<cid_record.mdbl_fit ; 
				fout<<","<<cid_record.mint_mono_intensity<<","<<cid_record.mdbl_cid_score<<",";
				fout<<cid_record.mdbl_cid_score_p_value<<","<<contains_oxonium; 
				fout<<","<<protein_seq<< "," << cid_record.mdbl_seq_mass<<","<<peptide_seq<<","<<contains_site<<","<<cid_record.mdbl_glycan_mass<<","; 
				fout<<glycan_comp<< "," <<cid_record.mdbl_mass_error <<","; 	
				fout<<cid_record.mdbl_y1_mz<<","<<cid_record.mdbl_parent_scan_time<<","<<cid_record.mdbl_cid_scan_time ; 
				fout<<"\n" ; 
			}
			fout.close() ; 


		}		
		
		
		void CIDScoringResults::LoadResultsCID(char *inpFile)
		{
			std::ifstream fin(inpFile, std::ios::binary) ; 
			mvect_cid_records.clear(); 
			

			// TO DO
		}

		
	}
}