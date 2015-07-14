#include "clscidscoring.h"
#include "clsPeakProcessor.h"
#include "clshorntransform.h"


namespace GlypID
{
	namespace CIDScoring
	{
		clsCIDScoring::clsCIDScoring(void)
		{
			try
			{
				mobj_cid_scoring = new Engine::MS2CIDScoring::CIDScoring() ;  
				mdbl_cc_mass  = 1.00727638 ;
			}
			catch (char *mesg)
			{
				System::String *exception_msg = new System::String(mesg) ; 
				throw new System::Exception(exception_msg) ; 
			}

		}

		clsCIDScoring::~clsCIDScoring(void)
		{
		}

		void clsCIDScoring::SetOptions(short max_charge, double max_mw, double max_fit, double min_s2n, double cc_mass, 
			double delete_threshold_theoretical_intensity, double min_theoretical_intensity_for_score, short num_peaks_for_shoulder, 
			bool use_caching, bool o16_o18_media, bool check_against_charge_1)
		{
			mobj_transform->SetOptions(max_charge, max_mw, max_fit, min_s2n, cc_mass, delete_threshold_theoretical_intensity,
				min_theoretical_intensity_for_score, num_peaks_for_shoulder, check_against_charge_1, use_caching, o16_o18_media) ; 
		}

		void clsCIDScoring::SetIsotopeFitOptions(System::String *str_averagine, System::String *str_tag, bool thrash_or_not, 
			bool complete_fit)
		{
			char averagine_formula[512] ;
			char tag_formula[512] ; 
			averagine_formula[0] = '\0' ; 
			tag_formula[0] = '\0' ; 

			GlypID::Utils::GetStr(str_averagine, averagine_formula) ; 
			if (str_tag != NULL)
			{
				GlypID::Utils::GetStr(str_tag, tag_formula) ; 
			}

			mobj_transform->SetIsotopeFitOptions(averagine_formula, tag_formula, thrash_or_not, complete_fit) ; 
		}	

		void clsCIDScoring::SetScanDetails(int msNScan, int parent_scan, int msNLevel, int parentLevel, double parentMz )
		{
			mint_msn_scan = msNScan ; 
			mint_parent_scan = parent_scan ; 
			mint_msn_level = msNLevel; 
			mint_parent_level = parentLevel ; 
			mdbl_parent_mz = parentMz ; 
		}

		void clsCIDScoring::SearchForY1Ion(GlypID::Sequence::clsSequence* (&peptides) __gc[],GlypID::Peaks::clsPeak* (&cid_peaks) __gc[], GlypID::CIDScoring::clsCIDScoringScanResults* (&cid_scoring_results) __gc[])
		{
			Engine::PeakProcessing::PeakData peakData ; 
			GlypID::Utils::SetPeaks(peakData, cid_peaks) ; 
			std::vector <Engine::MS2CIDScoring::CIDInformationRecord> cid_records ; 
			std::vector <Engine::SequenceManager::Sequence> vect_nglyco_peps ; 
			std::string temp = " "  ; 
			for(int k=0; k < cid_scoring_results->Length; k++)
			{
				Engine::MS2CIDScoring::CIDInformationRecord this_record ; 
				this_record.AddInfoToCIDRecord(mint_msn_scan, mint_parent_scan, mint_msn_level, mint_parent_level, mdbl_parent_mz, cid_scoring_results[k]->mdbl_mono_mz, 0, cid_scoring_results[k]->mint_mono_intensity, 
					cid_scoring_results[k]->mshort_cs, cid_scoring_results[k]->mdbl_fit, cid_scoring_results[k]->mdbl_mono_mw, cid_scoring_results[k]->mdbl_average_mw,
					cid_scoring_results[k]->mdbl_cid_score, 0,  false, 0 , 0,0) ;  
				
				for (int i=0 ; i < peptides->Length; i++)
				{					
					char seq_string[512] ; 
					char seq_name[512] ; 
					seq_string[0] = '\0' ; 
					seq_name[0] = '\0' ;
					GlypID::Utils::GetStr(peptides[i]->sequence, seq_string); 
					GlypID::Utils::GetStr(peptides[i]->proteinName, seq_name) ; 
					this_record.AddSearchInfoToCIDRecord(seq_name, 0, 0, temp,  seq_string, "",false, 0, false) ;
					
					double y1mz = mobj_cid_scoring->DetermineY1IonThroughPeptideSearchingCID(peakData, this_record); 
					if(y1mz >0)
					{
						// found a y1, could be a candidate
						cid_records.push_back(this_record) ; 						
					}
				}
			}

			cid_scoring_results->Clear(); 
			cid_scoring_results = new GlypID::CIDScoring::clsCIDScoringScanResults* __gc [cid_records.size()] ; 
			for (int i=0; i < cid_records.size() ; i++)
			{
				cid_scoring_results[i] = new clsCIDScoringScanResults() ; 
				cid_scoring_results[i]->Set(cid_records.at(i)) ; 
			}
		}
			
		bool clsCIDScoring::ScoreCIDSpectra(GlypID::Peaks::clsPeak* (&cid_peaks) __gc[], float (&mzs) __gc [], float (&intensities) __gc [], GlypID::HornTransform::clsHornTransformResults* (&precursor_records) __gc[], GlypID::CIDScoring::clsCIDScoringScanResults* (&cid_scoring_results) __gc[]) 
		{
			mint_percent_done = 0 ;
			std::vector<double> vectMzs ;
			std::vector<double> vectIntensities ;
			int numPoints = mzs->Length ; 

			if (mzs->Length == 0)
				return false  ; 

			// mzs should be in sorted order
			double minMZ = mzs[0] ; 
			double maxMZ = mzs[numPoints-1] ; 

			Engine::MS2CIDScoring::CIDInformationRecord cidRecord ; 
			Engine::PeakProcessing::PeakData peakData ; 
			Engine::PeakProcessing::PeakData peakDataForOxonium ; 
			GlypID::Utils::SetPeaks(peakDataForOxonium, cid_peaks) ; 
			
			GlypID::Utils::SetPeaks(peakData, cid_peaks) ; 
			
			double score = 0; 
			double p_value = 1 ; 
			bool found_a_score = false ; 
	
			cid_scoring_results = new GlypID::CIDScoring::clsCIDScoringScanResults* __gc [precursor_records->Length] ; 
			for (int k = 0 ; k < precursor_records->Length; k++)
			{	
				score = mobj_cid_scoring->CalculateCIDScore(peakData, precursor_records[k]->mshort_cs) ; 	
				p_value = mobj_cid_scoring->CalculateCIDScorePValue((int) score) ; 
				
				bool found_oxonium = mobj_cid_scoring->CheckForOxoniumIons(peakData) ; 

				found_a_score = true ; 
				double mono_mz = (precursor_records[k]->mdbl_mono_mw/precursor_records[k]->mshort_cs) + mdbl_cc_mass  ; 
				cidRecord.AddInfoToCIDRecord(mint_msn_scan, mint_parent_scan, mint_msn_level, mint_parent_level, mdbl_parent_mz, mono_mz, 0, precursor_records[k]->mint_abundance, 
					precursor_records[k]->mshort_cs, precursor_records[k]->mdbl_fit, precursor_records[k]->mdbl_mono_mw, precursor_records[k]->mdbl_average_mw,
					score, p_value, found_oxonium, 0 , 0, 0) ;  
				cid_scoring_results[k] = new clsCIDScoringScanResults() ; 
				cid_scoring_results[k]->Set(cidRecord) ; 
			
			}
			

			mint_percent_done = 100 ; 
			return found_a_score; 
		}	


		
	}
}