#include "clshcdscoring.h"
#include "clsPeakProcessor.h"
#include "clshorntransform.h"

namespace GlypID
{
	namespace HCDScoring
	{
		clsHCDScoring::clsHCDScoring(void)
		{
			try
			{				 
				mobj_hcd_scoring = new Engine::MS2HCDScoring::HCDScoring() ; 
			}
			catch (char *mesg)
			{
				System::String *exception_msg = new System::String(mesg) ; 
				throw new System::Exception(exception_msg) ; 
			}

		}

		clsHCDScoring::~clsHCDScoring(void)
		{
		}

		void clsHCDScoring::SetOptions(short max_charge, double max_mw, double max_fit, double min_s2n, double cc_mass, 
			double delete_threshold_theoretical_intensity, double min_theoretical_intensity_for_score, short num_peaks_for_shoulder, 
			bool use_caching, bool o16_o18_media, bool check_against_charge_1)
		{
			mobj_transform->SetOptions(max_charge, max_mw, max_fit, min_s2n, cc_mass, delete_threshold_theoretical_intensity,
				min_theoretical_intensity_for_score, num_peaks_for_shoulder, check_against_charge_1, use_caching, o16_o18_media) ; 
		}

		void clsHCDScoring::SetIsotopeFitOptions(System::String *str_averagine, System::String *str_tag, bool thrash_or_not, 
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

		void clsHCDScoring::SetScanDetails(int msNScan, int parent_scan, int msNLevel, int parentLevel, double parentMz )
		{
			mint_msn_scan = msNScan ; 
			mint_parent_scan = parent_scan ; 
			mint_msn_level = msNLevel; 
			mint_parent_level = parentLevel ; 
			mdbl_parent_mz = parentMz ; 

		}

		void clsHCDScoring::SearchForY1Ion(GlypID::Peaks::clsPeak* (&hcd_peaks) __gc[], GlypID::ETDScoring::clsETDScoringScanResults* (&etd_scoring_results) __gc[],
			GlypID::HornTransform::clsHornTransformResults* (&precursor_records) __gc[]) 
		{
			Engine::PeakProcessing::PeakData peakData ; 
			GlypID::Utils::SetPeaks(peakData, hcd_peaks) ; 			
			
			for (int i = 0 ; i < precursor_records->Length ; i++)
			{
				for(int k=0; k < etd_scoring_results->Length; k++)
				{
					Engine::MS2CIDScoring::CIDInformationRecord this_record ; 
					this_record.AddInfoToCIDRecord(mint_msn_scan, mint_parent_scan, mint_msn_level, mint_parent_level, mdbl_parent_mz, precursor_records[i]->mdbl_mz, 0, precursor_records[i]->mint_abundance, 
						precursor_records[i]->mshort_cs, precursor_records[i]->mdbl_fit, precursor_records[i]->mdbl_mono_mw,  precursor_records[i]->mdbl_average_mw, 0, 0,  false, 0 , 0,0) ;  
					
					char seq_string[512] ; 
					char seq_name[512] ; 
					char glycan[512] ; 
					char site_position[512] ; 
					seq_string[0] = '\0' ; 
					seq_name[0] = '\0' ;
					glycan[0] = '\0' ; 
					site_position[0] = '\0' ; 

					GlypID::Utils::GetStr(etd_scoring_results[k]->mstr_pep_seq, seq_string); 
					GlypID::Utils::GetStr(etd_scoring_results[k]->mstr_pro_seq_name, seq_name) ; 
					GlypID::Utils::GetStr(etd_scoring_results[k]->mstr_glycan_composition, glycan) ; 
					GlypID::Utils::GetStr(etd_scoring_results[k]->mstr_glyco_site, site_position) ; 

					this_record.AddSearchInfoToCIDRecord(seq_name,
						etd_scoring_results[k]->mdbl_seq_mass, 
						etd_scoring_results[k]->mdbl_glycan_mass,
						glycan, seq_string, site_position, true, 0, 0) ;

					double y1mz  = mobj_hcd_scoring->DetermineY1IonThroughPeptideSearching(peakData, this_record) ; 
					etd_scoring_results[k]->mdbl_y1_mz = y1mz ;  
				}
			}
		}
			
		double clsHCDScoring::ScoreHCDSpectra(GlypID::Peaks::clsPeak* (&hcd_peaks) __gc[], float (&mzs) __gc [], float (&intensities) __gc [],  GlypID::HornTransform::clsHornTransformResults* (&precursor_records) __gc[], GlypID::HCDScoring::clsHCDScoringScanResults* (&hcd_scoring_results) __gc[]) 
		{
			mint_percent_done = 0 ;
			std::vector<double> vectMzs ;
			std::vector<double> vectIntensities ;
			std::vector<int> vectPeakIndices ; 
			int numPoints = mzs->Length ; 

			if (mzs->Length == 0)
				return 0  ; 

			// mzs should be in sorted order
			double minMZ = mzs[0] ; 
			double maxMZ = mzs[numPoints-1] ; 

			Engine::MS2HCDScoring::HCDInformationRecord hcdRecord ; 
			Engine::PeakProcessing::PeakData peakData ; 
			Engine::MS2HCDScoring::GLYCAN_TYPE glycanType ; 
			GlypID::Utils::SetPeaks(peakData, hcd_peaks) ; 

			//Options
			mobj_hcd_scoring->SetOptions(mobj_scoring_parameters->get_MinNumPeaksToConsider(), mobj_scoring_parameters->get_MinHCDMz(), mobj_scoring_parameters->get_MaxHCDMz()) ; 
			
			// Score			
			double	score = mobj_hcd_scoring->DetermineGlycanType(peakData, vectPeakIndices, glycanType) ; 
			hcd_scoring_results = new GlypID::HCDScoring::clsHCDScoringScanResults* __gc [precursor_records->Length] ; 
			for (int k = 0 ; k < precursor_records->Length; k++)
			{						
				hcdRecord.AddInfoToHCDRecord(mint_msn_scan, mint_parent_scan, mint_msn_level, mint_parent_level, mdbl_parent_mz, precursor_records[k]->mdbl_mz, 0, precursor_records[k]->mint_abundance, 
					precursor_records[k]->mshort_cs, precursor_records[k]->mdbl_fit, precursor_records[k]->mdbl_mono_mw, score, vectPeakIndices, glycanType) ;  
				
				hcd_scoring_results[k] = new clsHCDScoringScanResults() ; 
				hcd_scoring_results[k]->Set(hcdRecord) ; 
			}			
			

			mint_percent_done = 100 ; 
			return score ; 

		}	
	}
}