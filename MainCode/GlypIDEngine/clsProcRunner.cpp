 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from Decon2LS code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)


#include<fstream>
#include ".\clsprocrunner.h"
#include "Readers\ReaderFactory.h"
#include "GlypIDEngineUtils.h"
#include "PeakProcessor/PeakProcessor.h"
#include "clsTransformResults.h"
#include "HornTransform/MassTransform.h"
#include "Utilities/SavGolSmoother.h"
#include "MS2HCDScore/HCDScoring.h"
#include "MS2CIDScore/CIDScoring.h"
#include "MS2ETDScore/ETDScoring.h"
#include "GlycanCompositionManager/GlycanProcessor.h"
#include "Writers/HCDScoringResults.h"
#include "Writers/CIDHCDCombinedScoringResults.h"
#include "Writers/CIDETDCombinedScoringResults.h"
#include "MS2HCDSCore/HCDInformationRecord.h"
#include "MS2CIDScore/CIDInformationRecord.h"
#include "MS2CIDHCDCombinedScore/CIDHCDCombinedInformationRecord.h"
//#include "MS2CIDETDCombinedScore/CIDETDCombinedInformationRecord.h"
#include "SequenceManager/Sequence.h"
#include "Readers/FastaFile.h" 

#include <time.h>
#include <map>

namespace GlypID
{
	clsProcRunner::clsProcRunner(void)
	{
		mint_percent_done = 0;
		menm_state = IDLE ;
		mstr_file_name = NULL ;
		mdbl_cc_mass  = 1.00727638 ;			

		mbln_search_fasta = false ; 

		mobj_peak_parameters = new GlypID::Peaks::clsPeakProcessorParameters();
		mobj_transform_parameters = new GlypID::HornTransform::clsHornTransformParameters();		
		mobj_transform_results = NULL ;
		mobj_hcd_scoring_results = NULL ; 
		mobj_cid_etd_combined_scoring_results = NULL ; 
		mobj_cid_hcd_combined_scoring_results = NULL ; 
	}

	clsProcRunner::~clsProcRunner(void)
	{
		if (mobj_transform_results != NULL)
		{
			delete mobj_transform_results ; 
			mobj_transform_results = NULL ; 
		}			
		if (mobj_hcd_scoring_results != NULL)
		{
			delete mobj_hcd_scoring_results ; 
			mobj_hcd_scoring_results = NULL ; 
		}
		if (mobj_cid_hcd_combined_scoring_results != NULL)
		{
			delete mobj_cid_hcd_combined_scoring_results ; 
			mobj_cid_hcd_combined_scoring_results = NULL ; 
		}
		if (mobj_cid_etd_combined_scoring_results != NULL)
		{
			delete mobj_cid_etd_combined_scoring_results ; 
			mobj_cid_etd_combined_scoring_results = NULL ; 
		}
	}

	GlypID::Results::clsTransformResults* clsProcRunner::CreateTransformResults(System::String *file_name, GlypID::Readers::FileType file_type,
		GlypID::Peaks::clsPeakProcessorParameters *peak_parameters,
		GlypID::HornTransform::clsHornTransformParameters *transform_parameters)
	{
		//--------- Function to create Horn Transform Results for entire dataset -------- //
		if (menm_state == RUNNING)
		{
			throw new System::Exception(S"Process already running in clsProcRunner. Cannot run two processes with same object");
		}

		if (file_name == NULL || file_name == S"")
		{
			throw new System::Exception(S"Please enter a file name to process");
		}

		if (peak_parameters == NULL)
		{
			throw new System::Exception(S"Please specify peak processing parameters.");
		}
		if (transform_parameters == NULL)
		{
			throw new System::Exception(S"Please specify mass transform parameters.");
		}
		
		Engine::Readers::RawData __nogc *raw_data = NULL ; 
		Engine::Writers::LCMSTransformResults __nogc *lcms_results = NULL ; 
		Engine::PeakProcessing::PeakProcessor __nogc *peak_processor = NULL ;
		Engine::PeakProcessing::PeakProcessor __nogc *original_peak_processor = NULL ; 
		Engine::HornTransform::MassTransform __nogc *mass_transform = NULL ;
	

		GlypID::Results::clsTransformResults *transform_results = __gc new GlypID::Results::clsTransformResults() ;

		try
		{
			std::vector<Engine::HornTransform::IsotopeFitRecord> vect_transform_records ; 

			// while the thresholded parameter is already set in the clsPeakProcessParameters, we would
			// like to override that here if the data type is Finnigan or mzXMl because that data is threshold.			
			bool thresholded = true; 

			
			mint_percent_done = 0 ;
			menm_state = RUNNING ;

			// Create a RawData object and read through each scan and discover peaks.

			char file_name_ch[256] ;
			GlypID::Utils::GetStr(file_name, file_name_ch) ;
			// enumerations of file type are the same in Readers namespace			
			raw_data = Engine::Readers::ReaderFactory::GetRawData((Engine::Readers::FileType)file_type, file_name_ch) ;
			if (raw_data == NULL)
			{
				throw new System::Exception(System::String::Concat(S"Could not open raw file: ", file_name));
			}

			lcms_results = new Engine::Writers::LCMSTransformResults() ;
			peak_processor = new Engine::PeakProcessing::PeakProcessor() ;
			original_peak_processor = new Engine::PeakProcessing::PeakProcessor() ;
			
			
			// Set parameters for discovering peaks. intensity threshold is set below.
			peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
			
			
			// Set parameters for horn transform (or mass transform)			
			mass_transform = new Engine::HornTransform::MassTransform() ;
			mass_transform->SetElementalIsotopeComposition(*transform_parameters->get_ElementIsotopeComposition()->mobjAtomicInfo) ;
			if (transform_parameters != NULL)
			{
				char averagine_formula[512] ;
				char tag_formula[512] ;
				mass_transform->SetOptions(transform_parameters->get_MaxCharge(), transform_parameters->get_MaxMW(), transform_parameters->get_MaxFit(),
					transform_parameters->get_MinS2N(), transform_parameters->get_CCMass(),transform_parameters->get_DeleteIntensityThreshold(),
					transform_parameters->get_MinIntensityForScore(), transform_parameters->get_NumPeaksForShoulder(),
					transform_parameters->get_CheckAllPatternsAgainstCharge1(), transform_parameters->get_UseMercuryCaching(), transform_parameters->get_O16O18Media()) ;

				averagine_formula[0] = '\0' ;
				tag_formula[0] = '\0' ;

				GlypID::Utils::GetStr(transform_parameters->get_AveragineFormula(), averagine_formula) ;
				if (transform_parameters->get_TagFormula() != NULL)
				{
					GlypID::Utils::GetStr(transform_parameters->get_TagFormula(), tag_formula) ;
				}
				mass_transform->SetIsotopeFitOptions(averagine_formula, tag_formula, transform_parameters->get_ThrashOrNot(), transform_parameters->get_CompleteFit()) ;
				mass_transform->SetIsotopeFitType((Engine::HornTransform::IsotopicFittingType) transform_parameters->get_IsotopeFitType()) ; 
			}	

			std::vector<double> vect_mzs ;
			std::vector<double> vect_intensities;
			std::vector<double> temp_vect_mzs ;
			std::vector<double> temp_vect_intensities;

			// Scan Range
			int min_scan = 1 ; 
			if (transform_parameters->get_UseScanRange() && transform_parameters->get_MinScan() > 1)
				min_scan = transform_parameters->get_MinScan() ;
			if (min_scan < raw_data->GetFirstScanNum())
				min_scan = raw_data->GetFirstScanNum() ; 
			int max_scan = raw_data->GetLastScanNum() ;		
			if (transform_parameters->get_UseScanRange() && transform_parameters->get_MaxScan() < max_scan)
				max_scan = transform_parameters->get_MaxScan() ;
			
			
			int raw_data_read_time = 0 ;
			int preprocessing_time = 0 ;
			int transform_time = 0 ;
			double scan_time = 0 ;
			double drift_time = 0 ;
			short scan_ms_level = 1 ;

			int average_time = 0 ;
			int current_time = 0 ;
			int tic_time = 0 ;
			int peak_discover_time = 0 ;
			int peak_save_time = 0 ;

			// Start processing each scan			
			for(int scan_num = min_scan ; scan_num <= max_scan && scan_num != -1 ; scan_num = raw_data->GetNextScanNum(scan_num))
			{				
				
				peak_processor->Clear() ;
				mint_current_scan = scan_num ;
				if (min_scan != max_scan)
					mint_percent_done = ((scan_num-min_scan)*100)/(max_scan - min_scan)  ;
				vect_mzs.clear() ; 
				vect_intensities.clear() ; 
				temp_vect_mzs.clear() ; 
				temp_vect_intensities.clear() ; 

				//Check if it needs to be deconcoluted
				if (raw_data->IsMSScan(scan_num))
				{
					scan_ms_level = 1 ; 
				}
				else 
				{
					scan_ms_level = 2 ; 
				}
				if (scan_ms_level != 1 && !transform_parameters->get_ProcessMSMS())
					continue ; //forget MSn scans

				//Get this scan first
				raw_data->GetRawData(&vect_mzs, &vect_intensities, scan_num) ;				
				scan_time = raw_data->GetScanTime(scan_num) ;			
				if (vect_mzs.size() == 0)
				{			
					lcms_results->AddInfoForScan(scan_num, 0, 0, 0, 0, 0, 0, scan_time, scan_ms_level) ;				
					continue ;
				}

				// get some characterteristics for that scan				
				double minMZ = vect_mzs[0] ;
				double maxMZ = vect_mzs[(int)vect_mzs.size()-1] ;				
				double thres = GlypID::Utils::GetAverage(vect_intensities, FLT_MAX) ;
				double background_intensity = GlypID::Utils::GetAverage(vect_intensities, (float)(5*thres)) ;			
				double bpi = 0, bp_mz = 0 ;				
				double tic_intensity = 0 ; 
				if (transform_parameters->get_UseMZRange())
				{
					tic_intensity = GlypID::Utils::GetTIC(transform_parameters->get_MinMZ(),
						transform_parameters->get_MaxMZ(), vect_mzs, vect_intensities, 
						(float)(background_intensity * peak_parameters->get_PeakBackgroundRatio()),
						bpi, bp_mz) ;
				}
				else
				{
					tic_intensity = GlypID::Utils::GetTIC(400.0, 2000.0, vect_mzs, vect_intensities, 
						(float)(background_intensity * peak_parameters->get_PeakBackgroundRatio()), bpi, bp_mz) ;
				}
				
				// Set peak processing noise floor				
				peak_processor->SetPeakIntensityThreshold(background_intensity * peak_parameters->get_PeakBackgroundRatio()) ;

				// Find peaks
				int numPeaks = 0 ; 
				if (transform_parameters->get_UseMZRange())
				{					
					numPeaks = peak_processor->DiscoverPeaks(&vect_mzs, &vect_intensities, transform_parameters->get_MinMZ(), transform_parameters->get_MaxMZ()) ;				
				}
				else
				{				
					numPeaks = peak_processor->DiscoverPeaks(&vect_mzs, &vect_intensities) ;				
				}

				// Record them in results				
				lcms_results->AddPeaksForScan(scan_num, peak_processor->mobj_peak_data->mvect_peak_tops) ;				
								

				// Set parameters for isotopic peaks (calculated from the peaks discovered)
				int numDeisotoped = 0 ;
				double min_peptide_intensity = background_intensity * transform_parameters->get_PeptideMinBackgroundRatio() ;
				if (transform_parameters->get_UseAbsolutePeptideIntensity())
				{
					if (min_peptide_intensity < transform_parameters->get_AbsolutePeptideIntensity())
						min_peptide_intensity = transform_parameters->get_AbsolutePeptideIntensity() ;
				}

				Engine::PeakProcessing::Peak currentPeak ;
				Engine::PeakProcessing::Peak originalPeak ; 
				Engine::HornTransform::IsotopeFitRecord transformRecord ;

				// this is important for each scan. It's a clearing up process
				peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;
				
				// Find the most intense peak
				bool found = peak_processor->mobj_peak_data->GetNextPeak(minMZ, maxMZ, currentPeak) ;
				double fwhm_SN = currentPeak.mdbl_FWHM ;

				mass_transform->Reset() ; 
				vect_transform_records.clear() ; 
				while(found)
				{
					if (currentPeak.mdbl_intensity < min_peptide_intensity)
						break ;

					bool found_transform = false ;				
					try
					{
						// Deconvolute
						found_transform = mass_transform->FindTransform(*peak_processor->mobj_peak_data, currentPeak, transformRecord, background_intensity) ;				
						if (found_transform && transformRecord.mshort_cs <= transform_parameters->get_MaxCharge())
						{
							// found an isotopic distribution and have convolved it so record it in results
							numDeisotoped++ ;
							transformRecord.mint_scan_num = scan_num ;
							vect_transform_records.push_back(transformRecord) ;
						}

						// Get next highest peak
						found = peak_processor->mobj_peak_data->GetNextPeak(minMZ, maxMZ, currentPeak) ;
					}
					catch (char *mesg1)
					{
						Console::WriteLine(mesg1) ; 
						Console::WriteLine(Convert::ToString(currentPeak.mint_peak_index)) ; 
					}
				}

				// Add all deconvolved isotopic distributions 
				lcms_results->AddTransforms(vect_transform_records) ; 						
				double signal_range = raw_data->GetSignalRange(scan_num) ; 
				lcms_results->AddInfoForScan(scan_num, bp_mz, bpi, tic_intensity, signal_range, numPeaks, numDeisotoped, scan_time, scan_ms_level) ; 
			}			

			transform_results->SetLCMSTransformResults(lcms_results) ; 
			
			mint_percent_done = 100 ;

		}
		catch (char *mesg)
		{
			if (lcms_results != NULL)
			{
				delete lcms_results ; 
				lcms_results = NULL ;
			}
			
			if (peak_processor != NULL) 
			{
				delete peak_processor ;
				peak_processor = NULL ; 
			}
			if( raw_data != NULL)
			{
				raw_data->Close() ; 
				delete raw_data ; 
			}

			if (mass_transform != NULL)
			{
				delete mass_transform ;
				mass_transform = NULL ; 
			}
			System::String *exception_msg = new System::String(mesg) ; 
			menm_state = enmProcessState::ERROR ;
			throw new System::Exception(exception_msg) ; 
		}		
		if (peak_processor != NULL) 
		{
			delete peak_processor ;
			peak_processor = NULL ; 
		}
		if( raw_data != NULL)
		{
			raw_data->Close() ; 
			delete raw_data ; 
		}

		if (mass_transform != NULL)
		{
			delete mass_transform ;
			mass_transform = NULL ; 
		}
		menm_state = COMPLETE ;
		return transform_results ;
	}	

	GlypID::Results::clsHCDScoringResults* clsProcRunner::CreateHCDScoringResults(System::String *file_name, GlypID::Readers::FileType file_type,
		System::String *protein_file_name, System::String *glycan_composition_file_name, 
		GlypID::Peaks::clsPeakProcessorParameters *peak_parameters,
		GlypID::HornTransform::clsHornTransformParameters *transform_parameters, GlypID::Scoring::clsScoringParameters *scoring_parameters, bool search_fasta)
	{
		// ---------- Function to create HCD scoring results for all HCD spectra ---------- //
		if (menm_state == RUNNING)
		{
			throw new System::Exception(S"Process already running in clsProcRunner. Cannot run two processes with same object");
		}

		if (file_name == NULL || file_name == S"")
		{
			throw new System::Exception(S"Please enter a file name to process");
		}

		if (peak_parameters == NULL)
		{
			throw new System::Exception(S"Please specify peak processing parameters.");
		}
		if (transform_parameters == NULL)
		{
			throw new System::Exception(S"Please specify mass transform parameters.");
		}

		if (search_fasta)
		{
			if (protein_file_name == NULL || protein_file_name == S"")
			{
				throw new System::Exception(S"Please enter a protein list file") ; 
			}
			if (glycan_composition_file_name == NULL || glycan_composition_file_name == S"")
			{
				throw new System::Exception(S"Please enter a glycan composition file") ; 
			}
		}
		
		
		Engine::Readers::RawData __nogc *raw_data = NULL ; 
		Engine::Readers::GlycanIo  __nogc *glycan_data = NULL ; 
		Engine::Writers::HCDScoringResults __nogc *hcd_scoring_results = NULL ; 		
		Engine::PeakProcessing::PeakProcessor __nogc *parent_peak_processor = NULL ;		
		Engine::PeakProcessing::PeakProcessor __nogc *hcd_peak_processor = NULL ;		
		Engine::GlycanCompositionManager::GlycanProcessor __nogc *glycan_processor = NULL ; 
		Engine::HornTransform::MassTransform __nogc *mass_transform = NULL ; 
		Engine::MS2HCDScoring::HCDScoring __nogc *hcd_scoring = NULL ; 
		GlypID::Results::clsHCDScoringResults *results_scoring = __gc new GlypID::Results::clsHCDScoringResults() ; 
		

		try
		{
			std::vector<Engine::SequenceManager::Sequence> vect_candidate_proteins ; 
			std::vector<Engine::GlycanTheoretical::GlycanComposition> vect_glycan_compositions ; 
			std::vector<Engine::MS2HCDScoring::HCDInformationRecord> vect_hcd_score_records ; 
			std::vector<Engine::HornTransform::IsotopeFitRecord> vect_transform_records ; 
			Engine::MS2HCDScoring::HCDInformationRecord this_hcd_record ; 
			Engine::HornTransform::IsotopeFitRecord this_transform_record ; 
		
			// while the thresholded parameter is already set in the clsPeakProcessParameters, we would
			// like to override that here if the data type is Finnigan or mzXMl because that data is threshold.			
			mint_percent_done = 0 ;
			menm_state = RUNNING ;

			// Create a RawData object and read through each scan and discover peaks.
			char file_name_ch[256] ;
			GlypID::Utils::GetStr(file_name, file_name_ch) ;
			// enumerations of file type are the same in Readers namespace			
			raw_data = Engine::Readers::ReaderFactory::GetRawData((Engine::Readers::FileType)file_type, file_name_ch) ;
			if (raw_data == NULL)
			{
				throw new System::Exception(System::String::Concat(S"Could not open raw file: ", file_name));
			}

			if (search_fasta)
			{
				// Now read in fasta file ; 
				char fasta_name_ch[256] ; 
				GlypID::Utils::GetStr(protein_file_name, fasta_name_ch) ;
				bool readfasta ; 
				readfasta = Engine::Readers::TFastaFile::ReadFastaFile(fasta_name_ch, vect_candidate_proteins );  

				// Load up glycan list
				char glycan_name_ch[256] ;
				GlypID::Utils::GetStr(glycan_composition_file_name, glycan_name_ch) ; 
				glycan_data = new Engine::Readers::GlycanIo(); 
				glycan_data->LoadGlycanListFromFile(glycan_name_ch, vect_glycan_compositions) ; 
			}
			
			hcd_scoring_results = new Engine::Writers::HCDScoringResults() ; 	
			parent_peak_processor = new Engine::PeakProcessing::PeakProcessor() ;
			hcd_peak_processor = new Engine::PeakProcessing::PeakProcessor() ; 
			glycan_processor = new Engine::GlycanCompositionManager::GlycanProcessor() ; 
			mass_transform = new Engine::HornTransform::MassTransform() ;
			hcd_scoring = new Engine::MS2HCDScoring::HCDScoring() ; 
						
			// Set parameters for discovering peaks. intensity threshold is set below.
			bool thresholded = true ; 
			parent_peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
			hcd_peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
					
			// Set parameters for horn transform (or mass transform) of the precursor	
			mass_transform->SetElementalIsotopeComposition(*transform_parameters->get_ElementIsotopeComposition()->mobjAtomicInfo) ;
			if (transform_parameters != NULL)
			{
				char averagine_formula[512] ;
				char tag_formula[512] ;
				mass_transform->SetOptions(transform_parameters->get_MaxCharge(), transform_parameters->get_MaxMW(), transform_parameters->get_MaxFit(),
					transform_parameters->get_MinS2N(), transform_parameters->get_CCMass(),transform_parameters->get_DeleteIntensityThreshold(),
					transform_parameters->get_MinIntensityForScore(), transform_parameters->get_NumPeaksForShoulder(),
					transform_parameters->get_CheckAllPatternsAgainstCharge1(), transform_parameters->get_UseMercuryCaching(), transform_parameters->get_O16O18Media()) ;

				averagine_formula[0] = '\0' ;
				tag_formula[0] = '\0' ;

				GlypID::Utils::GetStr(transform_parameters->get_AveragineFormula(), averagine_formula) ;
				if (transform_parameters->get_TagFormula() != NULL)
				{
					GlypID::Utils::GetStr(transform_parameters->get_TagFormula(), tag_formula) ;
				}
				mass_transform->SetIsotopeFitOptions(averagine_formula, tag_formula, transform_parameters->get_ThrashOrNot(), transform_parameters->get_CompleteFit()) ;
				mass_transform->SetIsotopeFitType((Engine::HornTransform::IsotopicFittingType) transform_parameters->get_IsotopeFitType()) ; 
			}

			// Set parameters for scoring
			hcd_scoring->SetOptions(scoring_parameters->get_MinNumPeaksToConsider(), scoring_parameters->get_MinHCDMz(), scoring_parameters->get_MaxHCDMz()) ; 

	
			// Scan Range
			int min_scan = 1 ; 
			if (transform_parameters->get_UseScanRange() && transform_parameters->get_MinScan() >1)
				min_scan = transform_parameters->get_MinScan() ; 
			if (min_scan < raw_data->GetFirstScanNum())
				min_scan = raw_data->GetFirstScanNum() ; 			
			int max_scan = raw_data->GetLastScanNum() ; 
			if (transform_parameters->get_UseScanRange() && transform_parameters->get_MaxScan() < max_scan)
				max_scan = transform_parameters->get_MaxScan() ; 
			short scan_ms_level = 1 ;			
			int scan_num = min_scan ; 
			
			// Prelim peptide processing
			std::vector<double> peptide_masses ; 	
			if (search_fasta)
			{
				if(scoring_parameters->get_UsePPM())
					glycan_processor->SetPPMError(scoring_parameters->get_PPMTolerance()) ; 
				else
					glycan_processor->SetDAError(scoring_parameters->get_DaTolerance()) ; 
				glycan_processor->GetSequenceMassesFromProteins(vect_candidate_proteins, peptide_masses) ; 				
			}

			// precursor details
			int parent_scan = 0  ; 
			double parent_mz = -1.0 ; 
			std::vector<double> parent_vect_mzs ;
			std::vector<double> parent_vect_intensities;
			
			// hcd 
			std::vector<double> hcd_vect_mzs ; 
			std::vector<double> hcd_vect_intensities ; 
			
			//start process			 			
			while(scan_num <= max_scan)
			{					
				mint_percent_done =  (int)((scan_num - min_scan)/(max_scan - min_scan))  * 100 ; 
				//check if scan is HCD				
				if (raw_data->IsHCDScan(scan_num))					
				{					
					// get parent_scan and mz
					parent_scan = raw_data->GetParentScan(scan_num) ; 
					parent_mz = raw_data->GetParentMz(scan_num) ; 

					// get all mzs and intensities
					parent_vect_intensities.clear();
					parent_vect_mzs.clear() ; 
					raw_data->GetRawData(&parent_vect_mzs, &parent_vect_intensities, parent_scan);			
		
					// set noise floor
					double thres =  GlypID::Utils::GetAverage(parent_vect_intensities, FLT_MAX) ; 
					double background_intensity = GlypID::Utils::GetAverage(parent_vect_intensities, (float)(5*thres)) ;
					parent_peak_processor->SetPeakIntensityThreshold(background_intensity * peak_parameters->get_PeakBackgroundRatio()) ; 
					parent_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(parent_scan)) ;
			
					// discover peaks
					int numPeaks = parent_peak_processor->DiscoverPeaks(&parent_vect_mzs, &parent_vect_intensities) ; 
					parent_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;											

					//// Invariant : discovering the parent peak has to be done before transforming the precursor distribution.  Since, THRASH removes the distribution from the spectra, 
					// trying to discover the parent peak after will result in error/
					Engine::PeakProcessing::Peak parentPeak ; 
					double parent_match = parent_peak_processor->GetClosestPeakMz(parent_mz, parentPeak) ; 
					

					// mass transformm
					// Set  noise floor parameters for isotopic peaks (calculated from the peaks discovered)
					int numDeisotoped = 0 ;
					double min_peptide_intensity = background_intensity * transform_parameters->get_PeptideMinBackgroundRatio() ;
					if (transform_parameters->get_UseAbsolutePeptideIntensity())
					{
						if (min_peptide_intensity < transform_parameters->get_AbsolutePeptideIntensity())
								min_peptide_intensity = transform_parameters->get_AbsolutePeptideIntensity() ;
					}

					mass_transform->Reset() ; 
					vect_transform_records.clear() ; 
					bool found_transform = mass_transform->FindPrecursorTransform(*parent_peak_processor->mobj_peak_data, parent_mz, &vect_transform_records, background_intensity, min_peptide_intensity) ;				

					// Now process HCD spectra
					hcd_vect_intensities.clear();
					hcd_vect_mzs.clear() ; 
					raw_data->GetRawData(&hcd_vect_mzs, &hcd_vect_intensities, scan_num);			
		
					// set noise floor
					thres =  GlypID::Utils::GetAverage(hcd_vect_intensities, FLT_MAX) ; 
					background_intensity = 0 ; //GlypID::Utils::GetAverage(hcd_vect_intensities, (float)(5*thres)) ;
					hcd_peak_processor->SetPeakIntensityThreshold(background_intensity * peak_parameters->get_PeakBackgroundRatio()) ; 
					hcd_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(scan_num)) ;
					hcd_peak_processor->mobj_peak_data->Clear() ; 
			
					// discover peaks
					int numHCDPeaks = hcd_peak_processor->DiscoverPeaks(&hcd_vect_mzs, &hcd_vect_intensities, scoring_parameters->get_MinHCDMz(), scoring_parameters->get_MaxHCDMz()) ; 
					hcd_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;											

					if (!found_transform) 
					{
						// Look in header		
						short cs = raw_data->GetMonoChargeFromHeader(scan_num) ; 
						if (cs > 0)
						{
							//create dummy transform record with charge cs
							this_transform_record.mshort_cs = cs ; 
							this_transform_record.mdbl_mono_mw = (parent_mz - mdbl_cc_mass) * cs ; 
							this_transform_record.mdbl_fit = 1; 
							this_transform_record.mint_abundance = (int) parentPeak.mdbl_intensity ; 
							vect_transform_records.push_back(this_transform_record) ; 
						}
						else
						{
							// instrument couldnt find charge as well, so create dummy transforms with cs 2 and 3 based on user preferences
							if (scoring_parameters->get_AllocateDefaultCharges())
							{
								for( int charge = 2 ; charge <= 3 ; charge++)
								{
									this_transform_record.mshort_cs = charge ; 
									this_transform_record.mdbl_mono_mw = (parent_mz - mdbl_cc_mass) * charge ; 
									this_transform_record.mdbl_fit = 1; 
									if (parent_match == 0) // couldn't find parent peak so set intensity to background
										this_transform_record.mint_abundance = (int) min_peptide_intensity  ; 
									else
										this_transform_record.mint_abundance = (int) parentPeak.mdbl_intensity  ; 
									vect_transform_records.push_back(this_transform_record) ; 
								}
							}
						}
					}
					
					// get score
					std::vector <int> vectPeakIndices ; 
					Engine::MS2HCDScoring::GLYCAN_TYPE glycanType ; 					
					double hcd_score_value =hcd_scoring->DetermineGlycanType(*hcd_peak_processor->mobj_peak_data, vectPeakIndices, glycanType) ; 
					
					for (int k = 0; k < (int) vect_transform_records.size(); k++)
					{
						this_transform_record = vect_transform_records[k] ; 
						double mono_mz = (this_transform_record.mdbl_mono_mw/this_transform_record.mshort_cs) + mdbl_cc_mass  ; 
			
						double y1mz = 0.0 ; 
						double parent_elution_time = 0.0 ; 
						this_hcd_record.Clear() ; 
						this_hcd_record.AddInfoToHCDRecord(scan_num, parent_scan, scan_ms_level, 1, parent_mz, mono_mz, 0, this_transform_record.mint_abundance, 
							this_transform_record.mshort_cs, this_transform_record.mdbl_fit, this_transform_record.mdbl_mono_mw, hcd_score_value, vectPeakIndices, glycanType) ; 

						if (search_fasta)
						{							
							glycan_processor->SearchAndFillPeptideInformationSingleRecord(vect_glycan_compositions, vect_candidate_proteins, this_hcd_record, glycanType, scoring_parameters->get_FilterCombList()) ; 
							// TBD Anoop y1mz = cid_scoring->DetermineY1IonThroughPeptideSearchingCID(*cid_peak_processor->mobj_peak_data, this_cid_record) ; 											
							// TBD Anoop if (y1mz > 0)
								//parent_elution_time = raw_data->GetScanTime(this_cid_record.mint_parent_scan_num) ; 									
						}	

						vect_hcd_score_records.push_back(this_hcd_record) ; 
						//hcd_scoring_results->AddHCDRecord(this_hcd_record) ; 
					}	
				}
				scan_num ++ ; 
				Console::WriteLine(scan_num) ; 
			}
			hcd_scoring_results->AddMultipleHCDRecords(vect_hcd_score_records) ; 
			results_scoring->SetHCDScoringResults(hcd_scoring_results) ; 
			mint_percent_done = 100 ;
		}
		catch (char *mesg)
		{
			if (hcd_scoring_results != NULL)
			{
				delete hcd_scoring_results ; 
				hcd_scoring_results = NULL ;
			}			
			if (parent_peak_processor != NULL) 
			{
				delete parent_peak_processor ;
				parent_peak_processor = NULL ; 
			}
			if (hcd_peak_processor != NULL) 
			{
				delete hcd_peak_processor ;
				hcd_peak_processor = NULL ; 
			}
			if( raw_data != NULL)
			{
				raw_data->Close() ; 
				delete raw_data ; 
			}
			if (mass_transform != NULL)
			{
				delete mass_transform ;
				mass_transform = NULL ; 
			}
			if(hcd_scoring != NULL)
			{
				delete hcd_scoring ; 
				hcd_scoring = NULL ; 
			}
			System::String *exception_msg = new System::String(mesg) ; 
			menm_state = enmProcessState::ERROR ;
			throw new System::Exception(exception_msg) ; 
		}		
		if (parent_peak_processor != NULL) 
		{
			delete parent_peak_processor ;
			parent_peak_processor = NULL ; 
		}
		if (hcd_peak_processor != NULL)
		{
			delete hcd_peak_processor ; 
			hcd_peak_processor = NULL ; 
		}
		if( raw_data != NULL)
		{
			raw_data->Close() ; 
			delete raw_data ; 
		}
		if (mass_transform != NULL)
		{
			delete mass_transform ;
			mass_transform = NULL ; 
		}
		if (hcd_scoring != NULL)
		{
			delete hcd_scoring ; 
			hcd_scoring = NULL ; 
		}
		menm_state = COMPLETE ;
		return results_scoring ;	
	}
	
GlypID::Results::clsCIDScoringResults* clsProcRunner::CreateCIDScoringResults(System::String *file_name, GlypID::Readers::FileType file_type,
		System::String *protein_file_name, System::String *glycan_composition_file_name , 
		GlypID::Peaks::clsPeakProcessorParameters *peak_parameters,
		GlypID::HornTransform::clsHornTransformParameters *transform_parameters, GlypID::Scoring::clsScoringParameters *scoring_parameters, bool search_fasta)
	{

		// ---------- Function to create CID scoring results for all CID spectra ---------- //
		if (menm_state == RUNNING)
		{
			throw new System::Exception(S"Process already running in clsProcRunner. Cannot run two processes with same object");
		}

		if (file_name == NULL || file_name == S"")
		{
			throw new System::Exception(S"Please enter a file name to process");
		}

		if (peak_parameters == NULL)
		{
			throw new System::Exception(S"Please specify peak processing parameters.");
		}
		if (transform_parameters == NULL)
		{
			throw new System::Exception(S"Please specify mass transform parameters.");
		}

		if (search_fasta)
		{
			if (protein_file_name == NULL || protein_file_name == S"")
			{
				throw new System::Exception(S"Please enter a protein list file") ; 
			}
			if (glycan_composition_file_name == NULL || glycan_composition_file_name == S"")
			{
				throw new System::Exception(S"Please enter a glycan composition file") ; 
			}
		}
		
		Engine::Readers::RawData __nogc *raw_data = NULL ; 
		Engine::Readers::GlycanIo  __nogc *glycan_data = NULL ; 
		Engine::Writers::CIDScoringResults __nogc *cid_scoring_results = NULL ; 		
		Engine::PeakProcessing::PeakProcessor __nogc *parent_peak_processor = NULL ;		
		Engine::PeakProcessing::PeakProcessor __nogc *cid_peak_processor = NULL ;		
		Engine::GlycanCompositionManager::GlycanProcessor __nogc *glycan_processor = NULL ; 
		Engine::HornTransform::MassTransform __nogc *mass_transform = NULL ; 
		Engine::MS2CIDScoring::CIDScoring __nogc *cid_scoring = NULL ; 			
		GlypID::Results::clsCIDScoringResults *results_scoring = __gc new GlypID::Results::clsCIDScoringResults() ; 

		try
		{
			std::vector<Engine::SequenceManager::Sequence> vect_candidate_proteins ; 
			std::vector<Engine::GlycanTheoretical::GlycanComposition> vect_glycan_compositions ; 
			std::vector<Engine::SequenceManager::Sequence> vect_nglyco_peptides ; 

			std::vector<Engine::MS2CIDScoring::CIDInformationRecord> vect_cid_score_records ;			
			std::vector<Engine::HornTransform::IsotopeFitRecord> vect_transform_records ; 
			Engine::MS2CIDScoring::CIDInformationRecord this_cid_record ; 
			Engine::HornTransform::IsotopeFitRecord this_transform_record ; 
		
			// while the thresholded parameter is already set in the clsPeakProcessParameters, we would
			// like to override that here if the data type is Finnigan o                            
			mint_percent_done = 0 ;
			menm_state = RUNNING ;

			// Create a RawData object and read through each scan and discover peaks.
			char file_name_ch[256] ;
			GlypID::Utils::GetStr(file_name, file_name_ch) ;
			// enumerations of file type are the same in Readers namespace			
			raw_data = Engine::Readers::ReaderFactory::GetRawData((Engine::Readers::FileType)file_type, file_name_ch) ;
			if (raw_data == NULL)
			{
				throw new System::Exception(System::String::Concat(S"Could not open raw file: ", file_name));
			}
			
			if (search_fasta)
			{
				// Now read in fasta file ; 
				char fasta_name_ch[256] ; 
				GlypID::Utils::GetStr(protein_file_name, fasta_name_ch) ;
				bool readfasta ; 
				readfasta = Engine::Readers::TFastaFile::ReadFastaFile(fasta_name_ch, vect_candidate_proteins );  

				// Load up glycan list
				char glycan_name_ch[256] ;
				GlypID::Utils::GetStr(glycan_composition_file_name, glycan_name_ch) ; 
				glycan_data = new Engine::Readers::GlycanIo(); 
				glycan_data->LoadGlycanListFromFile(glycan_name_ch, vect_glycan_compositions) ; 
			}

			cid_scoring_results = new Engine::Writers::CIDScoringResults() ; 	
			parent_peak_processor = new Engine::PeakProcessing::PeakProcessor() ;
			cid_peak_processor = new Engine::PeakProcessing::PeakProcessor() ; 
			glycan_processor = new Engine::GlycanCompositionManager::GlycanProcessor() ; 
			mass_transform = new Engine::HornTransform::MassTransform() ;
			cid_scoring = new Engine::MS2CIDScoring::CIDScoring() ; 
						
			// Set parameters for discovering peaks. intensity threshold is set below.
			bool thresholded = true ; 
			parent_peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
			cid_peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
					
			// Set parameters for horn transform (or mass transform) of the precursor	
			mass_transform->SetElementalIsotopeComposition(*transform_parameters->get_ElementIsotopeComposition()->mobjAtomicInfo) ;
			if (transform_parameters != NULL)
			{
				char averagine_formula[512] ;
				char tag_formula[512] ;
				mass_transform->SetOptions(transform_parameters->get_MaxCharge(), transform_parameters->get_MaxMW(), transform_parameters->get_MaxFit(),
					transform_parameters->get_MinS2N(), transform_parameters->get_CCMass(),transform_parameters->get_DeleteIntensityThreshold(),
					transform_parameters->get_MinIntensityForScore(), transform_parameters->get_NumPeaksForShoulder(),
					transform_parameters->get_CheckAllPatternsAgainstCharge1(), transform_parameters->get_UseMercuryCaching(), transform_parameters->get_O16O18Media()) ;

				averagine_formula[0] = '\0' ;
				tag_formula[0] = '\0' ;

				GlypID::Utils::GetStr(transform_parameters->get_AveragineFormula(), averagine_formula) ;
				if (transform_parameters->get_TagFormula() != NULL)
				{
					GlypID::Utils::GetStr(transform_parameters->get_TagFormula(), tag_formula) ;
				}
				mass_transform->SetIsotopeFitOptions(averagine_formula, tag_formula, transform_parameters->get_ThrashOrNot(), transform_parameters->get_CompleteFit()) ;
				mass_transform->SetIsotopeFitType((Engine::HornTransform::IsotopicFittingType) transform_parameters->get_IsotopeFitType()) ; 
			}
				
			// Scan Range
			int min_scan = 1 ; 			
			if (transform_parameters->get_UseScanRange() && transform_parameters->get_MinScan() >1)
				min_scan = transform_parameters->get_MinScan() ; 
			if (min_scan < raw_data->GetFirstScanNum())
				min_scan = raw_data->GetFirstScanNum() ; 			
			int max_scan = raw_data->GetLastScanNum() ; 
			if (transform_parameters->get_UseScanRange() && transform_parameters->get_MaxScan() < max_scan)
				max_scan = transform_parameters->get_MaxScan() ; 		
			short scan_ms_level = 1 ;			
			int scan_num = min_scan ; 
			cid_scoring->InitPTable() ; 


			// Prelim peptide processing
			std::vector<double> peptide_masses ; 	
			if (search_fasta)
			{
				glycan_processor->GetNGlycopeptideSeqList(vect_nglyco_peptides, vect_candidate_proteins);
				vect_candidate_proteins.clear() ; 

				if(scoring_parameters->get_UsePPM())
					glycan_processor->SetPPMError(scoring_parameters->get_PPMTolerance()) ; 
				else
					glycan_processor->SetDAError(scoring_parameters->get_DaTolerance()) ; 
				//glycan_processor->GetSequenceMassesFromProteins(vect_candidate_proteins, peptide_masses) ; 				
				
			}


			// precursor details
			int parent_scan = 0  ; 
			double parent_mz = -1.0 ; 
			std::vector<double> parent_vect_mzs ;
			std::vector<double> parent_vect_intensities;
			
			// cid
			std::vector<double> cid_vect_mzs ; 
			std::vector<double> cid_vect_intensities ; 
			
			//start process				 
			while(scan_num <= max_scan)
			{	
				mint_percent_done = (int) ((scan_num-min_scan)/(max_scan - min_scan)) * 100 ; 
				//check if scan is CID				
				if (raw_data->IsCIDScan(scan_num))
				{						 
					// get parent_scan and mz
					parent_scan = raw_data->GetParentScan(scan_num) ; 
					parent_mz = raw_data->GetParentMz(scan_num) ; 

					// get all mzs and intensities
					parent_vect_intensities.clear();
					parent_vect_mzs.clear() ; 
					raw_data->GetRawData(&parent_vect_mzs, &parent_vect_intensities, parent_scan);			
		
					// set noise floor
					double thres =  GlypID::Utils::GetAverage(parent_vect_intensities, FLT_MAX) ; 
					double background_intensity = GlypID::Utils::GetAverage(parent_vect_intensities, (float)(5*thres)) ;
					parent_peak_processor->SetPeakIntensityThreshold(background_intensity * peak_parameters->get_PeakBackgroundRatio()) ; 
					parent_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(parent_scan)) ;
			
					// discover peaks
					int numPeaks = parent_peak_processor->DiscoverPeaks(&parent_vect_mzs, &parent_vect_intensities) ; 
					parent_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;

					// Invariant : discovering the parent peak has to be done before transforming the precursor distribution.  Since, THRASH removes the distribution from the spectra, 
					// trying to discover the parent peak after will result in error/
					Engine::PeakProcessing::Peak parentPeak ; 
					double parent_match = parent_peak_processor->GetClosestPeakMz(parent_mz, parentPeak) ; 

					// Invariant: Takeing the header information as backup, once again needs to be before THRASH
					Engine::PeakProcessing::Peak monoPeak ; 
					double mono_mz = raw_data->GetMonoMZFromHeader(scan_num) ; 
					double mono_mz_match = parent_peak_processor->GetClosestPeakMz(mono_mz, monoPeak) ; 
				

					// mass transformm
					// Set  noise floor parameters for isotopic peaks (calculated from the peaks discovered)
					int numDeisotoped = 0 ;
					double min_peptide_intensity = background_intensity * transform_parameters->get_PeptideMinBackgroundRatio() ;
					if (transform_parameters->get_UseAbsolutePeptideIntensity())
					{
						if (min_peptide_intensity < transform_parameters->get_AbsolutePeptideIntensity())
								min_peptide_intensity = transform_parameters->get_AbsolutePeptideIntensity() ;
					}

					mass_transform->Reset() ; 
					vect_transform_records.clear() ; 
					bool found_transform = mass_transform->FindPrecursorTransform(*parent_peak_processor->mobj_peak_data, parent_mz, &vect_transform_records, background_intensity, min_peptide_intensity) ;				

					// Now process CID spectra
					cid_vect_intensities.clear();
					cid_vect_mzs.clear() ; 
					raw_data->GetRawData(&cid_vect_mzs, &cid_vect_intensities, scan_num);			
		
					// set noise floor
					thres =  GlypID::Utils::GetAverage(cid_vect_intensities, FLT_MAX) ; 
					background_intensity = 0 ; //GlypID::Utils::GetAverage(hcd_vect_intensities, (float)(5*thres)) ;
					cid_peak_processor->SetPeakIntensityThreshold(background_intensity *peak_parameters->get_PeakBackgroundRatio()) ; 
					cid_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(scan_num)) ;
					cid_peak_processor->mobj_peak_data->Clear() ; 
			
					// discover peaks
					int numCIDPeaks = cid_peak_processor->DiscoverPeaks(&cid_vect_mzs, &cid_vect_intensities) ; 
					cid_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;

					if (!found_transform)
					{
						// Look in header		
						short cs = raw_data->GetMonoChargeFromHeader(scan_num) ; 
						if (cs > 0)
						{
							//create dummy transform record with charge cs
							this_transform_record.mshort_cs = cs ; 
							this_transform_record.mdbl_mono_mw = (mono_mz - mdbl_cc_mass) * cs ; 
							this_transform_record.mdbl_fit = 1; 
							if (mono_mz_match == 0) //couldn't find parent peak so set intensity to background
								this_transform_record.mint_abundance = (int) min_peptide_intensity ; 
							else
								this_transform_record.mint_abundance = (int) monoPeak.mdbl_intensity ; 
							vect_transform_records.push_back(this_transform_record) ; 
						}
						else
						{
							// instrument couldnt find charge as well, so create dummy transforms with cs 2 and 3 based on preference
							if (scoring_parameters->get_AllocateDefaultCharges())
							{
								for( int charge = 2 ; charge <= 3 ; charge++)
								{
									this_transform_record.mshort_cs = charge ; 
									this_transform_record.mdbl_mono_mw = (parent_mz - mdbl_cc_mass) * charge ; 
									this_transform_record.mdbl_fit = 1; 
									if (parent_match == 0) // couldn't find parent peak so set intensity to background
										this_transform_record.mint_abundance = (int) min_peptide_intensity  ; 
									else
										this_transform_record.mint_abundance = (int) parentPeak.mdbl_intensity  ; 
									vect_transform_records.push_back(this_transform_record) ; 
								}
							}
						}
					}
					
					// get score
					Engine::MS2HCDScoring::GLYCAN_TYPE glycanType = Engine::MS2HCDScoring::GLYCAN_TYPE::NA ; // Put a dummy type so that filtering of glycan compositions does not take place
					for (int k = 0; k < (int)vect_transform_records.size(); k++)
					{
						this_transform_record = vect_transform_records[k] ; 
						double cid_score_value = cid_scoring->CalculateCIDScore(*cid_peak_processor->mobj_peak_data, this_transform_record.mshort_cs ) ; 
						double p_value = cid_scoring->CalculateCIDScorePValue((int) cid_score_value) ; 
						bool contains_oxonium = cid_scoring->CheckForOxoniumIons(*cid_peak_processor->mobj_peak_data) ; 
						double mono_mz = (this_transform_record.mdbl_mono_mw/this_transform_record.mshort_cs) + mdbl_cc_mass  ; 
						double y1mz = 0.0 ; 
						double elution_time = raw_data->GetScanTime(scan_num) ; 
						double parent_time = raw_data->GetScanTime(this_cid_record.mint_parent_scan_num) ; 
						this_cid_record.Clear() ; 
						this_cid_record.AddInfoToCIDRecord(scan_num, parent_scan, scan_ms_level, 1, parent_mz, mono_mz, 0, this_transform_record.mint_abundance, 
							this_transform_record.mshort_cs, this_transform_record.mdbl_fit, this_transform_record.mdbl_mono_mw, this_transform_record.mdbl_average_mw, cid_score_value, p_value, contains_oxonium, y1mz, parent_time, elution_time) ; 


						if (search_fasta && cid_score_value >0)
						{
							if (scoring_parameters->get_IncludeBestHit())
							{
								glycan_processor->SearchAndFillPeptideInformationSingleRecord(vect_glycan_compositions, vect_nglyco_peptides, this_cid_record, glycanType,  scoring_parameters->get_FilterCombList()) ; 
								y1mz = cid_scoring->DetermineY1IonThroughPeptideSearchingCID(*cid_peak_processor->mobj_peak_data, this_cid_record) ; 											
								vect_cid_score_records.push_back(this_cid_record) ; 
							}
							else
							{
								std::vector<Engine::MS2CIDScoring::CIDInformationRecord> temp_cid_records ; 
								temp_cid_records.push_back(this_cid_record) ; 
								glycan_processor->SearchAndFillPeptideInformation(vect_glycan_compositions, vect_candidate_proteins, temp_cid_records, glycanType,false, false) ;

								//store into cid records
								vect_cid_score_records.insert(vect_cid_score_records.end(), temp_cid_records.begin(), temp_cid_records.end()); 
								temp_cid_records.clear() ; 
							}							
						}	
						else
						{
							vect_cid_score_records.push_back(this_cid_record) ; 
						}
						
					}			
				}
				scan_num ++ ; 
				Console::WriteLine(scan_num) ; 
			}
			//--------- Peptide Searching ---------//
			/*if (search_fasta)
			{
					glycan_processor->SearchAndFillPeptideInformation(vect_glycan_compositions, vect_candidate_proteins, vect_cid_score_records) ; 
			}*/
			
			//---------- Store and return --------- //
			cid_scoring_results->AddMultipleCIDRecords(vect_cid_score_records) ; 			
			results_scoring->SetCIDScoringResults(cid_scoring_results) ; 
			mint_percent_done = 100 ;
		}
		catch (char *mesg)
		{
			if (cid_scoring_results != NULL)
			{
				delete cid_scoring_results ; 
				cid_scoring_results = NULL ;
			}			
			if (parent_peak_processor != NULL) 
			{
				delete parent_peak_processor ;
				parent_peak_processor = NULL ; 
			}
			if (cid_peak_processor != NULL) 
			{
				delete cid_peak_processor ;
				cid_peak_processor = NULL ; 
			}
			if( raw_data != NULL)
			{
				raw_data->Close() ; 
				delete raw_data ; 
			}
			if (mass_transform != NULL)
			{
				delete mass_transform ;
				mass_transform = NULL ; 
			}
			if(cid_scoring != NULL)
			{
				delete cid_scoring ; 
				cid_scoring = NULL ; 
			}
			System::String *exception_msg = new System::String(mesg) ; 
			menm_state = enmProcessState::ERROR ;
			throw new System::Exception(exception_msg) ; 
		}		
		if (parent_peak_processor != NULL) 
		{
			delete parent_peak_processor ;
			parent_peak_processor = NULL ; 
		}
		if (cid_peak_processor != NULL)
		{
			delete cid_peak_processor ; 
			cid_peak_processor = NULL ; 
		}
		if( raw_data != NULL)
		{
			raw_data->Close() ; 
			delete raw_data ; 
		}
		if (mass_transform != NULL)
		{
			delete mass_transform ;
			mass_transform = NULL ; 
		}
		if (cid_scoring != NULL)
		{
			delete cid_scoring ; 
			cid_scoring = NULL ; 
		}
		menm_state = COMPLETE ;
		return results_scoring ;	
	}


	
	GlypID::Results::clsCIDHCDCombinedScoringResults* clsProcRunner::CreateCIDHCDScoringResults(System::String *file_name, GlypID::Readers::FileType file_type, 
		System::String *protein_file_name, System::String *glycan_composition_file_name , 
		GlypID::Peaks::clsPeakProcessorParameters *peak_parameters,
		GlypID::HornTransform::clsHornTransformParameters *transform_parameters, GlypID::Scoring::clsScoringParameters *scoring_parameters, bool search_fasta)
	{		
		if (menm_state == RUNNING)
		{
			throw new System::Exception(S"Process already running in clsProcRunner. Cannot run two processes with same object");
		}

		if (file_name == NULL || file_name == S"")
		{
			throw new System::Exception(S"Please enter a file name to process");
		}
		if (peak_parameters == NULL)
		{
			throw new System::Exception(S"Please specify peak processing parameters.");
		}
		if (transform_parameters == NULL)
		{
			throw new System::Exception(S"Please specify mass transform parameters.");
		}

		if (search_fasta)
		{
			if (protein_file_name == NULL || protein_file_name == S"")
			{
				throw new System::Exception(S"Please enter a protein list file") ; 
			}
			if (glycan_composition_file_name == NULL || glycan_composition_file_name == S"")
			{
				throw new System::Exception(S"Please enter a glycan composition file") ; 
			}
		}

		// Declarations
		Engine::Readers::RawData __nogc *raw_data = NULL ; 		
		Engine::Readers::GlycanIo  __nogc *glycan_data = NULL ; 
		Engine::PeakProcessing::PeakProcessor __nogc *parent_peak_processor = NULL ;	
		Engine::PeakProcessing::PeakProcessor __nogc *temp_parent_peak_processor = NULL ;	
		Engine::PeakProcessing::PeakProcessor __nogc *cid_peak_processor = NULL ;				
		Engine::PeakProcessing::PeakProcessor __nogc *hcd_peak_processor = NULL ;				
		Engine::HornTransform::MassTransform __nogc *mass_transform = NULL ; 
		Engine::MS2CIDScoring::CIDScoring __nogc *cid_scoring = NULL ; 				
		Engine::MS2HCDScoring::HCDScoring __nogc *hcd_scoring = NULL ; 
		Engine::GlycanCompositionManager::GlycanProcessor __nogc *glycan_processor = NULL ; 
		Engine::Writers::CIDHCDCombinedScoringResults __nogc *comb_scoring_results = NULL ; 
		Results::clsCIDHCDCombinedScoringResults *all_results = __gc new GlypID::Results::clsCIDHCDCombinedScoringResults() ; 

		try
		{
			std::vector<Engine::SequenceManager::Sequence> vect_candidate_proteins ; 
			std::vector<Engine::SequenceManager::Sequence> vect_nglyco_peptides ; 
			std::vector<Engine::GlycanTheoretical::GlycanComposition> vect_glycan_compositions ; 

			std::vector<Engine::MS2CIDScoring::CIDInformationRecord> vect_cid_score_records ; 
			std::vector<Engine::MS2HCDScoring::HCDInformationRecord> vect_hcd_score_records ; 
			std::vector<Engine::HornTransform::IsotopeFitRecord> vect_transform_records ; 
			std::vector<Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord> vect_combined_score_records ; 
			Engine::MS2CIDScoring::CIDInformationRecord this_cid_record ; 
			Engine::HornTransform::IsotopeFitRecord this_transform_record ; 
			Engine::MS2HCDScoring::HCDInformationRecord this_hcd_record ; 
			Engine::MS2CIDHCDCombinedScoring::CIDHCDCombinedInformationRecord this_combined_record ; 
			
		
			// while the thresholded parameter is already set in the clsPeakProcessParameters, we would
			// like to override that here if the data type is Finnigan o                            
			bool thresholded = true; 
			
			mint_percent_done = 0 ;
			menm_state = RUNNING ;

			// Create a RawData object and read data in
			char file_name_ch[256] ;
			GlypID::Utils::GetStr(file_name, file_name_ch) ;
			// enumerations of file type are the same in Readers namespace			
			raw_data = Engine::Readers::ReaderFactory::GetRawData((Engine::Readers::FileType)file_type, file_name_ch) ;
			if (raw_data == NULL)
			{
				throw new System::Exception(System::String::Concat(S"Could not open raw file: ", file_name));
			}
			
			if (search_fasta)
			{
				// Now read in fasta file ; 
				char fasta_name_ch[256] ; 
				GlypID::Utils::GetStr(protein_file_name, fasta_name_ch) ;
				bool readfasta ; 
				readfasta = Engine::Readers::TFastaFile::ReadFastaFile(fasta_name_ch, vect_candidate_proteins );  

				// Load up glycan list
				char glycan_name_ch[256] ;
				GlypID::Utils::GetStr(glycan_composition_file_name, glycan_name_ch) ; 
				glycan_data = new Engine::Readers::GlycanIo(); 
				glycan_data->LoadGlycanListFromFile(glycan_name_ch, vect_glycan_compositions) ; 
			}

			//Init
			parent_peak_processor = new Engine::PeakProcessing::PeakProcessor() ;
			cid_peak_processor = new Engine::PeakProcessing::PeakProcessor() ; 					
			hcd_peak_processor = new Engine::PeakProcessing::PeakProcessor() ; 			
			mass_transform = new Engine::HornTransform::MassTransform() ;
			hcd_scoring = new Engine::MS2HCDScoring::HCDScoring() ; 
			cid_scoring = new Engine::MS2CIDScoring::CIDScoring() ; 
			glycan_processor = new Engine::GlycanCompositionManager::GlycanProcessor() ; 
			comb_scoring_results = new Engine::Writers::CIDHCDCombinedScoringResults() ; 

			
			// Set parameters for discovering peaks. intensity threshold is set below.
			parent_peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
			cid_peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
			hcd_peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
			// Set parameters for horn transform (or mass transform) of the precursor	
			mass_transform->SetElementalIsotopeComposition(*transform_parameters->get_ElementIsotopeComposition()->mobjAtomicInfo) ;
			if (transform_parameters != NULL)
			{
				char averagine_formula[512] ;
				char tag_formula[512] ;
				mass_transform->SetOptions(transform_parameters->get_MaxCharge(), transform_parameters->get_MaxMW(), transform_parameters->get_MaxFit(),
					transform_parameters->get_MinS2N(), transform_parameters->get_CCMass(),transform_parameters->get_DeleteIntensityThreshold(),
					transform_parameters->get_MinIntensityForScore(), transform_parameters->get_NumPeaksForShoulder(),
					transform_parameters->get_CheckAllPatternsAgainstCharge1(), transform_parameters->get_UseMercuryCaching(), transform_parameters->get_O16O18Media()) ;

				averagine_formula[0] = '\0' ;
				tag_formula[0] = '\0' ;

				GlypID::Utils::GetStr(transform_parameters->get_AveragineFormula(), averagine_formula) ;
				if (transform_parameters->get_TagFormula() != NULL)
				{
					GlypID::Utils::GetStr(transform_parameters->get_TagFormula(), tag_formula) ;
				}
				mass_transform->SetIsotopeFitOptions(averagine_formula, tag_formula, transform_parameters->get_ThrashOrNot(), transform_parameters->get_CompleteFit()) ;
				mass_transform->SetIsotopeFitType((Engine::HornTransform::IsotopicFittingType) transform_parameters->get_IsotopeFitType()) ; 
			}			
			// Set parameters for scoring
			hcd_scoring->SetOptions(scoring_parameters->get_MinNumPeaksToConsider(), scoring_parameters->get_MinHCDMz(), scoring_parameters->get_MaxHCDMz()) ;			
			cid_scoring->InitPTable() ; 
	
			// Scan Range
			int min_scan = 1 ; 
			if (transform_parameters->get_UseScanRange() && transform_parameters->get_MinScan() >1)
				min_scan = transform_parameters->get_MinScan() ; 
			if (min_scan < raw_data->GetFirstScanNum())
				min_scan = raw_data->GetFirstScanNum() ; 			
			int max_scan = raw_data->GetLastScanNum() ; 
			if (transform_parameters->get_UseScanRange() && transform_parameters->get_MaxScan() < max_scan)
				max_scan = transform_parameters->get_MaxScan() ; 		
			
			short scan_ms_level = 1 ;			
			int scan_num = min_scan ; 

			// Prelim peptide processing
			std::vector<double> peptide_masses ; 	
			if (search_fasta)
			{
				glycan_processor->GetNGlycopeptideSeqList(vect_nglyco_peptides, vect_candidate_proteins);
				vect_candidate_proteins.clear() ; 
			//	glycan_processor->GetSequenceMassesFromProteins(vect_candidate_proteins, peptide_masses) ; 
				if(scoring_parameters->get_UsePPM())
					glycan_processor->SetPPMError(scoring_parameters->get_PPMTolerance()) ; 
				else
					throw new System::Exception(System::String::Concat(S"Search parameters should be in PPM for HCD-CID combined scoring: ", file_name));
					
					
			}


			// precursor details
			int parent_scan = 0  ; 
			double parent_mz = -1.0 ; 
			std::vector<double> parent_vect_mzs ;
			std::vector<double> parent_vect_intensities;
			std::map<double, int> parent_present ; 			
			typedef pair <double, int> parent_pair;
			std::map <double, int > ::const_iterator parent_present_iter ; 
			// msn
			std::vector<double> msn_vect_mzs ; 
			std::vector<double> msn_vect_intensities ; 	
			
			
			while(scan_num <= max_scan)
			{	
				mint_percent_done = (int)((scan_num - min_scan)/(max_scan - min_scan)) * 100 ; 				
				if (!raw_data->IsMSScan(scan_num))
				{
					// ------------ MSn peak processing -----------//
					msn_vect_intensities.clear();
					msn_vect_mzs.clear() ; 					

					raw_data->GetRawData(&msn_vect_mzs, &msn_vect_intensities, scan_num);					
					
					// get parent_scan and mz
					parent_scan = raw_data->GetParentScan(scan_num) ; 
					parent_mz = raw_data->GetParentMz(scan_num) ; 
					parent_present_iter = parent_present.find(parent_mz) ; 
					if (parent_present_iter != parent_present.end())
					{
						// parent is present in map so the CID should have been before
						// and parent should be already mass transformed
						if(raw_data->IsHCDScan(scan_num))
						{
							//----------- HCD processing -----------//
							// Anoop testing double hcd_background_intensity = GlypID::Utils::GetAverage(msn_vect_intensities, FLT_MAX) ;
							double hcd_background_intensity = GlypID::Utils::GetAverage(msn_vect_intensities, msn_vect_mzs, scoring_parameters->get_MinHCDMz(), scoring_parameters->get_MaxHCDMz()) ; 
							hcd_peak_processor->SetPeakIntensityThreshold(hcd_background_intensity) ; 
							hcd_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(scan_num)) ;
							int numHCDPeaks = hcd_peak_processor->DiscoverPeaks(&msn_vect_mzs, &msn_vect_intensities) ;  
							hcd_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;								
							// get HCD score
							std::vector <int> vectPeakIndices ; 
							Engine::MS2HCDScoring::GLYCAN_TYPE glycanType ; 							
							double hcd_score_value = hcd_scoring->DetermineGlycanType(*hcd_peak_processor->mobj_peak_data, vectPeakIndices, glycanType) ; 
							for (int k = 0; k < (int)vect_cid_score_records.size(); k++)
							{
								this_cid_record = vect_cid_score_records[k] ; 
								if (this_cid_record.mint_parent_scan_num == parent_scan && this_cid_record.mdbl_parent_mz == parent_mz)
								{								
									// Get y1 ion
									double y1_mz = 0.0 ; 									
									if (search_fasta)
									{
										if (hcd_score_value < 1) // Feb 2011, Slight speed up.
										{	
											double parent_time = raw_data->GetScanTime(this_cid_record.mint_parent_scan_num) ; 
											double cid_time = raw_data->GetScanTime(this_cid_record.mint_msn_scan_num) ;
											double hcd_time = raw_data->GetScanTime(scan_num); 
											if (scoring_parameters->IncludeBestHit)
											{
												glycan_processor->SearchAndFillPeptideInformationSingleRecord(vect_glycan_compositions,vect_nglyco_peptides, this_cid_record, glycanType, scoring_parameters->get_FilterCombList()) ; //, scoring_parameters->get_IncludeBestHit()) ; 
												y1_mz = hcd_scoring->DetermineY1IonThroughPeptideSearching(*hcd_peak_processor->mobj_peak_data, this_cid_record) ; 

												
												// Add everything to record
												this_combined_record.AddInfoToCombinedRecord(parent_scan, this_cid_record.mint_msn_scan_num, scan_num, parent_mz, this_cid_record.mdbl_mono_mz,
													this_cid_record.mint_parent_intensity, this_cid_record.mint_mono_intensity, this_cid_record.mshort_cs, this_cid_record.mdbl_fit, this_cid_record.mdbl_mono_mw,
													this_cid_record.mdbl_cid_score,this_cid_record.mdbl_cid_score_p_value, hcd_score_value, this_cid_record.mbln_contains_oxonium_ion, vectPeakIndices, glycanType, y1_mz, parent_time, cid_time, hcd_time); 
												bool found_site ; 									
												if (this_cid_record.mstr_glyco_site == "yes")
													found_site = 1; 
												else
													found_site = 0 ; 
												this_combined_record.AddSearchInfoToCombinedRecord(this_cid_record.mstr_pro_seq_name, this_cid_record.mdbl_seq_mass, this_cid_record.mdbl_glycan_mass,
													this_cid_record.mstr_glycan_composition, this_cid_record.mstr_pep_seq_name, this_cid_record.mstr_nglyco_site, this_cid_record.mdbl_mass_error, found_site) ; 									

												// Store record
												// Anoop Jan 2011: Added filtering here for both HCD and CID scoring
												if (((int)this_cid_record.mdbl_cid_score >= scoring_parameters->get_MinPathLength()) && (hcd_score_value < 1))
												{
													vect_combined_score_records.push_back(this_combined_record) ; 
												}
												
											}
											else
											{
												std::vector<Engine::MS2CIDScoring::CIDInformationRecord> temp_cid_records ; 
												temp_cid_records.push_back(this_cid_record) ; 
												glycan_processor->SearchAndFillPeptideInformation(vect_glycan_compositions, vect_nglyco_peptides, temp_cid_records, glycanType,true, false) ;
												
												for (int  ss =0 ; ss < temp_cid_records.size() ; ss++)
												{
													this_cid_record = temp_cid_records[ss] ; 
													y1_mz = hcd_scoring->DetermineY1IonThroughPeptideSearching(*hcd_peak_processor->mobj_peak_data, this_cid_record) ; 
													
													// Add everything to record
													this_combined_record.AddInfoToCombinedRecord(parent_scan, this_cid_record.mint_msn_scan_num, scan_num, parent_mz, this_cid_record.mdbl_mono_mz,
														this_cid_record.mint_parent_intensity, this_cid_record.mint_mono_intensity, this_cid_record.mshort_cs, this_cid_record.mdbl_fit, this_cid_record.mdbl_mono_mw,
														this_cid_record.mdbl_cid_score,this_cid_record.mdbl_cid_score_p_value, hcd_score_value, this_cid_record.mbln_contains_oxonium_ion, vectPeakIndices, glycanType, y1_mz, parent_time, cid_time, hcd_time); 
													bool found_site ; 									
													if (this_cid_record.mstr_glyco_site == "yes")
														found_site = 1; 
													else
														found_site = 0 ; 
													this_combined_record.AddSearchInfoToCombinedRecord(this_cid_record.mstr_pro_seq_name, this_cid_record.mdbl_seq_mass, this_cid_record.mdbl_glycan_mass,
														this_cid_record.mstr_glycan_composition, this_cid_record.mstr_pep_seq_name, this_cid_record.mstr_nglyco_site, this_cid_record.mdbl_mass_error, found_site) ; 									

													// Store record
													// Anoop Jan 2011: Added filtering here for both HCD and CID scoring
													if (((int)this_cid_record.mdbl_cid_score >= scoring_parameters->get_MinPathLength()) && (hcd_score_value < 1))
													{
														vect_combined_score_records.push_back(this_combined_record) ; 
													}

												}

												temp_cid_records.clear() ; 
											}
										}
									}
									else
									{
									  y1_mz = hcd_scoring->DetermineY1IonThroughCorrelation(*hcd_peak_processor->mobj_peak_data, this_cid_record.m_peak_data, 0.3, 2000) ; 
									}
										
									

								}
							}
							
						}
						else
						{
							// Ignore the scan, could be ETD among other things

						
												
						}
					}
					else
					{
						if (raw_data->IsCIDScan(scan_num))
						{
							// CID scan so process it along with the parent

							// ------------- Parent Processing ------------- //						
							parent_vect_intensities.clear();
							parent_vect_mzs.clear() ; 
														
							raw_data->GetRawData(&parent_vect_mzs, &parent_vect_intensities, parent_scan);					
							// set noise floor
							double thres =  GlypID::Utils::GetAverage(parent_vect_intensities, FLT_MAX) ; 
							double background_intensity = GlypID::Utils::GetAverage(parent_vect_intensities, (float)(5*thres)) ;
							parent_peak_processor->SetPeakIntensityThreshold(background_intensity * peak_parameters->get_PeakBackgroundRatio()) ; 
							parent_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(parent_scan)) ;			
							// discover peaks
							int numPeaks = parent_peak_processor->DiscoverPeaks(&parent_vect_mzs, &parent_vect_intensities) ; 
							parent_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;	
							// Invariant : discovering the parent peak has to be done before transforming the precursor distribution.  Since, THRASH removes the distribution from the spectra, 
							// trying to discover the parent peak after will result in error/
							Engine::PeakProcessing::Peak parentPeak ; 
							double parent_match = parent_peak_processor->GetClosestPeakMz(parent_mz, parentPeak) ; 					

							// Invariant : Keeping the recorded mono m/z as backup, again needs to be done before mass transforming
							Engine::PeakProcessing::Peak monoPeak ; 
							double mono_mz = raw_data->GetMonoMZFromHeader(scan_num) ; 
							double mono_mz_match = parent_peak_processor->GetClosestPeakMz(mono_mz, monoPeak) ; 


							// If summing is turned on, repeat the process							
							if (transform_parameters->get_SumSpectra())
							{
								// back up the spectra so that we can get the original mono intensity 
								temp_parent_peak_processor = new Engine::PeakProcessing::PeakProcessor() ; 
								temp_parent_peak_processor = parent_peak_processor ; 								
								parent_vect_intensities.clear();
								parent_vect_mzs.clear() ; 
								parent_peak_processor->mobj_peak_data->Clear() ; 
								raw_data->GetSummedSpectra(&parent_vect_mzs, &parent_vect_intensities, parent_scan, transform_parameters->get_NumScansToSumOver(), transform_parameters->get_MinMZ(), transform_parameters->get_MaxMZ() );					
								thres =  GlypID::Utils::GetAverage(parent_vect_intensities, FLT_MAX) ; 
								background_intensity = GlypID::Utils::GetAverage(parent_vect_intensities, (float)(5*thres)) ;
								parent_peak_processor->SetPeakIntensityThreshold(background_intensity * peak_parameters->get_PeakBackgroundRatio()) ; 
								parent_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(parent_scan)) ;			
								// discover peaks
								int numPeaks = parent_peak_processor->DiscoverPeaks(&parent_vect_mzs, &parent_vect_intensities) ; 
								parent_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;
							}
							

							// mass transform
							// Set  noise floor parameters for isotopic peaks (calculated from the peaks discovered)
							int numDeisotoped = 0 ;
							double min_peptide_intensity = background_intensity * transform_parameters->get_PeptideMinBackgroundRatio() ;
							if (transform_parameters->get_UseAbsolutePeptideIntensity())
							{
								if (min_peptide_intensity < transform_parameters->get_AbsolutePeptideIntensity())
										min_peptide_intensity = transform_parameters->get_AbsolutePeptideIntensity() ;
							}
							mass_transform->Reset() ; 
							vect_transform_records.clear() ; 
							bool found_transform = mass_transform->FindPrecursorTransform(*parent_peak_processor->mobj_peak_data, parent_mz, &vect_transform_records, background_intensity, min_peptide_intensity) ;				
							if (!found_transform)
							{
								// this could be from a bad spectra
								// Look in header		
 								short cs = raw_data->GetMonoChargeFromHeader(scan_num) ; 
								if (cs > 0)
								{
									//create transform record with  instrument recorded charge cs and mass
									this_transform_record.mshort_cs = cs ; 									
									this_transform_record.mdbl_mono_mw = (mono_mz - mdbl_cc_mass) * cs ; 
									this_transform_record.mdbl_average_mw = 0 ; // since no composition
									this_transform_record.mdbl_fit = 1; 
									if (mono_mz_match == 0) //couldnt find the parent, so set intensity to noise level
										this_transform_record.mint_abundance = (int) min_peptide_intensity ; 
									else
										this_transform_record.mint_abundance = (int) monoPeak.mdbl_intensity ; 
									vect_transform_records.push_back(this_transform_record) ; 
								}
								else
								{
									// instrument couldnt find charge as well, so create dummy transforms with cs 2 and 3, if user has agreed to do this
									if (scoring_parameters->get_AllocateDefaultCharges())
									{
										for( int charge = 2 ; charge <= 3 ; charge++)
										{
											this_transform_record.mshort_cs = charge ; 
											this_transform_record.mdbl_mono_mw = (parent_mz - mdbl_cc_mass) * charge ; 
											this_transform_record.mdbl_average_mw = 0 ; 
											this_transform_record.mdbl_fit = 1; 
											if (parent_match == 0) //couldnt find the parent, so set intensity to noise level
												this_transform_record.mint_abundance = (int) min_peptide_intensity ; 
											else
												this_transform_record.mint_abundance = (int) parentPeak.mdbl_intensity  ; 
											vect_transform_records.push_back(this_transform_record) ; 
										}
									}
								}
							}
							else
							{
								// If summing is enabled, replace the detected mono intensity with parent intensity
								if(transform_parameters->get_SumSpectra())
								{
									
									for (int k = 0; k < (int)vect_transform_records.size(); k++)
									{		
										Engine::PeakProcessing::Peak tempPeak ;
										double mono_match = temp_parent_peak_processor->GetClosestPeakMz(vect_transform_records[k].mdbl_mz, tempPeak) ; 
										if (tempPeak.mdbl_intensity >  min_peptide_intensity)
										{
											vect_transform_records[k].mint_abundance = (int) tempPeak.mdbl_intensity ;
											vect_transform_records[k].mint_mono_intensity =(int)tempPeak.mdbl_intensity ;
										}
										else
										{
											// In case if precursor is in the noise floor. 
											vect_transform_records[k].mint_abundance = (int) min_peptide_intensity ; 
											vect_transform_records[k].mint_mono_intensity = (int) min_peptide_intensity ;
										}
									}
								}
							}

							//------------- CID processing -------------//
							double msn_thres =  GlypID::Utils::GetAverage(msn_vect_intensities, FLT_MAX) ;
							background_intensity = 0 ; //GlypID::Utils::GetAverage(hcd_vect_intensities, (float)(5*thres)) ;
							cid_peak_processor->SetPeakIntensityThreshold(background_intensity * peak_parameters->get_PeakBackgroundRatio()) ; 							
							cid_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(scan_num)) ;
							int numCIDPeaks = cid_peak_processor->DiscoverPeaks(&msn_vect_mzs, &msn_vect_intensities) ; 
							cid_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;		
							// store in map
							parent_present.insert(parent_pair (parent_mz, 1)) ; 																
							for (int k = 0; k < (int)vect_transform_records.size(); k++)
							{
								this_transform_record = vect_transform_records[k] ; 
								// get score
								double cid_score_value = cid_scoring->CalculateCIDScore(*cid_peak_processor->mobj_peak_data, this_transform_record.mshort_cs ) ; 
								// get p value
								double p_value = cid_scoring->CalculateCIDScorePValue((int) cid_score_value) ; 

								// see if oxonium is present
								bool contains_oxonium = cid_scoring->CheckForOxoniumIons(*cid_peak_processor->mobj_peak_data) ; 
								//get mono mz
								double mono_mz = (this_transform_record.mdbl_mono_mw/this_transform_record.mshort_cs) + mdbl_cc_mass  ; 

								
								this_cid_record.Clear() ; 
								this_cid_record.AddInfoToCIDRecord(scan_num, parent_scan, scan_ms_level, 1, parent_mz, mono_mz, 0, this_transform_record.mint_abundance, 
									this_transform_record.mshort_cs, this_transform_record.mdbl_fit, this_transform_record.mdbl_mono_mw, this_transform_record.mdbl_average_mw, 
									cid_score_value, p_value, contains_oxonium, 0 , 0,0) ; 
								this_cid_record.AddPeaks(*cid_peak_processor->mobj_peak_data) ; 
								//store into cid record
								vect_cid_score_records.push_back(this_cid_record) ; 
								
							}		
						}
					}					
				}
				else
				{
					parent_present.clear() ; 
				}
				scan_num++ ; 
				Console::WriteLine(scan_num) ;
				
			}

		
			//----------- Store and return ----------//			
			comb_scoring_results->AddMultipleCombinedRecords(vect_combined_score_records) ; 
			all_results->SetCombinedScoringResults(comb_scoring_results) ; 

			mint_percent_done = 100 ; 
		}
		catch (char *mesg)
		{
			if (comb_scoring_results != NULL)
			{
				delete comb_scoring_results ; 
				comb_scoring_results = NULL ;
			}			
			if (parent_peak_processor != NULL) 
			{
				delete parent_peak_processor ;
				parent_peak_processor = NULL ; 
			}
			if (cid_peak_processor != NULL) 
			{
				delete cid_peak_processor ;
				cid_peak_processor = NULL ; 
			}
			if (hcd_peak_processor != NULL) 
			{
				delete hcd_peak_processor ;
				hcd_peak_processor = NULL ; 
			}
			if( raw_data != NULL)
			{
				raw_data->Close() ; 
				delete raw_data ; 
			}
			if (mass_transform != NULL)
			{
				delete mass_transform ;
				mass_transform = NULL ; 
			}
			if(hcd_scoring != NULL)
			{
				delete hcd_scoring ; 
				hcd_scoring = NULL ; 
			}
			if(cid_scoring != NULL)
			{
				delete cid_scoring ; 
				cid_scoring = NULL ; 
			}
			if (glycan_data != NULL)
			{
				delete glycan_data ; 
				glycan_data = NULL ; 
			}
			if (glycan_processor != NULL)
			{
				delete glycan_processor  ;
				glycan_processor = NULL ; 
			}
			System::String *exception_msg = new System::String(mesg) ; 
			menm_state = enmProcessState::ERROR ;
			throw new System::Exception(exception_msg) ; 
		}			
		if (parent_peak_processor != NULL) 
		{
			delete parent_peak_processor ;
			parent_peak_processor = NULL ; 
		}
		if (cid_peak_processor != NULL) 
		{
			delete cid_peak_processor ;
			cid_peak_processor = NULL ; 
		}
		if (hcd_peak_processor != NULL) 
		{
			delete hcd_peak_processor ;
			hcd_peak_processor = NULL ; 
		}
		if( raw_data != NULL)
		{
			raw_data->Close() ; 
			delete raw_data ; 
		}
		if (mass_transform != NULL)
		{
			delete mass_transform ;
			mass_transform = NULL ; 
		}
		if(hcd_scoring != NULL)
		{
			delete hcd_scoring ; 
			hcd_scoring = NULL ; 
		}
		if(cid_scoring != NULL)
		{
			delete cid_scoring ; 
			cid_scoring = NULL ; 
		}
		if (glycan_data != NULL)
		{
			delete glycan_data ; 
			glycan_data = NULL ; 
		}
		if (glycan_processor != NULL)
		{
			delete glycan_processor  ;
			glycan_processor = NULL ; 
		}

		menm_state = COMPLETE ;
		return all_results ;
	}


	GlypID::Results::clsCIDETDCombinedScoringResults*  clsProcRunner::CreateCIDETDScoringResults(System::String *file_name, GlypID::Readers::FileType file_type, 
		System::String *protein_file_name, System::String *glycan_composition_file_name , 
		GlypID::Peaks::clsPeakProcessorParameters *peak_parameters,
		GlypID::HornTransform::clsHornTransformParameters *transform_parameters, GlypID::Scoring::clsScoringParameters *scoring_parameters, bool search_fasta)
	{
		GlypID::Results::clsCIDETDCombinedScoringResults *all_results = __gc new GlypID::Results::clsCIDETDCombinedScoringResults() ; 



	/*	if (menm_state == RUNNING)
		{
			throw new System::Exception(S"Process already running in clsProcRunner. Cannot run two processes with same object");
		}

		if (file_name == NULL || file_name == S"")
		{
			throw new System::Exception(S"Please enter a file name to process");
		}
		if (peak_parameters == NULL)
		{
			throw new System::Exception(S"Please specify peak processing parameters.");
		}
		if (transform_parameters == NULL)
		{
			throw new System::Exception(S"Please specify mass transform parameters.");
		}

		if (search_fasta)
		{
			if (protein_file_name == NULL || protein_file_name == S"")
			{
				throw new System::Exception(S"Please enter a protein list file") ; 
			}
			if (glycan_composition_file_name == NULL || glycan_composition_file_name == S"")
			{
				throw new System::Exception(S"Please enter a glycan composition file") ; 
			}
		}

		// Declarations
		Engine::Readers::RawData __nogc *raw_data = NULL ; 		
		Engine::Readers::GlycanIo  __nogc *glycan_data = NULL ; 
		Engine::PeakProcessing::PeakProcessor __nogc *parent_peak_processor = NULL ;	
		Engine::PeakProcessing::PeakProcessor __nogc *temp_parent_peak_processor = NULL ;	
		Engine::PeakProcessing::PeakProcessor __nogc *cid_peak_processor = NULL ;				
		Engine::PeakProcessing::PeakProcessor __nogc *etd_peak_processor = NULL ;				
		Engine::HornTransform::MassTransform __nogc *mass_transform = NULL ; 
		Engine::MS2CIDScoring::CIDScoring __nogc *cid_scoring = NULL ; 				
		Engine::MS2ETDScoring::ETDScoring __nogc *etd_scoring = NULL ; 
		Engine::GlycanCompositionManager::GlycanProcessor __nogc *glycan_processor = NULL ; 

		MwtWinDll::MolecularWeightCalculator *mwt =  __gc new MwtWinDll::MolecularWeightCalculator()  ;
		MwtWinDll::MWPeptideClass::udtFragmentationSpectrumOptionsType *udtFragSpectrumOptions ; 

									

		Engine::Writers::CIDETDCombinedScoringResults __nogc *comb_scoring_results = NULL ; 
		
		try
		{
			std::vector<Engine::SequenceManager::Sequence> vect_candidate_proteins ; 
			std::vector<Engine::GlycanTheoretical::GlycanComposition> vect_glycan_compositions ; 

			std::vector<Engine::MS2CIDScoring::CIDInformationRecord> vect_cid_score_records ; 			
			std::vector<Engine::MS2CIDETDCombinedScoring::CIDETDCombinedInformationRecord> vect_comb_records ; 
			std::vector<Engine::HornTransform::IsotopeFitRecord> vect_transform_records ; 			
			Engine::MS2CIDScoring::CIDInformationRecord this_cid_record ; 
			Engine::HornTransform::IsotopeFitRecord this_transform_record ; 
			Engine::MS2CIDETDCombinedScoring::CIDETDCombinedInformationRecord this_combined_record ; 

			// while the thresholded parameter is already set in the clsPeakProcessParameters, we would
			// like to override that here if the data type is Finnigan o                            
			bool thresholded = true; 
			
			mint_percent_done = 0 ;
			menm_state = RUNNING ;

			// Create a RawData object and read data in
			char file_name_ch[256] ;
			GlypID::Utils::GetStr(file_name, file_name_ch) ;
			// enumerations of file type are the same in Readers namespace			
			raw_data = Engine::Readers::ReaderFactory::GetRawData((Engine::Readers::FileType)file_type, file_name_ch) ;
			if (raw_data == NULL)
			{
				throw new System::Exception(System::String::Concat(S"Could not open raw file: ", file_name));
			}
			
			if (search_fasta)
			{
				// Now read in fasta file ; 
				char fasta_name_ch[256] ; 
				GlypID::Utils::GetStr(protein_file_name, fasta_name_ch) ;
				bool readfasta ; 
				readfasta = Engine::Readers::TFastaFile::ReadFastaFile(fasta_name_ch, vect_candidate_proteins );  

				// Load up glycan list
				char glycan_name_ch[256] ;
				GlypID::Utils::GetStr(glycan_composition_file_name, glycan_name_ch) ; 
				glycan_data = new Engine::Readers::GlycanIo(); 
				glycan_data->LoadGlycanListFromFile(glycan_name_ch, vect_glycan_compositions) ; 
			}

			//Init
			parent_peak_processor = new Engine::PeakProcessing::PeakProcessor() ;
			cid_peak_processor = new Engine::PeakProcessing::PeakProcessor() ; 					
			etd_peak_processor = new Engine::PeakProcessing::PeakProcessor() ; 			
			mass_transform = new Engine::HornTransform::MassTransform() ;
			etd_scoring = new Engine::MS2ETDScoring::ETDScoring() ; 
			cid_scoring = new Engine::MS2CIDScoring::CIDScoring() ; 
			glycan_processor = new Engine::GlycanCompositionManager::GlycanProcessor() ; 
			comb_scoring_results = new Engine::Writers::CIDETDCombinedScoringResults() ; 

			// Set parameters for discovering peaks. intensity threshold is set below.
			parent_peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
			cid_peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
			etd_peak_processor->SetOptions(peak_parameters->get_SignalToNoiseThreshold(), 0, thresholded, (Engine::PeakProcessing::PEAK_FIT_TYPE)peak_parameters->get_PeakFitType()) ;
			
			// Set parameters for horn transform (or mass transform) of the precursor	
			mass_transform->SetElementalIsotopeComposition(*transform_parameters->get_ElementIsotopeComposition()->mobjAtomicInfo) ;
			if (transform_parameters != NULL)
			{
				char averagine_formula[512] ;
				char tag_formula[512] ;
				mass_transform->SetOptions(transform_parameters->get_MaxCharge(), transform_parameters->get_MaxMW(), transform_parameters->get_MaxFit(),
					transform_parameters->get_MinS2N(), transform_parameters->get_CCMass(),transform_parameters->get_DeleteIntensityThreshold(),
					transform_parameters->get_MinIntensityForScore(), transform_parameters->get_NumPeaksForShoulder(),
					transform_parameters->get_CheckAllPatternsAgainstCharge1(), transform_parameters->get_UseMercuryCaching(), transform_parameters->get_O16O18Media()) ;

				averagine_formula[0] = '\0' ;
				tag_formula[0] = '\0' ;

				GlypID::Utils::GetStr(transform_parameters->get_AveragineFormula(), averagine_formula) ;
				if (transform_parameters->get_TagFormula() != NULL)
				{
					GlypID::Utils::GetStr(transform_parameters->get_TagFormula(), tag_formula) ;
				}
				mass_transform->SetIsotopeFitOptions(averagine_formula, tag_formula, transform_parameters->get_ThrashOrNot(), transform_parameters->get_CompleteFit()) ;
				mass_transform->SetIsotopeFitType((Engine::HornTransform::IsotopicFittingType) transform_parameters->get_IsotopeFitType()) ; 
			}			
	
			// Set parameters for scoring
			// blah etd_scoring->SetOptions(scoring_parameters->get_MinNumPeaksToConsider(), scoring_parameters->get_MinHCDMz(), scoring_parameters->get_MaxHCDMz()) ;			
			cid_scoring->InitPTable() ; 

			// Set parameters for MolecularWeightCalculatoer DLL
			mwt->SetElementMode(MwtWinDll::MWElementAndMassRoutines::emElementModeConstants::emIsotopicMass) ; 
            //udtFragSpectrumOptions->Initialize();
            udtFragSpectrumOptions = &mwt->Peptide->GetFragmentationSpectrumOptions();
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itAIon].ShowIon = false;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itBIon].ShowIon = true;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itYIon].ShowIon = true;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itCIon].ShowIon = true;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itZIon].ShowIon = true;

			    
			// Scan Range
			int min_scan = 1 ; 
			if (transform_parameters->get_UseScanRange() && transform_parameters->get_MinScan() >1)
				min_scan = transform_parameters->get_MinScan() ; 
			if (min_scan < raw_data->GetFirstScanNum())
				min_scan = raw_data->GetFirstScanNum() ; 			
			int max_scan = raw_data->GetLastScanNum() ; 
			if (transform_parameters->get_UseScanRange() && transform_parameters->get_MaxScan() < max_scan)
				max_scan = transform_parameters->get_MaxScan() ; 		
			
			short scan_ms_level = 1 ;			
			int scan_num = min_scan ; 

			// Prelim peptide processing
			std::vector<double> peptide_masses ; 	
			if (search_fasta)
			{
				glycan_processor->GetSequenceMassesFromProteins(vect_candidate_proteins, peptide_masses) ; 
				if(scoring_parameters->get_UsePPM())
					glycan_processor->SetPPMError(scoring_parameters->get_PPMTolerance()) ; 
				else
					throw new System::Exception(System::String::Concat(S"Search parameters should be in PPM for HCD-CID combined scoring: ", file_name));
			}
			
			// precursor details
			int parent_scan = 0  ; 
			double parent_mz = -1.0 ; 
			std::vector<double> parent_vect_mzs ;
			std::vector<double> parent_vect_intensities;
			std::map<double, int> parent_present ; 			
			typedef pair <double, int> parent_pair;
			std::map <double, int > ::const_iterator parent_present_iter ; 
			// msn
			std::vector<double> msn_vect_mzs ; 
			std::vector<double> msn_vect_intensities ; 	
			

			// Start process
			// Natural assumption ETD follows CID					
			while(scan_num <= max_scan)
			{	
				mint_percent_done = (int)((scan_num - min_scan)/(max_scan - min_scan)) * 100 ; 				
				if (!raw_data->IsMSScan(scan_num))
				{
					// ------------ MSn peak processing -----------//
					msn_vect_intensities.clear();
					msn_vect_mzs.clear() ; 					

					raw_data->GetRawData(&msn_vect_mzs, &msn_vect_intensities, scan_num);					
					
					// get parent_scan and mz
					parent_scan = raw_data->GetParentScan(scan_num) ; 
					parent_mz = raw_data->GetParentMz(scan_num) ; 
					parent_present_iter = parent_present.find(parent_mz) ; 
					if (parent_present_iter != parent_present.end())
					{
						// parent is present in map so the CID should have been before
						// and parent should be already mass transformed
						if(raw_data->IsETDScan(scan_num))
						{
							 //---------- ETD scoring ---------//

							// Get peaks
							double msn_thres =  GlypID::Utils::GetAverage(msn_vect_intensities, FLT_MAX) ;
							double background_intensity = 0 ; //GlypID::Utils::GetAverage(msn_vect_intensities, (float)(2*msn_thres)) ;
							etd_peak_processor->Clear() ; 
							this_combined_record.Clear() ; 
							etd_peak_processor->SetPeakIntensityThreshold(background_intensity) ; 							
							etd_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(scan_num)) ;
							
							int numETDPeaks = etd_peak_processor->DiscoverPeaks(&msn_vect_mzs, &msn_vect_intensities) ; 
							etd_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;	

							std::vector <double> obs_peak_mzs ; 
							std::vector <double> obs_peak_intensities ;
							for (int i=0 ; i <200 ; i++) // Need to parameterize this
							{
								Engine::PeakProcessing::Peak thisPeak ; 
								bool isPeak = etd_peak_processor->mobj_peak_data->GetNextPeak(300, 2000, thisPeak) ; 
								if (isPeak)
								{
									obs_peak_mzs.push_back(thisPeak.mdbl_mz) ; 
									obs_peak_intensities.push_back(thisPeak.mdbl_intensity) ; 
								}
							}


							// Scoring
							double max_etd_score = 0; 
							for (int k = 0; k < (int)vect_cid_score_records.size(); k++)
							{
								this_cid_record = vect_cid_score_records[k] ; 
								double glycan_mass = this_cid_record.mdbl_glycan_mass ; 
								if (glycan_mass > 0 )
								{
									// glycopeptide (s) has been matched
									if (this_cid_record.mint_parent_scan_num == parent_scan && this_cid_record.mdbl_parent_mz == parent_mz)
									{								
											
										double elution_time = 0.0 ; 
										elution_time = raw_data->GetScanTime(this_cid_record.mint_parent_scan_num) ; 

										// Update peak viewing options
										switch(this_cid_record.mshort_cs)
										{
										case 1:  
											udtFragSpectrumOptions->DoubleChargeIonsShow = false;
											udtFragSpectrumOptions->TripleChargeIonsShow = false;
											break ; 
										case 2:
											udtFragSpectrumOptions->DoubleChargeIonsShow = true;
											udtFragSpectrumOptions->TripleChargeIonsShow = false;
											break ;
										default:
											udtFragSpectrumOptions->DoubleChargeIonsShow = true;
											udtFragSpectrumOptions->TripleChargeIonsShow = true;
											break ; 
										}

										// Add modification to list of modifications in MWT
										System::String *symbol = "^" ; 
										char chr_symbol [20] ; 
										GlypID::Utils::GetStr(symbol, chr_symbol);
										System::String *comment = "Glycan";
										bool var = false ; 
										

										mwt->Peptide->SetModificationSymbol(symbol, glycan_mass, &var, &comment) ; 					
										// Indicate modification in sequence
										std::string sequence = "" ; 
										etd_scoring->AddGlycanModificationToSequence(chr_symbol, this_cid_record.mstr_pep_seq_name, sequence) ; 
										// Get theoretical spectra
										//mwt->Peptide->SetSequence(this_cid_record.mstr_pep_seq_name.c_str(), 
										mwt->Peptide->SetSequence(sequence.c_str(), 
											MwtWinDll::MWPeptideClass::ntgNTerminusGroupConstants::ntgHydrogen,
											MwtWinDll::MWPeptideClass::ctgCTerminusGroupConstants::ctgHydroxyl,
											false);										
										MwtWinDll::MWPeptideClass::udtFragmentationSpectrumDataType udtFrag[]; 
										mwt->Peptide->GetFragmentationMasses(&udtFrag) ; 
										std::vector<double> theor_mzs  ; 
										std::vector<double> theor_intensities ; 
										for (int j = 0 ; j < udtFrag->Length ; j++)
										{
											double thmz = udtFrag[j].Mass ; 
											double thintn = udtFrag[j].Intensity;
											theor_mzs.push_back(thmz) ; 
											theor_intensities.push_back(thintn) ; 
										}										

										// compare theoretical and observed
										//double etd_score = etd_scoring->CalculateETDScore(*etd_peak_processor->mobj_peak_data, theor_mzs, theor_intensities) ; 										
										double etd_score = etd_scoring->CalculateETDScore2(obs_peak_mzs, obs_peak_intensities, theor_mzs, theor_intensities) ;  
										/*Console::Write(etd_score) ; 
										Console::Write("\t") ; 
										Console::WriteLine(this_cid_record.mstr_pep_seq_name.c_str())		;*/		
										// Choose the best match if matching has been done
									/*    if (etd_score > max_etd_score)
										{
											max_etd_score = etd_score ; 
											this_combined_record.Clear(); 
											this_combined_record.AddInfoToCombinedRecord(this_cid_record.mint_parent_scan_num, this_cid_record.mint_msn_scan_num, scan_num, 
												this_cid_record.mdbl_parent_mz, this_cid_record.mdbl_mono_mz, this_cid_record.mint_parent_intensity, this_cid_record.mint_mono_intensity, 
												this_cid_record.mshort_cs, this_cid_record.mdbl_fit, this_cid_record.mdbl_mono_mw, this_cid_record.mdbl_cid_score, this_cid_record.mdbl_cid_score_p_value, 
												etd_score, 0, this_cid_record.mbln_contains_oxonium_ion, this_cid_record.mdbl_y1_mz, this_cid_record.mdbl_parent_scan_time) ;
											this_combined_record.AddSearchInfoToCombinedRecord(this_cid_record.mstr_pro_seq_name, this_cid_record.mdbl_seq_mass, this_cid_record.mdbl_glycan_mass, 
												this_cid_record.mstr_glycan_composition, this_cid_record.mstr_pep_seq_name, this_cid_record.mdbl_mass_error, true); 

										}
									}
									else
									{
										this_combined_record.Clear() ; 
										this_combined_record.AddInfoToCombinedRecord(this_cid_record.mint_parent_scan_num, this_cid_record.mint_msn_scan_num, scan_num, 
												this_cid_record.mdbl_parent_mz, this_cid_record.mdbl_mono_mz, this_cid_record.mint_parent_intensity, this_cid_record.mint_mono_intensity, 
												this_cid_record.mshort_cs, this_cid_record.mdbl_fit, this_cid_record.mdbl_mono_mw, this_cid_record.mdbl_cid_score, this_cid_record.mdbl_cid_score_p_value, 
												0, 0, this_cid_record.mbln_contains_oxonium_ion, this_cid_record.mdbl_y1_mz, this_cid_record.mdbl_parent_scan_time) ;
										this_combined_record.AddSearchInfoToCombinedRecord(this_cid_record.mstr_pro_seq_name, this_cid_record.mdbl_seq_mass, this_cid_record.mdbl_glycan_mass, 
												this_cid_record.mstr_glycan_composition, this_cid_record.mstr_pep_seq_name, this_cid_record.mdbl_mass_error, true); 

									}
								}
							}
							
							vect_cid_score_records.clear() ; 

							// FDR calculuation
							int num_false_positives = 0 ;
							int num_compositions = vect_glycan_compositions.size();
							for(int i=0; i < num_compositions ; i++)
							{
								double glycan_mass = (double) vect_glycan_compositions[i].averageMass() ; 
								System::String *symbol = "^" ; 
								char chr_symbol [20] ; 
								GlypID::Utils::GetStr(symbol, chr_symbol);
								System::String *comment = "Glycan";
								bool var = false ; 
										
								mwt->Peptide->SetModificationSymbol(symbol, glycan_mass, &var, &comment) ; 					
								// Indicate modification in sequence
								std::string sequence = "" ; 
								etd_scoring->AddGlycanModificationToSequence(chr_symbol, this_combined_record.mstr_pep_seq_name, sequence) ; 
								// Get theoretical spectra										
								mwt->Peptide->SetSequence(sequence.c_str(), 
									MwtWinDll::MWPeptideClass::ntgNTerminusGroupConstants::ntgHydrogen,
									MwtWinDll::MWPeptideClass::ctgCTerminusGroupConstants::ctgHydroxyl,
									false);										
								MwtWinDll::MWPeptideClass::udtFragmentationSpectrumDataType udtFrag[]; 
								mwt->Peptide->GetFragmentationMasses(&udtFrag) ; 
								std::vector<double> theor_mzs  ; 
								std::vector<double> theor_intensities ; 
								for (int j = 0 ; j < udtFrag->Length ; j++)
								{
									double thmz = udtFrag[j].Mass ; 
									double thintn = udtFrag[j].Intensity;
									theor_mzs.push_back(thmz) ; 
									theor_intensities.push_back(thintn) ; 
								}										

								// compare theoretical and observed
								//double etd_score = etd_scoring->CalculateETDScore(*etd_peak_processor->mobj_peak_data, theor_mzs, theor_intensities) ; 										
								double etd_score = etd_scoring->CalculateETDScore2(obs_peak_mzs, obs_peak_intensities, theor_mzs, theor_intensities) ;  
								if (etd_score > this_combined_record.mdbl_etd_score)
								{
									num_false_positives++ ; 
								}
							}
							double fdr = ((double)num_false_positives/num_compositions)*100  ; 
							this_combined_record.mdbl_etd_score_fdr = fdr ; 
						
							// Store record							
							if (((int)this_combined_record.mdbl_cid_score >= scoring_parameters->get_MinPathLength()) || this_combined_record.mdbl_etd_score_fdr <1)
							{
								vect_comb_records.push_back(this_combined_record) ; 
							}						
						}
					}
					else
					{
						if (raw_data->IsCIDScan(scan_num))
						{
							// CID scan, deisotope partne first, then score cid and do glycopeptide matching

							// ------------- Parent Processing ------------- //						
							parent_vect_intensities.clear();
							parent_vect_mzs.clear() ; 														
							raw_data->GetRawData(&parent_vect_mzs, &parent_vect_intensities, parent_scan);												
							// set noise floor
							double thres =  GlypID::Utils::GetAverage(parent_vect_intensities, FLT_MAX) ; 
							double background_intensity = GlypID::Utils::GetAverage(parent_vect_intensities, (float)(5*thres)) ;
							parent_peak_processor->SetPeakIntensityThreshold(background_intensity * peak_parameters->get_PeakBackgroundRatio()) ; 
							parent_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(parent_scan)) ;			
							// discover peaks
							int numPeaks = parent_peak_processor->DiscoverPeaks(&parent_vect_mzs, &parent_vect_intensities) ; 
							parent_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;	
							// Invariant : discovering the parent peak has to be done before transforming the precursor distribution.  Since, THRASH removes the distribution from the spectra, 
							// trying to discover the parent peak after will result in error/
							Engine::PeakProcessing::Peak parentPeak ; 
							double parent_match = parent_peak_processor->GetClosestPeakMz(parent_mz, parentPeak) ; 					

							// Invariant : Keeping the recorded mono m/z as backup, again needs to be done before mass transforming
							Engine::PeakProcessing::Peak monoPeak ; 
							double mono_mz = raw_data->GetMonoMZFromHeader(scan_num) ; 
							double mono_mz_match = parent_peak_processor->GetClosestPeakMz(mono_mz, monoPeak) ; 


							// If summing is turned on, repeat the peak finding process							
							if (transform_parameters->SumSpectra)
							{
								// back up the spectra so that we can get the original mono intensity 
								temp_parent_peak_processor = new Engine::PeakProcessing::PeakProcessor() ; 
								temp_parent_peak_processor = parent_peak_processor ; 								
								parent_vect_intensities.clear();
								parent_vect_mzs.clear() ; 
								parent_peak_processor->mobj_peak_data->Clear() ; 
								raw_data->GetSummedSpectra(&parent_vect_mzs, &parent_vect_intensities, parent_scan, transform_parameters->get_NumScansToSumOver(), transform_parameters->get_MinMZ(), transform_parameters->get_MaxMZ() );					
								thres =  GlypID::Utils::GetAverage(parent_vect_intensities, FLT_MAX) ; 
								background_intensity = GlypID::Utils::GetAverage(parent_vect_intensities, (float)(5*thres)) ;
								parent_peak_processor->SetPeakIntensityThreshold(background_intensity * peak_parameters->get_PeakBackgroundRatio()) ; 
								parent_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(parent_scan)) ;			
								// discover peaks
								int numPeaks = parent_peak_processor->DiscoverPeaks(&parent_vect_mzs, &parent_vect_intensities) ; 
								parent_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;
							}
							

							// 	mass transform
							// Set  noise floor parameters for isotopic peaks (calculated from the peaks discovered)
							int numDeisotoped = 0 ;
							double min_peptide_intensity = background_intensity * transform_parameters->get_PeptideMinBackgroundRatio() ;
							if (transform_parameters->get_UseAbsolutePeptideIntensity())
							{
								if (min_peptide_intensity < transform_parameters->get_AbsolutePeptideIntensity())
										min_peptide_intensity = transform_parameters->get_AbsolutePeptideIntensity() ;
							}
							mass_transform->Reset() ; 
							vect_transform_records.clear() ; 
							bool found_transform = mass_transform->FindPrecursorTransform(*parent_peak_processor->mobj_peak_data, parent_mz, &vect_transform_records, background_intensity, min_peptide_intensity) ;
							if (!found_transform)
							{
								// this could be from a bad spectra
								// Look in header		
 								short cs = raw_data->GetMonoChargeFromHeader(scan_num) ; 
								if (cs > 0)
								{
									//create transform record with  instrument recorded charge cs and mass
									this_transform_record.mshort_cs = cs ; 									
									this_transform_record.mdbl_mono_mw = (mono_mz - mdbl_cc_mass) * cs ; 
									this_transform_record.mdbl_fit = 1; 
									if (mono_mz_match == 0) //couldnt find the parent, so set intensity to noise level
										this_transform_record.mint_abundance = (int) min_peptide_intensity ; 
									else
										this_transform_record.mint_abundance = (int) monoPeak.mdbl_intensity ; 
									vect_transform_records.push_back(this_transform_record) ; 
								}
								else
								{
									// instrument couldnt find charge as well, so create dummy transforms with cs 2 and 3, if user has agreed to do this
									if (scoring_parameters->get_AllocateDefaultCharges())
									{
										for( int charge = 2 ; charge <= 3 ; charge++)
										{
											this_transform_record.mshort_cs = charge ; 
											this_transform_record.mdbl_mono_mw = (parent_mz - mdbl_cc_mass) * charge ; 
											this_transform_record.mdbl_fit = 1; 
											if (parent_match == 0) //couldnt find the parent, so set intensity to noise level
												this_transform_record.mint_abundance = (int) min_peptide_intensity ; 
											else
												this_transform_record.mint_abundance = (int) parentPeak.mdbl_intensity  ; 
											vect_transform_records.push_back(this_transform_record) ; 
										}
									}
								}
							}
							else
							{
								// If summing is enabled, replace the detected mono intensity with parent intensity
								if(transform_parameters->get_SumSpectra())
								{
									for (int k = 0; k < (int)vect_transform_records.size(); k++)
									{		
										Engine::PeakProcessing::Peak tempPeak ;
										double mono_match = temp_parent_peak_processor->GetClosestPeakMz(vect_transform_records[k].mdbl_mz, tempPeak) ; 
										if (tempPeak.mdbl_intensity >  min_peptide_intensity)
										{
											vect_transform_records[k].mint_abundance = (int) tempPeak.mdbl_intensity ;
											vect_transform_records[k].mint_mono_intensity =(int)tempPeak.mdbl_intensity ;
										}
										else
										{
											// In case if precursor is in the noise floor. 
											vect_transform_records[k].mint_abundance = (int) min_peptide_intensity ; 
											vect_transform_records[k].mint_mono_intensity = (int) min_peptide_intensity ;
										}
									}
								}


							}

							// ------------------ CID Scoring --------------- //
							double msn_thres =  GlypID::Utils::GetAverage(msn_vect_intensities, FLT_MAX) ;
							background_intensity = 0 ; //GlypID::Utils::GetAverage(hcd_vect_intensities, (float)(5*thres)) ;
							cid_peak_processor->SetPeakIntensityThreshold(background_intensity * mobj_peak_parameters->get_PeakBackgroundRatio()) ; 							
							cid_peak_processor->SetPeaksProfileType(raw_data->IsProfileScan(scan_num)) ;
							int numCIDPeaks = cid_peak_processor->DiscoverPeaks(&msn_vect_mzs, &msn_vect_intensities) ; 
							cid_peak_processor->mobj_peak_data->InitializeUnprocessedPeakData() ;		
							// store in map
							parent_present.insert(parent_pair (parent_mz, 1)) ; 																
							for (int k = 0; k < (int)vect_transform_records.size(); k++)
							{
								this_transform_record = vect_transform_records[k] ; 
								// get score
								double cid_score_value = cid_scoring->CalculateCIDScore(*cid_peak_processor->mobj_peak_data, this_transform_record.mshort_cs ) ; 
								// get p value
								double p_value = cid_scoring->CalculateCIDScorePValue((int) cid_score_value) ; 

								// see if oxonium is present
								bool contains_oxonium = cid_scoring->CheckForOxoniumIons(*cid_peak_processor->mobj_peak_data) ; 
								//get mono mz
								double mono_mz = (this_transform_record.mdbl_mono_mw/this_transform_record.mshort_cs) + mdbl_cc_mass  ;

								this_cid_record.Clear() ; 
								this_cid_record.AddInfoToCIDRecord(scan_num, parent_scan, scan_ms_level, 1, parent_mz, mono_mz, 0, this_transform_record.mint_abundance, 
									this_transform_record.mshort_cs, this_transform_record.mdbl_fit, this_transform_record.mdbl_mono_mw, cid_score_value, p_value, contains_oxonium, 0 , 0, 0) ; 
								this_cid_record.AddPeaks(*cid_peak_processor->mobj_peak_data) ; 

								
								

								// Do glycopeptide searching 
								//but make sure all matches within tolerance are taken ito consideration
								std::vector<Engine::MS2CIDScoring::CIDInformationRecord> temp_cid_records ; 
								temp_cid_records.push_back(this_cid_record) ; 
								Engine::MS2HCDScoring::GLYCAN_TYPE type = Engine::MS2HCDScoring::GLYCAN_TYPE::NA ; 
								glycan_processor->SearchAndFillPeptideInformation(vect_glycan_compositions, vect_candidate_proteins, temp_cid_records, type, false, false) ; 

								//store into cid records
								vect_cid_score_records.insert(vect_cid_score_records.end(), temp_cid_records.begin(), temp_cid_records.end()); 
								temp_cid_records.clear() ; 

								
							}		

						}
					}					
				}
				else
				{
					parent_present.clear() ; 
				}
				scan_num++ ; 
				Console::WriteLine(scan_num) ;
				
			}

			//----------- Store and return ----------//			
			comb_scoring_results->AddMultipleCombinedCIDETDRecords(vect_comb_records) ; 
			all_results->SetCombinedScoringResults(comb_scoring_results) ; 
			
			mint_percent_done = 100 ; 
		}
		catch (char *mesg)
		{
			// blah if (comb_scoring_results != NULL)
			/*{
				delete comb_scoring_results ; 
				comb_scoring_results = NULL ;
			}			
		/*	if (parent_peak_processor != NULL) 
			{
				delete parent_peak_processor ;
				parent_peak_processor = NULL ; 
			}
			if (cid_peak_processor != NULL) 
			{
				delete cid_peak_processor ;
				cid_peak_processor = NULL ; 
			}
			if (etd_peak_processor != NULL) 
			{
				delete etd_peak_processor ;
				etd_peak_processor = NULL ; 
			}
			if( raw_data != NULL)
			{
				raw_data->Close() ; 
				delete raw_data ; 
			}
			if (mass_transform != NULL)
			{
				delete mass_transform ;
				mass_transform = NULL ; 
			}
			if(etd_scoring != NULL)
			{
				delete etd_scoring ; 
				etd_scoring = NULL ; 
			}
			if(cid_scoring != NULL)
			{
				delete cid_scoring ; 
				cid_scoring = NULL ; 
			}
			if (glycan_data != NULL)
			{
				delete glycan_data ; 
				glycan_data = NULL ; 
			}
			if (glycan_processor != NULL)
			{
				delete glycan_processor  ;
				glycan_processor = NULL ; 
			}
			System::String *exception_msg = new System::String(mesg) ; 
			menm_state = enmProcessState::ERROR ;
			throw new System::Exception(exception_msg) ; 
		}			
		if (parent_peak_processor != NULL) 
		{
			delete parent_peak_processor ;
			parent_peak_processor = NULL ; 
		}
		if (cid_peak_processor != NULL) 
		{
			delete cid_peak_processor ;
			cid_peak_processor = NULL ; 
		}
		if (etd_peak_processor != NULL) 
		{
			delete etd_peak_processor ;
			etd_peak_processor = NULL ; 
		}
		if( raw_data != NULL)
		{
			raw_data->Close() ; 
			delete raw_data ; 
		}
		if (mass_transform != NULL)
		{
			delete mass_transform ;
			mass_transform = NULL ; 
		}
		if(etd_scoring != NULL)
		{
			delete etd_scoring ; 
			etd_scoring = NULL ; 
		}
		if(cid_scoring != NULL)
		{
			delete cid_scoring ; 
			cid_scoring = NULL ; 
		}
		if (glycan_data != NULL)
		{
			delete glycan_data ; 
			glycan_data = NULL ; 
		}
		if (glycan_processor != NULL)
		{
			delete glycan_processor  ;
			glycan_processor = NULL ; 
		}

		menm_state = COMPLETE ;*/
		return all_results ;
	}




	void clsProcRunner::CreateTransformResults()
	{
		if (mstr_file_name == NULL)
		{
			throw new System::Exception(S"File name is not set.") ;
		}
		if (mobj_peak_parameters == NULL)
		{
			throw new System::Exception(S"Peak parameters not set.") ;
		}
		if (mobj_transform_parameters == NULL)
		{
			throw new System::Exception(S"Horn Transform parameters not set.") ;
		}
		mobj_transform_results = CreateTransformResults(mstr_file_name, menm_file_type,  mobj_peak_parameters, mobj_transform_parameters) ;	
	}	
	
	void clsProcRunner::CreateCIDScoringResults()
	{
		if (mstr_file_name == NULL)
		{
			throw new System::Exception(S"File name is not set.") ; 
		}
		if (mobj_peak_parameters == NULL)
		{
			throw new System::Exception(S"Peak parameters not set.") ;
		}
		if (mobj_transform_parameters == NULL)
		{
			throw new System::Exception(S"Horn Transform parameters not set.") ;
		}
		if(mobj_scoring_parameters == NULL)
		{
			throw new System::Exception(S"Scoring parameters not set.") ; 
		}
		if (mbln_search_fasta)
		{
			if (mstr_protein_list_name == NULL)
			{
				throw new System::Exception(S"Protein List name not set.") ; 
			}
			if (mstr_glycan_composition_list_name == NULL)
			{
				throw new System::Exception(S"Glycan Composition list name not set.") ; 
			}
		}
		mobj_cid_scoring_results = CreateCIDScoringResults(mstr_file_name, menm_file_type, mstr_protein_list_name, mstr_glycan_composition_list_name, mobj_peak_parameters, mobj_transform_parameters, mobj_scoring_parameters, mbln_search_fasta) ; 

	}

	void clsProcRunner::CreateHCDScoringResults()
	{
		if (mstr_file_name == NULL)
		{
			throw new System::Exception(S"File name is not set.") ;
		}
		if (mobj_peak_parameters == NULL)
		{
			throw new System::Exception(S"Peak parameters not set.") ;
		}
		if (mobj_transform_parameters == NULL)
		{
			throw new System::Exception(S"Horn Transform parameters not set.") ;
		}
		if(mobj_scoring_parameters == NULL)
		{
			throw new System::Exception(S"Scoring parameters not set.") ; 
		}
		if (mbln_search_fasta)
		{
			if (mstr_protein_list_name == NULL)
			{
				throw new System::Exception(S"Protein List name not set.") ; 
			}
			if (mstr_glycan_composition_list_name == NULL)
			{
				throw new System::Exception(S"Glycan Composition list name not set.") ; 
			}
		}
		mobj_hcd_scoring_results = CreateHCDScoringResults(mstr_file_name, menm_file_type, mstr_protein_list_name, mstr_glycan_composition_list_name, mobj_peak_parameters, mobj_transform_parameters, mobj_scoring_parameters, mbln_search_fasta) ; 
	}

	void clsProcRunner::CreateCIDHCDScoringResults()
	{
		if (mstr_file_name == NULL)
		{
			throw new System::Exception(S"File name is not set.") ;
		}
		if (mobj_peak_parameters == NULL)
		{
			throw new System::Exception(S"Peak parameters not set.") ;
		}
		if (mbln_search_fasta)
		{
			if (mstr_protein_list_name == NULL)
			{
				throw new System::Exception(S"Protein List name not set.") ; 
			}
			if (mstr_glycan_composition_list_name == NULL)
			{
				throw new System::Exception(S"Glycan Composition list name not set.") ; 
			}
		}
		if (mobj_transform_parameters == NULL)
		{
			throw new System::Exception(S"Horn Transform parameters not set.") ;
		}
		if(mobj_scoring_parameters == NULL)
		{
			throw new System::Exception(S"Scoring parameters not set.") ; 
		}
		mobj_cid_hcd_combined_scoring_results = CreateCIDHCDScoringResults(mstr_file_name, menm_file_type, mstr_protein_list_name, mstr_glycan_composition_list_name, mobj_peak_parameters, mobj_transform_parameters, mobj_scoring_parameters, mbln_search_fasta) ; 
	}	

	void clsProcRunner::CreateCIDETDScoringResults()
	{

		if (mstr_file_name == NULL)
		{
			throw new System::Exception(S"File name is not set.") ;
		}
		if (mobj_peak_parameters == NULL)
		{
			throw new System::Exception(S"Peak parameters not set.") ;
		}
		if (mbln_search_fasta)
		{
			if (mstr_protein_list_name == NULL)
			{
				throw new System::Exception(S"Protein List name not set.") ; 
			}
			if (mstr_glycan_composition_list_name == NULL)
			{
				throw new System::Exception(S"Glycan Composition list name not set.") ; 
			}
		}
		if (mobj_transform_parameters == NULL)
		{
			throw new System::Exception(S"Horn Transform parameters not set.") ;
		}
		if(mobj_scoring_parameters == NULL)
		{
			throw new System::Exception(S"Scoring parameters not set.") ; 
		}

		mobj_cid_etd_combined_scoring_results = CreateCIDETDScoringResults(mstr_file_name, menm_file_type, mstr_protein_list_name, mstr_glycan_composition_list_name, mobj_peak_parameters, mobj_transform_parameters, mobj_scoring_parameters, mbln_search_fasta) ; 

	}
}