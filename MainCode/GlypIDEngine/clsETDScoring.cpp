#include "clsETDScoring.h"
#include "clsPeakProcessor.h"
#include "clshorntransform.h"

namespace GlypID
{
	namespace ETDScoring
	{	
			
		clsETDScoring::clsETDScoring(void)
		{
			try
			{
				mobj_etd_scoring = new Engine::MS2ETDScoring::ETDScoring() ;  
				mobj_mwt =  __gc new MwtWinDll::MolecularWeightCalculator()  ;
			

				
			}
			catch (char *mesg)
			{
				System::String *exception_msg = new System::String(mesg) ; 
				throw new System::Exception(exception_msg) ; 
			}

		}

		clsETDScoring::~clsETDScoring(void)
		{
		}

		void clsETDScoring::SetOptions(short max_charge, double max_mw, double max_fit, double min_s2n, double cc_mass, 
			double delete_threshold_theoretical_intensity, double min_theoretical_intensity_for_score, short num_peaks_for_shoulder, 
			bool use_caching, bool o16_o18_media, bool check_against_charge_1)
		{
			mobj_transform->SetOptions(max_charge, max_mw, max_fit, min_s2n, cc_mass, delete_threshold_theoretical_intensity,
				min_theoretical_intensity_for_score, num_peaks_for_shoulder, check_against_charge_1, use_caching, o16_o18_media) ; 
		}

		void clsETDScoring::SetIsotopeFitOptions(System::String *str_averagine, System::String *str_tag, bool thrash_or_not, 
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

		void clsETDScoring::SetScanDetails(int msNScan, int parent_scan, int msNLevel, int parentLevel, double parentMz )
		{
			mint_msn_scan = msNScan ; 
			mint_parent_scan = parent_scan ; 
			mint_msn_level = msNLevel; 
			mint_parent_level = parentLevel ; 
			mdbl_parent_mz = parentMz ; 
		}

		void clsETDScoring::PreprocessETDSpectra(GlypID::Peaks::clsPeak* (etd_peaks) __gc[],  Engine::PeakProcessing::PeakData& etd_proc_peakData, double parent_mz) 
		{
			if (etd_peaks->Length == 0)
				return ; 


			Engine::PeakProcessing::PeakData peakData ; 	
			Engine::PeakProcessing::PeakData tpeakData ; 
			GlypID::Utils::SetPeaks(tpeakData, etd_peaks) ; 
			etd_proc_peakData.Clear() ; 

			// -------- Preprocessing ------//

			// Method 1 : top n peaks
			// AM L Nov 11 : change to top x% intensity
			/*int num_peaks_to_select = 100 + 4 ; //peakData.GetNumPeaks() * 0.1  + 4;// just to   account for the parent peak and shoulders

			for (int i=0 ; i < num_peaks_to_select ; i++) // Need to parameterize this
			{
				Engine::PeakProcessing::Peak thisPeak ; 
				bool isPeak = peakData.GetNextPeak(0, parent_mz, thisPeak) ; 
				
				// AM Nov 2011 : change to remove precursor from list to be counted. 
				bool not_parent = true ; 
				if (abs(thisPeak.mdbl_mz -parent_mz) < 1)
					not_parent = false ; 

				if (isPeak && not_parent)
				{
					etd_proc_peakData.AddPeak(thisPeak) ; 
				}
			}*/


			// Method 2 : top n in m peaks
			int num_peaks_per_bin = 1 ; 
			double bin_size = 5 ; //m/z

			double start_mz = 0;
			double stop_mz = start_mz + bin_size; 

			while(stop_mz < parent_mz - bin_size)
			{
				for (int i= 0; i < num_peaks_per_bin ; i++)
				{
					Engine::PeakProcessing::Peak thisPeak ; 
					bool isPeak = tpeakData.GetNextPeak(start_mz, stop_mz, thisPeak) ; 
					if (isPeak)
						peakData.AddPeak(thisPeak) ; 
					else 
						break ;
				}
				start_mz = stop_mz ; 
				stop_mz += bin_size ; 
			}
			peakData.InitializeUnprocessedPeakData() ; 

			int num_peaks_to_select = 20 ; 

			for (int i=0 ; i < num_peaks_to_select ; i++) // Need to parameterize this
			{
				Engine::PeakProcessing::Peak thisPeak ; 
				bool isPeak = peakData.GetNextPeak(0, parent_mz, thisPeak) ; 
				
				// AM Nov 2011 : change to remove precursor from list to be counted. 
				bool not_parent = true ; 
				if (abs(thisPeak.mdbl_mz -parent_mz) < 1)
					not_parent = false ; 

				if (isPeak && not_parent)
				{
					etd_proc_peakData.AddPeak(thisPeak) ; 
				}
			}
		}


		double clsETDScoring::ScoreETDSpectraOLinked(GlypID::Peaks::clsPeak* (&etd_peaks) __gc[],  GlypID::ETDScoring::clsETDScoringScanResults* etd_scoring_result) 
		{
			std::vector<double> vectMzs ;
			std::vector<double> vectIntensities ;
			double etd_score  = 0 ; 
			MwtWinDll::MWPeptideClass::udtFragmentationSpectrumOptionsType *udtFragSpectrumOptions ; 	
			udtFragSpectrumOptions = &mobj_mwt->Peptide->GetFragmentationSpectrumOptions();
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itAIon].ShowIon = false;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itBIon].ShowIon = true;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itBIon].NeutralLossWater = true ; 
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itYIon].ShowIon = true;			
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itYIon].NeutralLossWater = true ; 
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itCIon].ShowIon = true;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itZIon].ShowIon = true;
			udtFragSpectrumOptions->IntensityOptions.IonType[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itZIon] = 100 ; 
			udtFragSpectrumOptions->IntensityOptions.IonType[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itCIon] = 100 ; 
			udtFragSpectrumOptions->IntensityOptions.IonType[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itBIon] = 50 ; 
			udtFragSpectrumOptions->IntensityOptions.IonType[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itYIon] = 50 ;
			udtFragSpectrumOptions->IntensityOptions.NeutralLoss = 25 ; 
			udtFragSpectrumOptions->IntensityOptions.BYIonShoulder = 0 ; 


			if (etd_peaks->Length == 0)
				return false  ; 

			std::vector <double> obs_peak_mzs ; 
			std::vector <double> obs_peak_intensities ;
			Engine::PeakProcessing::PeakData peakData ; 	

			// ---- Preprocessing --- //
			PreprocessETDSpectra(etd_peaks, peakData, etd_scoring_result->mdbl_parent_mz) ; 		


			Engine::PeakProcessing::Peak thisPeak ; 
			for (int i= 0; i < peakData.GetNumPeaks() ; i++)
			{
				peakData.GetPeak(i, thisPeak) ; 
				obs_peak_mzs.push_back(thisPeak.mdbl_mz) ; 
				obs_peak_intensities.push_back(thisPeak.mdbl_intensity) ;
			}			
			

			Engine::MS2ETDScoring::ETDInformationRecord etdRecord ; 			
			double max_etd_score = 0 ; 
			std::string maxSite = "" ; 

			double glycan_mass = etd_scoring_result->mdbl_glycan_mass; 
			if (glycan_mass > 0)
			{
				switch(etd_scoring_result->mshort_cs)
				{
					case 1:  
						udtFragSpectrumOptions->DoubleChargeIonsShow = false;
						udtFragSpectrumOptions->TripleChargeIonsShow = false;					 
						break ; 
					case 2:
						udtFragSpectrumOptions->DoubleChargeIonsShow = true;
						udtFragSpectrumOptions->TripleChargeIonsShow = false;
						udtFragSpectrumOptions->DoubleChargeIonsThreshold = 0 ; 
						break ;
					default:
						udtFragSpectrumOptions->DoubleChargeIonsShow = true;
						udtFragSpectrumOptions->TripleChargeIonsShow = true;
						udtFragSpectrumOptions->DoubleChargeIonsThreshold = 0 ;
						udtFragSpectrumOptions->TripleChargeIonsThreshold = 0 ; 
						break ;
				}

				mobj_mwt->Peptide->SetFragmentationSpectrumOptions(udtFragSpectrumOptions) ; 

				// Add glycan modification to list of modifications in MWT
			
				System::String *symbol = "^" ; 
				char chr_symbol [20] ; 
				GlypID::Utils::GetStr(symbol, chr_symbol);
				System::String *comment = "Glycan";
				bool var = false ; 
				double mdbl_H2O = 18.01056; 
				double mod_mass = glycan_mass - mdbl_H2O ; 
				mobj_mwt->Peptide->SetModificationSymbol(symbol, mod_mass , &var, &comment) ; 	

				// Processing each Oglycosite separetely
				Engine::GlycanCompositionManager::GlycanProcessor gp_processor ;
				char peptide[512]; 	
				peptide[0] = '\0'; 
				GlypID::Utils::GetStr(etd_scoring_result->mstr_pep_seq, peptide);
				int nextO = etd_scoring_result->mstr_oglyco_site->IndexOf("/") ; 
				int prevO = 0 ; 
				while (prevO < etd_scoring_result->mstr_oglyco_site->Length  )
				{	
					if (nextO == -1)
					{
						nextO = etd_scoring_result->mstr_oglyco_site->Length ;
					}
					
					//  Add glycan modification to the sequence
					char site [512] ; 
					site[0] = '\0' ; 
					char s_or_t ; 
					GlypID::Utils::GetStr(etd_scoring_result->mstr_oglyco_site, site) ;
					char extract_site[20]; 
					std::string tempSite ; 
					strncpy(extract_site, site+prevO+1, nextO-prevO-1) ; 	
					extract_site[nextO-prevO-1] = '\0' ; 
					s_or_t = site[prevO] ; 
					std::string sequencewithGlyc = "";
					int pos = atoi(extract_site) ; 
					mobj_etd_scoring->AddOGlycanModificationToSequence(chr_symbol, peptide, sequencewithGlyc, pos) ; 
					std::string sequence = "" ; 
					mobj_etd_scoring->AddCarboModificationToSequence(sequencewithGlyc, sequence); // Carbo modification

					// Get theoretical spectra					
					mobj_mwt->Peptide->SetSequence(sequence.c_str(), 
						MwtWinDll::MWPeptideClass::ntgNTerminusGroupConstants::ntgHydrogen,
						MwtWinDll::MWPeptideClass::ctgCTerminusGroupConstants::ctgHydroxyl,
						false);		

					double temp_m = mobj_mwt->Peptide->GetPeptideMass() ; 

					MwtWinDll::MWPeptideClass::udtFragmentationSpectrumDataType udtFrag[]; 
					mobj_mwt->Peptide->GetFragmentationMasses(&udtFrag) ; 
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
					etd_score = mobj_etd_scoring->CalculateETDScore2(obs_peak_mzs, obs_peak_intensities, theor_mzs, theor_intensities) ; 

					if (etd_score > max_etd_score)
					{
						max_etd_score = etd_score ; 
						etd_scoring_result->mstr_pep_seq = sequencewithGlyc.c_str() ; 
						char chr_pos[10];
						itoa(pos, chr_pos, 10) ; 
						maxSite = strcat(&s_or_t, chr_pos) ; 
					}

					// Get next site
					prevO = nextO+1;
					if (prevO < etd_scoring_result->mstr_oglyco_site->Length)
						nextO = etd_scoring_result->mstr_oglyco_site->IndexOf("/", prevO+1) ; 
				}
			}
			etd_scoring_result->mstr_oglyco_site = maxSite.c_str() ; 
			return max_etd_score ; 
		}
		
			
		double clsETDScoring::ScoreETDSpectra(GlypID::Peaks::clsPeak* (&etd_peaks) __gc[],  GlypID::ETDScoring::clsETDScoringScanResults* etd_scoring_result) 
		{
			
			std::vector<double> vectMzs ;
			std::vector<double> vectIntensities ;
			double etd_score  = 0 ; 
			MwtWinDll::MWPeptideClass::udtFragmentationSpectrumOptionsType *udtFragSpectrumOptions ; 	

		
	

			//udtFragSpectrumOptions->Initialize();
            udtFragSpectrumOptions = &mobj_mwt->Peptide->GetFragmentationSpectrumOptions();
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itAIon].ShowIon = false;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itBIon].ShowIon = true;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itBIon].NeutralLossWater = true ; 
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itYIon].ShowIon = true;			
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itYIon].NeutralLossWater = true ; 
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itCIon].ShowIon = true;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itZIon].ShowIon = true;
			udtFragSpectrumOptions->IntensityOptions.IonType[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itZIon] = 100 ; 
			udtFragSpectrumOptions->IntensityOptions.IonType[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itCIon] = 100 ; 
			udtFragSpectrumOptions->IntensityOptions.IonType[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itBIon] = 50 ; 
			udtFragSpectrumOptions->IntensityOptions.IonType[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itYIon] = 50 ;
			udtFragSpectrumOptions->IntensityOptions.NeutralLoss = 25 ; 
			udtFragSpectrumOptions->IntensityOptions.BYIonShoulder = 0 ; 
			
			//MwudtFragSpectrumOptions->IonTypeOpti
			
			

			if (etd_peaks->Length == 0)
				return false  ; 

			std::vector <double> obs_peak_mzs ; 
			std::vector <double> obs_peak_intensities ;
			Engine::PeakProcessing::PeakData peakData ; 	

			// ---- Preprocessing --- //
			PreprocessETDSpectra(etd_peaks, peakData, etd_scoring_result->mdbl_parent_mz) ; 		


			Engine::PeakProcessing::Peak thisPeak ; 
			for (int i= 0; i < peakData.GetNumPeaks() ; i++)
			{
				peakData.GetPeak(i, thisPeak) ; 
				obs_peak_mzs.push_back(thisPeak.mdbl_mz) ; 
				obs_peak_intensities.push_back(thisPeak.mdbl_intensity) ;
			}
			
			

			Engine::MS2ETDScoring::ETDInformationRecord etdRecord ; 
			
			double max_etd_score = 0 ; 
			double glycan_mass = etd_scoring_result->mdbl_glycan_mass; 
			if (glycan_mass > 0)
			{
				switch(etd_scoring_result->mshort_cs)
				{
					case 1:  
						udtFragSpectrumOptions->DoubleChargeIonsShow = false;
						udtFragSpectrumOptions->TripleChargeIonsShow = false;					 
						break ; 
					case 2:
						udtFragSpectrumOptions->DoubleChargeIonsShow = true;
						udtFragSpectrumOptions->TripleChargeIonsShow = false;
						udtFragSpectrumOptions->DoubleChargeIonsThreshold = 0 ; 
						break ;
					default:
						udtFragSpectrumOptions->DoubleChargeIonsShow = true;
						udtFragSpectrumOptions->TripleChargeIonsShow = true;
						udtFragSpectrumOptions->DoubleChargeIonsThreshold = 0 ;
						udtFragSpectrumOptions->TripleChargeIonsThreshold = 0 ; 
						break ;
				}

				mobj_mwt->Peptide->SetFragmentationSpectrumOptions(udtFragSpectrumOptions) ; 

				// Add glycan modification to list of modifications in MWT
			
				System::String *symbol = "^" ; 
				char chr_symbol [20] ; 
				GlypID::Utils::GetStr(symbol, chr_symbol);
				System::String *comment = "Glycan";
				bool var = false ; 
				double mdbl_H2O = 18.01056; 
				double mod_mass = glycan_mass - mdbl_H2O ; 
				mobj_mwt->Peptide->SetModificationSymbol(symbol, mod_mass , &var, &comment) ; 	

				
				// July 2012, Doing multiple glycosylation sites seperately.
				Engine::GlycanCompositionManager::GlycanProcessor gp_processor ;
				
				System::String *site_Positions = etd_scoring_result->mstr_glyco_site; 
				char peptide[512]; 	
				peptide[0] = '\0'; 
				GlypID::Utils::GetStr(etd_scoring_result->mstr_pep_seq, peptide); 
				int nextN = gp_processor.GetMotifPositionWithinSequence(peptide, 0);
				int positionN=0 ; 				
				while( nextN != -1)
				{
					// Indicate modifications in sequence
					std::string sequencewithGlyc = "" ; 					
					int pos = mobj_etd_scoring->AddGlycanModificationToSequence(chr_symbol, peptide , sequencewithGlyc, nextN) ; 		// Glycan modification					
					std::string sequence = "" ; 
					mobj_etd_scoring->AddCarboModificationToSequence(sequencewithGlyc, sequence); // Carbo modification

						// Get theoretical spectra					
					mobj_mwt->Peptide->SetSequence(sequence.c_str(), 
						MwtWinDll::MWPeptideClass::ntgNTerminusGroupConstants::ntgHydrogen,
						MwtWinDll::MWPeptideClass::ctgCTerminusGroupConstants::ctgHydroxyl,
						false);		

					double temp_m = mobj_mwt->Peptide->GetPeptideMass() ; 

					MwtWinDll::MWPeptideClass::udtFragmentationSpectrumDataType udtFrag[]; 
					mobj_mwt->Peptide->GetFragmentationMasses(&udtFrag) ; 
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
					//etd_score = mobj_etd_scoring->CalculateETDScore(obs_peak_mzs , obs_peak_intensities, theor_mzs, theor_intensities) ; 										
					etd_score = mobj_etd_scoring->CalculateETDScore2(obs_peak_mzs, obs_peak_intensities, theor_mzs, theor_intensities) ; 

					// Normalizing
					double sum_score = 0; 
					/*for(int shift = -15; shift < 16; shift++)
					{
						std::vector<double> shiftmz; 
						for (int jj = 0 ; jj < obs_peak_mzs.size() ; jj++)
							shiftmz.push_back(obs_peak_mzs.at(jj) + shift) ; 
						double this_score = mobj_etd_scoring->CalculateETDScore2(shiftmz, obs_peak_intensities, theor_mzs, theor_intensities) ; 
						sum_score  += this_score ; 

					}*/

					//etd_score = etd_score/sum_score ; 

					if (etd_score > max_etd_score)
					{
						max_etd_score = etd_score ; 
						etd_scoring_result->mstr_pep_seq = sequencewithGlyc.c_str() ; 
					}

					nextN = gp_processor.GetMotifPositionWithinSequence(peptide, nextN+1) ; 
				}
			}	
			return max_etd_score; 
		}	



		double clsETDScoring::ScoreMS3CIDSpectra(GlypID::Peaks::clsPeak* (&etd_peaks) __gc[],  GlypID::ETDScoring::clsETDScoringScanResults* etd_scoring_result) 
		{
			
			std::vector<double> vectMzs ;
			std::vector<double> vectIntensities ;
			double etd_score  = 0 ; 
			MwtWinDll::MWPeptideClass::udtFragmentationSpectrumOptionsType *udtFragSpectrumOptions ; 
			

		/// ANOOP Change May 5, 2012. to score b and y ms3. 


			//udtFragSpectrumOptions->Initialize();
            udtFragSpectrumOptions = &mobj_mwt->Peptide->GetFragmentationSpectrumOptions();
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itAIon].ShowIon = false;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itBIon].ShowIon = true;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itYIon].ShowIon = true;			
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itCIon].ShowIon = false;
			udtFragSpectrumOptions->IonTypeOptions[(int)MwtWinDll::MWPeptideClass::itIonTypeConstants::itZIon].ShowIon = false;

			if (etd_peaks->Length == 0)
				return false  ; 

			Engine::PeakProcessing::PeakData peakData ; 			
			GlypID::Utils::SetPeaks(peakData, etd_peaks) ; 

			std::vector <double> obs_peak_mzs ; 
			std::vector <double> obs_peak_intensities ;
			// AM L Nov 11 : change to top x% intensity
			int num_peaks_to_select = 100 + 4 ; //peakData.GetNumPeaks() * 0.1  + 4;// just to   account for the parent peak and shoulders

			for (int i=0 ; i < num_peaks_to_select ; i++) // Need to parameterize this
			{
				Engine::PeakProcessing::Peak thisPeak ; 
				bool isPeak = peakData.GetNextPeak(300, 2000, thisPeak) ; 
				
				// AM Nov 2011 : change to remove precursor from list to be counted. 
				bool not_parent = true ; 
				if (abs(thisPeak.mdbl_mz - etd_scoring_result->mdbl_parent_mz) < 1)
					not_parent = false ; 

				if (isPeak && not_parent)
				{
					obs_peak_mzs.push_back(thisPeak.mdbl_mz) ; 
					obs_peak_intensities.push_back(thisPeak.mdbl_intensity) ; 
				}
			}

			Engine::MS2ETDScoring::ETDInformationRecord etdRecord ; 
			
			double max_score = 0 ; 
			double glycan_mass = etd_scoring_result->mdbl_glycan_mass; 
			if (glycan_mass > 0)
			{
				switch(etd_scoring_result->mshort_cs)
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
				double mdbl_H2O = 18.01056; 
				double mod_mass = glycan_mass - mdbl_H2O ; 
				mobj_mwt->Peptide->SetModificationSymbol(symbol, mod_mass , &var, &comment) ; 					

				// Indicate modification in sequence
				std::string sequencewithGlyc = "" ; 
				char peptide[512]; 
				//peptide[0] = '\0';
				GlypID::Utils::GetStr(etd_scoring_result->mstr_pep_seq, peptide); 
				int pos = mobj_etd_scoring->AddGlycanModificationToSequence(chr_symbol, peptide , sequencewithGlyc, 0 ) ; 		
				std::string sequence = "" ; 
				mobj_etd_scoring->AddCarboModificationToSequence(sequencewithGlyc, sequence); 

				// Get theoretical spectra					
				mobj_mwt->Peptide->SetSequence(sequence.c_str(), 
					MwtWinDll::MWPeptideClass::ntgNTerminusGroupConstants::ntgHydrogen,
					MwtWinDll::MWPeptideClass::ctgCTerminusGroupConstants::ctgHydroxyl,
					false);		
				MwtWinDll::MWPeptideClass::udtFragmentationSpectrumDataType udtFrag[]; 
				mobj_mwt->Peptide->GetFragmentationMasses(&udtFrag) ; 
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
				//etd_score = mobj_etd_scoring->CalculateETDScore(obs_peak_mzs , obs_peak_intensities, theor_mzs, theor_intensities) ; 										
				etd_score = mobj_etd_scoring->CalculateETDScore2(obs_peak_mzs, obs_peak_intensities, theor_mzs, theor_intensities) ; 
			}	
			return etd_score; 
		}	


	}
}