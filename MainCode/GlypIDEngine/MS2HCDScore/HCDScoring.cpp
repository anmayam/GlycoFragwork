// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University

#include "./hcdscoring.h"

#include "../Utilities/system.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctype.h>
#include <io.h>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <fcntl.h>

namespace Engine
{
	namespace MS2HCDScoring
	{		

		HCDScoring::HCDScoring(void)
		{
			mdbl_max_mz = 700 ;
			mdbl_min_mz = 0 ; 
			mint_min_number_peaks  = 0 ; 
			marr_theoretical_peaks[6] = 657.23 ; // NeuAc + Gal/Mannose + GlcNac
			marr_theoretical_peaks[5] = 366.14; // GlcNAc + Gal/Mannose
			marr_theoretical_peaks[4] = 292.1; //NeuAc
			marr_theoretical_peaks[3] = 274.09; // NeuAc-water loss
			marr_theoretical_peaks[2] = 204.09; // GlcNAc
			marr_theoretical_peaks[1] = 163.06 ; // Gal/Mannose //prev 168.06
			marr_theoretical_peaks[0] = 138.05 ; // Oxonium ion at HexNac

			//Distribution of intensities for hm
			marr_theoretical_hm[0] = 0 ; 
			marr_theoretical_hm[1] = 0.99 ; 
			marr_theoretical_hm[2]= 0 ;
			marr_theoretical_hm[3] = 0;
			marr_theoretical_hm[4] = 0;
			marr_theoretical_hm[5] = 0.3;
			marr_theoretical_hm[6] = 0 ; 

			//Distribution of intensities for ca
			marr_theoretical_ca[0] = 0 ; 
			marr_theoretical_ca[1] = 0.3 ; 
			marr_theoretical_ca[2]= 0.5 ;
			marr_theoretical_ca[3] = 0;
			marr_theoretical_ca[4] = 0;
			marr_theoretical_ca[5] = 0.99;
			marr_theoretical_ca[6] = 0 ; 
			
			//Distribution of intensities for cs
			marr_theoretical_cs[0] = 0.7 ; 
			marr_theoretical_cs[1] = 0.7 ; 
			marr_theoretical_cs[2]= 0.7;
			marr_theoretical_cs[3] = 0.7;
			marr_theoretical_cs[4] = 0.8;
			marr_theoretical_cs[5] = 0.7;
			marr_theoretical_cs[6] = 0.7 ; 

			//Distribution of intensities for hybrid
			marr_theoretical_hybrid[0] = 0.6 ; 
			marr_theoretical_hybrid[1] = 0.99 ; 
			marr_theoretical_hybrid[2]= 0.6;
			marr_theoretical_hybrid[3] = 0.6;
			marr_theoretical_hybrid[4] = 0.6;
			marr_theoretical_hybrid[5] = 0.6;
			marr_theoretical_hybrid[6] = 0.6 ; 

			// CC Mass
			mdbl_cc_mass = 1.00727638;		    

		}

		HCDScoring::~HCDScoring(void)
		{
		}

		void HCDScoring::GetOptions(int &min_num_peaks, double &min_mz, double &max_mz)
		{
			max_mz = mdbl_max_mz ; 
			min_mz = mdbl_min_mz ; 
			min_num_peaks  = mint_min_number_peaks ; 
		}

		void HCDScoring::SetOptions(int min_num_peaks, double min_mz, double max_mz)
		{
			mint_min_number_peaks = min_num_peaks ; 
			mdbl_min_mz = min_mz ; 
			mdbl_max_mz = max_mz ; 
		}

		
		double HCDScoring::CalculateHCDScore(Engine::PeakProcessing::PeakData &pk_data, std::vector <int> &identified_peak_indices)
		{
			// simple scoring scheme
			int num_hcd_present = 0 ; 
			double hcd_score = 0 ; 

			std::vector<double> *mzs = pk_data.mptr_vect_mzs ; 
			std::vector<double> *intensities = pk_data.mptr_vect_intensities ; 

			Engine::PeakProcessing::Peak currentPeak ; 
			mint_peak_indices.clear() ; 

			bool found_peak = pk_data.GetNextPeak(mdbl_min_mz, mdbl_max_mz, currentPeak) ; 
			while(found_peak)
			{
				for (int i= 0; i < 7; i++)
				{
					double this_mz = marr_theoretical_peaks[i] ; 
					if(abs(currentPeak.mdbl_mz-this_mz) < 0.1)
					{
						num_hcd_present++ ; 
						mint_peak_indices.push_back(currentPeak.mint_peak_index); 
					}

				}
				pk_data.RemovePeak(currentPeak) ; 
				found_peak = pk_data.GetNextPeak(mdbl_min_mz, mdbl_max_mz, currentPeak) ;
				
			}
			if (num_hcd_present > mint_min_number_peaks)
			{
				for (int i=0 ; i < (int) mint_peak_indices.size() ; i++)
				{
					identified_peak_indices.push_back(mint_peak_indices[i]) ; 
				}
			}
			hcd_score = (double)num_hcd_present/(double)6 ; 
			if (hcd_score >1)
				hcd_score = 1; 

			return hcd_score ; 
		}

		double HCDScoring::CalculateHCDScorePValue(Engine::PeakProcessing::PeakData &pk_data, std::vector <int> &identified_peak_indices)
		{
			int num_hcd_present = 0 ; 
			double hcd_score = 0 ; 
			double del_mz_bin = 0.1 ; 
			double del_mz_range = mdbl_max_mz - mdbl_min_mz ; 

			double p = (2*del_mz_bin)/del_mz_range ; 
			double p_value = 0 ; 

			std::vector<double> *mzs = pk_data.mptr_vect_mzs ; 
			std::vector<double> *intensities = pk_data.mptr_vect_intensities ; 

			Engine::PeakProcessing::Peak currentPeak ; 
			mint_peak_indices.clear() ; 

			bool found_peak = pk_data.GetNextPeak(mdbl_min_mz, mdbl_max_mz, currentPeak) ; 
			while(found_peak)
			{
				for (int i= 0; i < 7; i++)
				{
					double this_mz = marr_theoretical_peaks[i] ; 
					if(abs(currentPeak.mdbl_mz-this_mz) < del_mz_bin)
					{
						num_hcd_present++ ; 
						mint_peak_indices.push_back(currentPeak.mint_peak_index); 
					}

				}
				pk_data.RemovePeak(currentPeak) ; 
				found_peak = pk_data.GetNextPeak(mdbl_min_mz, mdbl_max_mz, currentPeak) ;
				
			}

			// Store all peaks away
			if (num_hcd_present > mint_min_number_peaks)
			{
				for (int i=0 ; i < (int) mint_peak_indices.size() ; i++)
				{
					identified_peak_indices.push_back(mint_peak_indices[i]) ; 
				}
			

				// Start calculating p value			
				long N_factorial = Engine::Utilities::factorial(7) ; 

				for (int x = num_hcd_present;  x <= 7;  x++)
				{
					long x_factorial = Engine::Utilities::factorial((long) x) ; 
					long N_x_factorial = Engine::Utilities::factorial((long) 7 - x) ; 
					double N_choose_x = (double) (N_factorial/(x_factorial *N_x_factorial)) ; 
					double pX = N_choose_x * pow(p, x) * pow(1-p,7-x) ; 
					p_value += pX ; 
				}	
			}
			else
				p_value = 0 ;  

			return p_value ; 
		}		

		double HCDScoring::DetermineGlycanType(Engine::PeakProcessing::PeakData &pk_data, std::vector <int> &identified_peak_indices, Engine::MS2HCDScoring::GLYCAN_TYPE &glycan_type) 
		{
			// Function that attempts to categorize the glycan type by calculating cross correlation

			//First find peaks and get p value
			int num_hcd_present = 0 ; 
			double hcd_score = 0 ; 
			double del_mz_bin = 0.015 ; 
			double del_mz_range = mdbl_max_mz - mdbl_min_mz ; 
			double p = (2*del_mz_bin)/del_mz_range ; 
			double p_value = 0 ; 
			std::vector<double> *mzs = pk_data.mptr_vect_mzs ; 
			std::vector<double> *intensities = pk_data.mptr_vect_intensities ; 
			Engine::PeakProcessing::Peak currentPeak ;
			float obs_intensities[7] ; 
			double maxintensity = 0  ; 
			for (int j = 0 ; j < 7 ; j++)
				obs_intensities[j] = 0 ; 		

			//bool found_peak = pk_data.GetNextPeak(mdbl_min_mz, mdbl_max_mz, currentPeak) ; 
			/*while(found_peak)
			{
				for (int i= 0; i < 7; i++)
				{
					double this_mz = marr_theoretical_peaks[i] ; 
					if(abs(currentPeak.mdbl_mz-this_mz) < del_mz_bin)
					{
						num_hcd_present++ ; 
						mint_peak_indices.push_back(currentPeak.mint_peak_index); 
						obs_intensities[i] = (float)currentPeak.mdbl_intensity ; 
						if(currentPeak.mdbl_intensity > maxintensity)
							maxintensity = currentPeak.mdbl_intensity ; 
						
					}
				}
				pk_data.RemovePeak(currentPeak) ; 
				found_peak = pk_data.GetNextPeak(mdbl_min_mz, mdbl_max_mz, currentPeak) ;				
			}*/


			mint_peak_indices.clear() ; 
			for(int i = 0; i < pk_data.GetNumPeaks() ; i++)
			{
				pk_data.GetPeak(i, currentPeak) ; 
				for (int j= 0; j < 7; j++)
				{
					double this_mz = marr_theoretical_peaks[j] ; 
					if(abs(currentPeak.mdbl_mz-this_mz) < del_mz_bin)
					{						
						if (currentPeak.mdbl_intensity > obs_intensities[j]) //Anoop jan 2011, make sure you get teh highest intensity
						{
							num_hcd_present++ ; 
							mint_peak_indices.push_back(currentPeak.mint_peak_index); 						
							obs_intensities[j] = (float)currentPeak.mdbl_intensity ; 
							if(currentPeak.mdbl_intensity > maxintensity)
								maxintensity = currentPeak.mdbl_intensity ; 
						}
						
					}
				}
			}

			// Normalizing
			for (int j = 0 ; j < 7 ; j++)
				obs_intensities[j] = (float) (obs_intensities[j]/maxintensity); 

			// Store all peaks away
			if (num_hcd_present >=mint_min_number_peaks)
			{
				for (int i=0 ; i < (int) mint_peak_indices.size() ; i++)
				{
					identified_peak_indices.push_back(mint_peak_indices[i]) ; 
				}	

				// Start calculating p value			
				long N_factorial = Engine::Utilities::factorial(7) ; 

				if (num_hcd_present > 7)
					num_hcd_present = 7 ;

				for (int x = num_hcd_present;  x <= 7;  x++)
				{
					long x_factorial = Engine::Utilities::factorial((long) x) ; 
					long N_x_factorial = Engine::Utilities::factorial((long) 7 - x) ; 
					double N_choose_x = (double) (N_factorial/(x_factorial *N_x_factorial)) ; 
					double pX = N_choose_x * pow(p, x) * pow(1-p,7-x) ; 
					p_value += pX ; 
				}	
				
				//Now do a correlation with each class
				float corr_HM = 0 ; 
				float corr_CA = 0 ; 
				float corr_CS = 0 ; 
				float corr_Hybrid =0  ; 
				float max_corr = 0  ; 

				//Anoop Jan 2011, performing the glycan type detection in two stages, 
				//if NeuAc peak  (or its water equivalent) is present, a distinction between hybrid and sialylated structures is made
				// else only look for high mannose and complex_asialylated
				if ((obs_intensities[4] > 0) || (obs_intensities[3]>0))
				{
					corr_CS = Engine::Utilities::correlation(&obs_intensities[0], &marr_theoretical_cs[0], 0) ; 
					corr_Hybrid = Engine::Utilities::correlation(&obs_intensities[0], &marr_theoretical_hybrid[0], 0) ; 
					corr_HM = 0 ; 
					corr_CA = 0 ; 
				}
				else
				{
					corr_CS = 0 ; 
					corr_Hybrid = 0 ; 
					corr_HM = Engine::Utilities::correlation(&obs_intensities[0], &marr_theoretical_hm[0], 0 ) ; 
				    corr_CA = Engine::Utilities::correlation(&obs_intensities[0], &marr_theoretical_ca[0], 0) ; 
				}
				

				max_corr = Engine::Utilities::max(Engine::Utilities::max(Engine::Utilities::max(corr_HM, corr_CA), corr_CS), corr_Hybrid) ; 

				if(max_corr == corr_HM)
						glycan_type = Engine::MS2HCDScoring::GLYCAN_TYPE::HIGH_MANNOSE ;
				if(max_corr == corr_CA)
						glycan_type = Engine::MS2HCDScoring::GLYCAN_TYPE::COMPLEX_ASIALYLATED ; 
				if(max_corr == corr_CS)
						glycan_type = Engine::MS2HCDScoring::GLYCAN_TYPE::COMPLEX_SIALYLATED ; 
				if(max_corr == corr_Hybrid)
						glycan_type = Engine::MS2HCDScoring::GLYCAN_TYPE::HYBRID ; 
				
			}
			else
			{
				p_value = 1 ;  
				glycan_type = Engine::MS2HCDScoring::GLYCAN_TYPE::NA; 
			}

			return p_value ; 
		}

	

		double HCDScoring::DetermineY1IonThroughPeptideSearching(Engine::PeakProcessing::PeakData &hcd_data, Engine::MS2CIDScoring::CIDInformationRecord &cid_record) 
		{
			//----------- Function to return y1 ion by looking for peptide + GlucNac peak in HCD ---------- //
			double seq_mass, seq_plus_glcnac_mass, seq_plus_glcnac_mz ; 
			Engine::PeakProcessing::Peak peak ; 
			double max_peak_mz = 0.0 ; 
			double max_peak_intensity = 0.0; 
			int charge = (int) cid_record.mshort_cs ; 

			// Check if record already has a peptide associated.  The difference between this and CID
			// y1 determination is that this is done on HCD peak data which has a higher chance of 
			// observing y1
			if (cid_record.mdbl_seq_mass > 0 )
			{
				seq_mass = cid_record.mdbl_seq_mass ; 
				seq_plus_glcnac_mass = (double) seq_mass + marr_theoretical_peaks[2] ; 

				for (int cs = charge-1 ; cs >= 1 ; cs--)
				{
					seq_plus_glcnac_mz = seq_plus_glcnac_mass/cs  + mdbl_cc_mass ; 
					bool found = hcd_data.GetClosestPeakFromAll(seq_plus_glcnac_mz, MZ_ERROR, peak) ; 
					if (found) 
					{
						if (peak.mdbl_intensity> max_peak_intensity)
						{
							max_peak_intensity = peak.mdbl_intensity; 
 							max_peak_mz = peak.mdbl_mz ; 
						}

					}
				}
			}
			else
			{
				max_peak_mz = 0 ; 
			}

			return max_peak_mz ; 

		}
		
		double HCDScoring::DetermineY1IonThroughCorrelation(Engine::PeakProcessing::PeakData &hcd_peaks, Engine::PeakProcessing::PeakData &cid_peaks, double min_mz, double max_mz) 
		{
			//---------- Function to determine the y1 ion through comparison of hcd and cid peak intensities ---------- //

			// y1 is returned as
			//  - max peak in HCD of all common peaks between CID and HCD after 700
			//  - if no common peak is found, the max peak in the HCD post all the mono sacharide m/zs are returned i.e. > 700
			
			double set_precision = 100 ; 
			int num_hcd_points = hcd_peaks.GetNumPeaks() ; 
			int num_cid_points = cid_peaks.GetNumPeaks() ; 
			double hcd_intensity = 0.0 ; 
			double hcd_mz = 0.0 ; 
			double max_peak_intensity = 0.0 ; 
			double max_peak_mz = 0.0 ; 			
			typedef std::pair <int, double> this_hcd_index; 
			typedef std::pair <int, double> this_cid_index ; 
			std::map <int, double> hcd_index_intensity ; 
			std::map <int, double> cid_index_intensity ; 
			std::map <int, double>::iterator hcd_iter ; 
			std::map <int, double>::iterator cid_iter ; 

			Engine::PeakProcessing::Peak peak; 
			bool found_peak ; 			
			double temp_intensity ; 		
			
			//----------- HCD ----------//
			// First bin up hcd fragments
			int i = 0 ; 
			double delta_mz = 0.5 ; 
			int num_bins = (max_mz - min_mz)/delta_mz ; 
			for (int j = 0 ; j < hcd_peaks.GetNumPeaks(); j++)			
			{
				hcd_peaks.GetPeak(j, peak) ;			
				int this_index = (int)((peak.mdbl_mz-min_mz)/delta_mz) ;
				hcd_iter  = hcd_index_intensity.find(this_index) ; 
				if (hcd_iter == hcd_index_intensity.end()) 
				{
					hcd_index_intensity.insert(this_hcd_index(this_index, peak.mdbl_intensity)) ; 			
				}
				else
				{
					temp_intensity = hcd_index_intensity[hcd_iter->first] ; 
					hcd_index_intensity[hcd_iter->first] = temp_intensity + peak.mdbl_intensity ; 
				}			
				i++ ; 
			}

			//----------- CID comparison ----------//
			// Second bin up all cid fragments
			i = 0 ; 
			for (int j = 0 ; j < cid_peaks.GetNumPeaks(); j++)			
			{
				cid_peaks.GetPeak(j, peak) ; 
				int this_index = (int)((peak.mdbl_mz-min_mz)/delta_mz) ;
				cid_iter  = cid_index_intensity.find(this_index) ; 
				if (cid_iter == cid_index_intensity.end()) 
				{
					// new guy insert in cid map
					cid_index_intensity.insert(this_cid_index(this_index, peak.mdbl_intensity)) ; 			
					hcd_iter = hcd_index_intensity.find(this_index);
					if (hcd_iter == hcd_index_intensity.end())
					{
						//do nothing since not in HCD
						i++ ; 
					}
					else
					{
						// in HCD now so see if it's a prominent peak
						hcd_intensity = hcd_iter->second  ; 
						
						if (hcd_intensity > max_peak_intensity && peak.mdbl_mz > 700 )
						{
							max_peak_intensity = hcd_intensity ; 
							max_peak_mz = peak.mdbl_mz ; 
						}
					}
				}
				else
				{
					// guy is already present in cid_map so add this guy to already existing guy and recomputer hcd_cid product 
					temp_intensity = cid_index_intensity[cid_iter->first] ; 
					cid_index_intensity[cid_iter->first] = temp_intensity + peak.mdbl_intensity ; 
					hcd_iter = hcd_index_intensity.find(this_index);
					if (hcd_iter == hcd_index_intensity.end())
					{
						//do nothing
						i++ ; 
					}
					else
					{
						// See if the intensity is highest in the HCD scan
						hcd_intensity = hcd_iter->second ; //* cid_index_intensity[cid_iter->first] ; 
						if (hcd_intensity > max_peak_intensity && peak.mdbl_mz > 700 )
						{
							max_peak_intensity = hcd_intensity ; 
							max_peak_mz = peak.mdbl_mz ; 
						}
					}
				}			
			}

			// If no max peak was found, return largest peak in HCD outside of max_peak_mz
			if (max_peak_mz == 0)
			{
				found_peak = hcd_peaks.GetNextPeak(mdbl_max_mz, max_mz, peak);
				max_peak_mz =  peak.mdbl_mz ; 
			}




			//// Multiply intensities from cid_spectra and choose the max peak
			//i = 0 ; 
			//found_peak = cid_peaks.GetNextPeak(min_mz, max_mz, peak) ;			
			//while(found_peak)
			//{		
			//	int this_index = (int)((peak.mdbl_mz-min_mz)/delta_mz) ; 
			//	index_iter_map  = index_intensity_map.find(this_index) ; 
			//	if (index_iter_map == index_intensity_map.end()) 
			//	{
			//		//absent do nothing
			//		i++ ; 
			//	}
			//	else
			//	{
			//		// present, multiply both and compare to max peak
			//		temp_intensity = index_intensity_map[index_iter_map->first] ; 
			//		index_intensity_map[index_iter_map->first] = temp_intensity * peak.mdbl_intensity ; 
			//		temp_intensity = index_intensity_map[index_iter_map->first] ; 
			//		if (temp_intensity > max_peak_intensity)
			//		{
			//			max_peak_mz = peak.mdbl_mz ; 
			//			max_peak_intensity = temp_intensity ; 
			//		}
			//	}
			//	cid_peaks.RemovePeak(peak) ; 
			//	found_peak = cid_peaks.GetNextPeak(min_mz, max_mz, peak) ; 
			//}


			/*	

			//Create map from hcd_spectra
			int i = 0 ; 
			found_peak = hcd_peaks.GetNextPeak(min_mz, max_mz, peak) ; 
			while(found_peak)			
			{				
				double temp_mz1 = floor(peak.mdbl_mz*set_precision+0.5); 				
				temp_mz = temp_mz1/set_precision;
				if (peak.mdbl_intensity >max_peak_intensity)
				{
					max_peak_intensity = peak.mdbl_intensity ; 
				}
				mz_intensity_map.insert(Mz_Int_Pair(temp_mz, peak.mdbl_intensity)) ; 			
				hcd_peaks.RemovePeak(peak) ; 
				found_peak = hcd_peaks.GetNextPeak(min_mz, max_mz, peak) ; 
				i++ ; 				
			}		*/

			
			



			/*// Multiply intensities from cid_spectra and choose the max peak
			i = 0 ; 
			found_peak = cid_peaks.GetNextPeak(min_mz, max_mz, peak) ;			
			while(found_peak)
			{				
				double temp_mz1 = floor(peak.mdbl_mz*set_precision+0.5);
				temp_mz = temp_mz1/set_precision;
				iter_map = mz_intensity_map.find(temp_mz) ; 					
				if (iter_map == mz_intensity_map.end())
				{
					// absent, do nothing	
					i++ ; 
					
				}
				else
				{
					//present, multiple both and compare to max peak
					hcd_mz = (double) iter_map->first ; 									
					temp_intensity = mz_intensity_map[iter_map->first]; 
					mz_intensity_map[iter_map->first] = temp_intensity * peak.mdbl_intensity  ;
					temp_intensity = mz_intensity_map[iter_map->first]; 					
					if (temp_intensity > max_peak_intensity)
					{
						max_peak_mz = hcd_mz ; 
						max_peak_intensity = temp_intensity ; 
					}
				}
				cid_peaks.RemovePeak(peak) ; 
				found_peak = cid_peaks.GetNextPeak(min_mz, max_mz, peak) ; 
			}*/

			// Return 
			return max_peak_mz ;
		}
	}
}