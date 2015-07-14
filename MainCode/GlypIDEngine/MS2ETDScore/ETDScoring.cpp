// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University

#include "./ETDscoring.h"
#include "./../Utilities/DeconException.h" 
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctype.h>
#include <io.h>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vcclr.h>

namespace Engine
{
	namespace MS2ETDScoring
	{		

		ETDScoring::ETDScoring(void)
		{
			//mdbl_GlcNac_Mass = 204.09;	
			//mdbl_cc_mass = 1.00727638;			
			//mobj_Mwt_Win = new MwtWinDll::MolecularWeightCalculator() ; 
			//std::string test = "TDNESALPVTVLSAEDIAK" ; 
			// 
			//

			//mobj_Mwt_Win->Peptide->SetSequence(test.c_str());
			//MwtWinDll::MWPeptideClass::udtFragmentationSpectrumDataType udtTest[]  ; 
			//mobj_Mwt_Win->Peptide->GetFragmentationMasses(&udtTest) ; 

			mdbl_bin_start = 0 ; 
			mdbl_bin_size = 1 ; 
			mint_num_sections = 10 ; 
			mdbl_min_mz = 0 ; 
			mdbl_max_mz = 2000  ; 

				
		}

		ETDScoring::~ETDScoring(void)
		{
		}
		

		void ETDScoring::AddOGlycanModificationToSequence(char *sym, std::string inp_sequence, std::string &out_sequence, int position) 
		{
			// Function to insert O-glycan symbol into the sequence
			// position is 1 based
			out_sequence = inp_sequence ; 
			out_sequence.insert(position, sym) ; 

		}
		

		int ETDScoring::AddGlycanModificationToSequence(char *sym, std::string inp_sequence, std::string &out_sequence, int search_start_position)
		{
			//--------- Function that replaces the frist glycosational site N from the search_start_position with N"symbol"  ---------//
			// Eg. replace "...VLNXS/TLP.." with "..VLN^XS/TLP.."
			// But checks if the glycosylation site happens to be at the end, due to misscleavage
			// Returns position on the N within the peptide
			char *seq = new char[inp_sequence.size() + 1]; 
			seq[inp_sequence.size()] = 0 ; 
			memcpy(seq, inp_sequence.c_str(), inp_sequence.size()) ;
			int pos = 0 ; 
			bool found = false ; 
			for (int i = search_start_position ; i < (int)inp_sequence.size()-2 ; i++)
			{
				if (seq[i] == 'N')
				{
					pos = i ; 
					if(seq[i+2] == 'S' || seq[i+2] == 'T')
					{
						if(seq[i+1] != 'P')
						{
							pos = i ; 
							found = true ; 
							break ; 
						}
					}
				}
			}

			if (found)
			{
				out_sequence = inp_sequence ; 
				out_sequence.insert(pos+1, sym) ; 
			}
			else
			{
				if (pos == inp_sequence.size() -1)
				{
					// miscleaved peptide
					out_sequence = inp_sequence ; 
					out_sequence.insert(pos, sym) ; 
				}
				// check once
				if ((seq[inp_sequence.size()-2] == 'N') && (seq[inp_sequence.size()-1] != 'P'))
				{
					pos = inp_sequence.size()-2;
					// miscleaved peptide
					out_sequence = inp_sequence ; 
					out_sequence.insert(pos+1, sym) ; 
				}

			}

			return pos ; 
		}

	
		void ETDScoring::AddCarboModificationToSequence(std::string inp_sequence, std::string &out_sequence) 
		{
			//--------- Function that replaces the C with C! corresponding to carbamidomethylation ---------//
			// Eg. replace "...LCPC.." with "..LC!PC!.."


			size_t found ; 
			found = inp_sequence.find('C') ; 
			
			while (found != string::npos)
			{
				inp_sequence.insert(found+1, "!");
				found = inp_sequence.find('C', found+1);
			}

			out_sequence = inp_sequence; 


			/*char *seq = new char[inp_sequence.size() + 1]; 
			seq[inp_sequence.size()] = 0 ; 
			memcpy(seq, inp_sequence.c_str(), inp_sequence.size()) ;
			int pos = 0 ; 
			bool found = false ; 
			int count =0; 
			for (int i = 0 ; i < (int)inp_sequence.size()-2 ; i++)
			{
				if (seq[i] == 'C')
				{
					out_sequence[count] = 'C') ; 
					count++; 
					out_sequence[count] = '!') ; 
				}
					
			}*/

			

		}

		double ETDScoring::CalculateETDScore2(std::vector<double> &obs_mzs,std::vector<double> &obs_intensities, std::vector<double> &th_mzs, std::vector<double> &th_intensities)
		{
			double etd_score = 0 ;

			// Bin up and normalize observed fragment in a SEQUEST Fashion
			BinNormalizeTheoretical(th_intensities, th_mzs) ; 

			

			// Bin up experimental fragment
			//BinNormalizeTheoretical(*pk_data.mptr_vect_intensities, *pk_data.mptr_vect_mzs) ; 
			BinNormalizeTheoretical(obs_intensities, obs_mzs) ; 

			// 			etd_score = CalculateMorpheusScore(th_intensities, obs_intensities,  th_mzs) ; 
			etd_score = CalculateWeightedScore(th_intensities, obs_intensities, th_mzs) ; 

			return etd_score ; 
		}




		double ETDScoring::CalculateETDScore(std::vector<double> &obs_mzs,std::vector<double> &obs_intensities, std::vector<double> &th_mzs, std::vector<double> &th_intensities)
		{

			// Score theoretical against observed peak profiles using a SEQUEST Xcorr approach
			// Scoring taken from Yong F. Li's scripts
			
			double etd_score = 0 ; 
			// Bin up and normalize observed fragment in a SEQUEST Fashion
			BinNormalizeTheoretical(th_intensities, th_mzs) ; 

			// Bin up experimental fragment
			BinNormalizeTheoretical(obs_intensities, obs_mzs) ; 
			

			// Convolve the intensities
			etd_score = ConvolveIntensities(th_intensities, obs_intensities) ; 

			return etd_score ; 
		}

		double ETDScoring::CalculateWeightedScore(std::vector<double> &u, std::vector <double> &v, std::vector<double> &mzs) 
		{

			int num_common_peaks = 0 ;
			double total_intensity = 0 ; 
			double explained_intensity = 0 ; 
			double ui, vi ; 
			double common_mz = 0 ; 

			double w_score = 0 ;

			if (u.size()!= v.size())
			{
				throw "Unequal lengths" ; 
			}

			int N = u.size() ; 
			bool debug = true; 

			for (int i=0; i<N; i++)
			{
				total_intensity += v[i] ; 

				if (u[i] >0  && v[i] >0)
				{
					num_common_peaks++ ; 
					common_mz = mzs[i] ; 
					ui = u[i] ; 
					vi = v[i] ; 
					explained_intensity += v[i] ; 

					w_score += ui; 

				}
			}
			w_score = w_score + (explained_intensity/total_intensity) ; 

			return w_score ; 

			

		}

		double ETDScoring::CalculateMorpheusScore(std::vector<double> &u, std::vector <double> &v, std::vector <double> &mzs)
		{
			int num_common_peaks = 0 ;
			double total_intensity = 0 ; 
			double explained_intensity = 0 ; 
			double ui, vi ; 
			double common_mz = 0 ; 

			double m_score = 0 ;

			double m_score2 = 0 ;

			if (u.size()!= v.size())
			{
				throw "Unequal lengths" ; 
			}
			int N = u.size() ; 
			bool debug = true; 

			for (int i=0; i<N; i++)
			{
				total_intensity += v[i] ; 

				if (u[i] >0  && v[i] >0)
				{
					num_common_peaks++ ; 
					common_mz = mzs[i] ; 
					ui = u[i] ; 
					vi = v[i] ; 
					explained_intensity += v[i] ; 

					m_score2 += log(v[i]); 

				}
			}
			m_score = num_common_peaks + (explained_intensity/total_intensity) ; 

			return m_score ; 


		}




		double ETDScoring::ConvolveIntensities(std::vector<double> &u, std::vector<double> &v)
		{
			double score = 0 ; 
			if (u.size()!= v.size())
			{
				throw "Unequal lengths" ; 
			}

			int N = u.size() ; 
			double num_sum = 0 ; 
			double den_sum1 = 0 ; 
			double den_sum2 =0 ; 

			for (int i=0; i<N; i++)
			{
				num_sum += u[i]*v[i] ; 
			}

			for (int tau = -75; tau < 75 ; tau++)
			{
				den_sum1 =0 ;
				for (int i=0; i<N ; i++)
				{
					double term; 
					if ((i-tau <0) || (i-tau >= N))
					{
						term = 0 ; 
					}
					else
					{
						term = u[i] * v[i-tau]; 
					}
					den_sum1 += term; 
				}
				den_sum2 += den_sum1 ; 
			}
			score = 0.015 *num_sum/den_sum2 ; 
			return score ; 

		}

		


		void ETDScoring::BinNormalizeTheoretical(std::vector<double> &intensities, std::vector <double> &mzs)
		{
		
			// Init base sepctrum
			std::vector <double> base_mzs ; 
			std::vector <double> base_intensities ; 
			std::map <int, double> mz_index_intensity ; 
			std::map <int, double>::iterator map_iter ;
			typedef std::pair <int, double> map_pair; 

			int num_bins =  int ((mdbl_max_mz - mdbl_min_mz)/mdbl_bin_size) ; 
			for(int j= 0; j < num_bins ; j++)
			{
				mz_index_intensity.insert(map_pair(j,0)); 
			}

			bool debug = true ; 
			// Binning - data is binned into m/z bin	
			for(int j= 0 ; j < (int)intensities.size() ; j++)
			{
				if (abs(mzs.at(j)-646) < 1)
				{
					debug = false ; 
				}

				int this_index = (int)((mzs.at(j)-mdbl_min_mz)/mdbl_bin_size) ;
				map_iter  = mz_index_intensity.find(this_index) ; 
				if (map_iter != mz_index_intensity.end())
				{
					double temp_intensity = mz_index_intensity[map_iter->first] ; 
					mz_index_intensity[map_iter->first] = temp_intensity + intensities.at(j) ; 
				}
			}

			// Return
			intensities.clear() ; 
			mzs.clear() ; 
			for (map_iter = mz_index_intensity.begin(); map_iter != mz_index_intensity.end() ;++map_iter)
			{
				int index = map_iter->first ; 
				double this_mz = (index *mdbl_bin_size) + mdbl_min_mz ; 
				mzs.push_back(this_mz) ; 
				intensities.push_back(map_iter->second) ; 
			}
		}

		

		void ETDScoring::BinNormalize(std::vector<double> &intensities, std::vector<double> &mzs)
		{
			// ----------- Sequest-style binning and normalizing ---------- //

			// Init
			std::map <int, double> mz_index_intensity ; 			
			std::map <int, double> mz_index_max_intensity ; 
			std::map <int, double>::iterator map_iter ; 
			std::map <int, double>::iterator map_max_iter ; 
			typedef std::pair <int, double> map_pair; 

			int num_bins =  int ((mdbl_max_mz - mdbl_min_mz)/mdbl_bin_size) ; 
			for(int j= 0; j < num_bins ; j++)
			{
				mz_index_intensity.insert(map_pair(j,0)); 
			}

			// Binning - data is binned into 1 m/z bin					
			for (int j = 0 ; j < (int)intensities.size() ; j++)			
			{				
				int this_index = (int)((mzs.at(j)-mdbl_min_mz)/mdbl_bin_size) ;
				map_iter  = mz_index_intensity.find(this_index) ; 
				if (map_iter != mz_index_intensity.end()) 
				{					
					double temp_intensity = mz_index_intensity[map_iter->first] ; 
					mz_index_intensity[map_iter->first] = temp_intensity + sqrt(intensities.at(j)) ; 
				}		
			}

			// Normalizing - spectra is divided into 10 sections.  All intensities in a section are normalized so that the max is 50.
			double section_size = (mdbl_max_mz - mdbl_min_mz)/mint_num_sections ; 			
			int count = 0 ; 
			// Finding max in each section
			while(count < mint_num_sections)
			{
				int i_start = count*section_size  ; 
				int i_stop = i_start + section_size ; 
				double max_intensity = 0 ; 
				//buf gree max_intensity = (double) std::min_element(intensities.begin(), intensities.end()) ; 

				for (int i = i_start ; i < i_stop ; i++)
				{				
					if (mz_index_intensity[i] > max_intensity)
					{
						max_intensity = mz_index_intensity[i] ; 
					}					
				}
				mz_index_max_intensity.insert(map_pair(count, max_intensity)) ; 
				count++ ; 
			}
			// now normalize
			intensities.clear() ; 
			count = 0 ; 
			while(count < mint_num_sections)
			{
				int i = 0 ;
				double m = 0; 
				map_max_iter = mz_index_max_intensity.find(count) ; 
				if (map_max_iter != mz_index_max_intensity.end())
				{
					if (map_max_iter->second ==0)
						m = 1; 
					else
						m = map_max_iter->second ; 
				}else
				{
					m =1; 
				}
				int i_start = count*section_size  ; 
				int i_stop = i_start + section_size ; 
				for (int i = i_start ; i < i_stop ; i++)
				{
					mz_index_intensity[i] = mz_index_intensity[i]*50/m; 				    
				}
				count++ ; 
			}

			// Return 
			intensities.clear() ; 
			mzs.clear() ; 

			for (map_iter = mz_index_intensity.begin(); map_iter != mz_index_intensity.end() ;++map_iter)
			{
				int index = map_iter->first ; 
				double this_mz = (index *mdbl_bin_size) + mdbl_min_mz ; 
				mzs.push_back(this_mz) ; 
				intensities.push_back(map_iter->second) ; 
			}


/*
			for (int i = 0; i < mz_index_intensity.size()-1 ; i++)
			{
				map_iter = mz_index_intensity.find(i) ;
				int index = map_iter->first; 
				double this_mz = (index *mdbl_bin_size) + mdbl_min_mz ; 
				mzs.push_back(this_mz) ; 
				intensities.push_back(map_iter->second) ; */
			//}
		
		}
		
	}
}

