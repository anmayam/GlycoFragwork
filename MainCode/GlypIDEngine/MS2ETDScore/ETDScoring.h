// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University


#include "../PeakProcessor/PeakProcessor.h"
#include "../PeakProcessor/PeakData.h"
#include "../Utilities/util.h"
#include "../Utilities/system.h" 
#include <vcclr.h>


namespace Engine
{
	namespace MS2ETDScoring
	{		

		class ETDScoring
		{
				// Binning start point
				double mdbl_bin_start ; 
				// Binning unit size
				double mdbl_bin_size ; 
				// number of sections to normalize
				double mint_num_sections ; 
				// min mz of fragments
				double mdbl_min_mz ; 
				//max mz of fragments
				double mdbl_max_mz ; 
		
			public:
				// charge carrier mass
				double mdbl_cc_mass; 
				// constructor
				ETDScoring(void); 
				//destructor
				~ETDScoring(void) ; 
				// Main scoring function
				double CalculateETDScore(std::vector<double> &obs_mzs,std::vector<double> &obs_intensities, std::vector<double> &th_mzs, std::vector<double> &th_intensities) ;
				double CalculateETDScore2(std::vector<double> &obs_mzs,std::vector<double> &obs_intensities, std::vector<double> &th_mzs, std::vector<double> &th_intensities) ; 
				// Function to modify sequence to include glycan modification
				int AddGlycanModificationToSequence(char *sym, std::string inp_sequence, std::string &out_sequence, int search_start_position) ; 				
				void AddCarboModificationToSequence(std::string inp_sequence, std::string &out_sequence) ; 
				void AddOGlycanModificationToSequence(char *sym, std::string inp_sequence, std::string &out_sequence, int position) ;

				// Function to normalize the spectra but by SEQUEST standards
				void BinNormalize(std::vector<double> &intensities, std::vector<double> &mzs) ; 
				void BinNormalizeTheoretical(std::vector<double> &intensities, std::vector <double> &mzs) ;
				// Function to convolve theoretical and observed intensities
				double ConvolveIntensities(std::vector<double> &u, std::vector<double> &v);
				double CalculateMorpheusScore(std::vector<double> &u, std::vector <double> &v, std::vector<double> &mzs) ;
				double CalculateWeightedScore(std::vector<double> &u, std::vector <double> &v, std::vector<double> &mzs) ;
		


				
		};
	}
}