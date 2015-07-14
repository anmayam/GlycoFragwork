#include "./CIDScoring.h"
#include "../Utilities/util.h" 
#include "../Utilities/system.h"


#include <iostream>
#include <vector>
#include <algorithm>



namespace Engine
{
	namespace MS2CIDScoring
	{		
		CIDScoring::CIDScoring(void)
		{
			mdbl_GlcNac_Mass = 204.09;	
			mdbl_cc_mass = 1.00727638;	
		}

		CIDScoring::~CIDScoring(void)
		{
		}

		void CIDScoring::InitPTable()
		{
			mobj_Spectrum.sampleMS2Distribution() ;
		}

	

		double CIDScoring::CalculateCIDScore(Engine::PeakProcessing::PeakData &pk_data, int charge)
		{				
			std::vector <float> mass_diff ; 			
			mobj_Spectrum.allMassDif(mass_diff, 1) ; 
			std::sort(mass_diff.begin(), mass_diff.end()) ; 
			float score = 0  ; 
			if (charge == 1)
				score = mobj_Spectrum.computeMs2Score(mass_diff, pk_data,  charge) ; 
			else
				score = mobj_Spectrum.computeMs2Score(mass_diff, pk_data,  charge - 1) ; 		

			return score ; 
		}

		bool CIDScoring::CheckForOxoniumIons(Engine::PeakProcessing::PeakData &pk_data) 
		{
			return mobj_Spectrum.containsOxoniumion(pk_data) ; 
		}		

		double CIDScoring::CalculateCIDScorePValue(int score)
		{
			double val1 = (double)mobj_Spectrum.MS2ScoreToPValue(score);
			double val = -log(val1);
			return val ; 
		}

		double CIDScoring::DetermineY1IonThroughPeptideSearchingCID(Engine::PeakProcessing::PeakData &pk_data, Engine::MS2CIDScoring::CIDInformationRecord &cid_record)
		{
			//----------- Function to return y1 ion by looking for peptide + GlucNac peak in CID ---------- //
			double seq_mass, seq_plus_glcnac_mass, seq_plus_glcnac_mz ; 
			Engine::PeakProcessing::Peak peak ; 
			double max_peak_mz = 0.0 ; 
			double max_peak_intensity = 0.0; 
			int charge = (int)cid_record.mshort_cs ; 

			// Check if record already has a peptide associated with precursor
			if (cid_record.mdbl_seq_mass >0)
			{
				seq_mass = cid_record.mdbl_seq_mass ; 
				seq_plus_glcnac_mass = (double) seq_mass + mdbl_GlcNac_Mass ; 

				for (int cs = charge-1 ; cs >= 1 ; cs--)
				{
					seq_plus_glcnac_mz = seq_plus_glcnac_mass/cs  + mdbl_cc_mass ; 
					bool found = pk_data.GetClosestPeakFromAll(seq_plus_glcnac_mz, MZ_ERROR, peak) ; 
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
				max_peak_mz = 0; 				
			}

			return max_peak_mz ; 
		}
	}
}
