// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)


#include "clshorntransform.h"
#include "clsPeakProcessor.h"
#include "GlypIDEngineUtils.h"

namespace GlypID
{
	namespace HornTransform
	{

		clsHornTransform::clsHornTransform(void)
		{
			try
			{
				mobj_transform = new Engine::HornTransform::MassTransform() ; 
			}
			catch (char *mesg)
			{
				System::String *exception_msg = new System::String(mesg) ; 
				throw new System::Exception(exception_msg) ; 
			}
		}

		clsHornTransform::~clsHornTransform(void)
		{
		}
		

		void clsHornTransform::SetOptions(short max_charge, double max_mw, double max_fit, double min_s2n, double cc_mass, 
			double delete_threshold_theoretical_intensity, double min_theoretical_intensity_for_score, short num_peaks_for_shoulder, 
			bool use_caching, bool o16_o18_media, bool check_against_charge_1)
		{
			mobj_transform->SetOptions(max_charge, max_mw, max_fit, min_s2n, cc_mass, delete_threshold_theoretical_intensity,
				min_theoretical_intensity_for_score, num_peaks_for_shoulder, check_against_charge_1, use_caching, o16_o18_media) ; 
		}

		void clsHornTransform::SetIsotopeFitOptions(System::String *str_averagine, System::String *str_tag, bool thrash_or_not, 
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

		void clsHornTransform::AllocateValuesToTransform(float parent_mz, int abundance, short (&charges) __gc[], GlypID::HornTransform::clsHornTransformResults* (&parentTransformResults) __gc [])
		{
			int num_charges = charges->Length ; 
			Engine::HornTransform::IsotopeFitRecord transformRecord ; 
			std::vector<Engine::HornTransform::IsotopeFitRecord> vectTransformRecord ; 
			double cc_mass  = 1.00727638 ;	

			for (int i = 0 ; i < num_charges ; i++)
			{
				transformRecord.mshort_cs = charges[i];                            
                transformRecord.mdbl_mono_mw = (parent_mz - cc_mass) * charges[i];
				transformRecord.mdbl_mz = parent_mz ; 
                transformRecord.mdbl_fit = 1;
                transformRecord.mint_abundance = abundance ;                             
                vectTransformRecord.push_back(transformRecord) ; 
			}
			
			// Done with the transform. Lets copy them all to the given memory structure.
			int numTransformed = (int) vectTransformRecord.size() ; 				
			parentTransformResults = new clsHornTransformResults* __gc [numTransformed] ; 
			for (int transformNum = 0 ; transformNum < numTransformed ; transformNum++)
			{
				parentTransformResults[transformNum] = new clsHornTransformResults() ; 
				parentTransformResults[transformNum]->Set(vectTransformRecord[transformNum]) ; 
			}
			mint_percent_done = 100 ; 
		}
      

		bool clsHornTransform::FindPrecursorTransform(float background_intensity, float min_peptide_intensity, 
			float (&parent_mzs) __gc[], float (&parent_intensities) __gc [], 
			GlypID::Peaks::clsPeak* (&parent_peaks) __gc [], float parent_mz, 
			GlypID::HornTransform::clsHornTransformResults* (&parentTransformResults) __gc [])
		{
			mint_percent_done = 0 ; 
			bool found = false ; 
			Engine::HornTransform::IsotopeFitRecord transformRecord ; 
			Engine::PeakProcessing::PeakData parentPeakData ; 
			std::vector<Engine::HornTransform::IsotopeFitRecord> vectTransformRecord ; 
			std::vector<double> vectMzs ;
			std::vector<double> vectIntensities ;
			
			if (parent_mzs->Length == 0)
				return false; 

			mobj_transform->Reset() ; 

		
			GlypID::Utils::SetData(vectMzs, parent_mzs) ; 
			GlypID::Utils::SetData(vectIntensities, parent_intensities) ; 
			GlypID::Utils::SetPeaks(parentPeakData, parent_peaks) ;
			parentPeakData.mptr_vect_mzs = &vectMzs ; 
			parentPeakData.mptr_vect_intensities = &vectIntensities ;

			found = mobj_transform->FindPrecursorTransform(parentPeakData, parent_mz, &vectTransformRecord, background_intensity,
				min_peptide_intensity) ; 

			if (found)
			{				
				// Done with the transform. Lets copy them all to the given memory structure.
				int numTransformed = (int) vectTransformRecord.size() ; 				
				parentTransformResults = new clsHornTransformResults* __gc [numTransformed] ; 
				for (int transformNum = 0 ; transformNum < numTransformed ; transformNum++)
				{
					parentTransformResults[transformNum] = new clsHornTransformResults() ; 
					parentTransformResults[transformNum]->Set(vectTransformRecord[transformNum]) ; 
				}
				mint_percent_done = 100 ; 
			}
		

			return found ; 
		}

		void clsHornTransform::PerformTransform(float background_intensity, float min_peptide_intensity, float (&mzs) __gc [], float (&intensities) __gc [], 
			GlypID::Peaks::clsPeak* (&peaks) __gc [], GlypID::HornTransform::clsHornTransformResults* (&transformResults) __gc [])
		{
			mint_percent_done = 0 ;
			std::vector<double> vectMzs ;
			std::vector<double> vectIntensities ;
			int numPoints = mzs->Length ; 

			if (mzs->Length == 0)
				return ; 

			mobj_transform->Reset() ; 

			// mzs should be in sorted order
			double minMZ = mzs[0] ; 
			double maxMZ = mzs[numPoints-1] ; 
			GlypID::Utils::SetData(vectMzs, mzs) ; 
			GlypID::Utils::SetData(vectIntensities, intensities) ; 

			Engine::PeakProcessing::PeakData peakData ; 
			Engine::PeakProcessing::Peak currentPeak ; 
			Engine::HornTransform::IsotopeFitRecord transformRecord ; 
			GlypID::Utils::SetPeaks(peakData, peaks) ; 
			peakData.mptr_vect_mzs = &vectMzs ; 
			peakData.mptr_vect_intensities = &vectIntensities ; 

			if (mobj_transform_parameters->get_UseMZRange())
			{
				minMZ = mobj_transform_parameters->get_MinMZ() ; 
				maxMZ = mobj_transform_parameters->get_MaxMZ() ; 
			}

			bool found = peakData.GetNextPeak(minMZ, maxMZ, currentPeak) ; 
			double fwhm_SN = currentPeak.mdbl_FWHM ; 

		
			std::vector<Engine::HornTransform::IsotopeFitRecord> vectTransformRecord ; 

			int numTotalPeaks = peakData.GetNumPeaks() ; 
			mstr_status_mesg = S"Performing Horn Transform on peaks" ; 
			//mobj_transform->mbln_debug = true ;
			while(found)
			{
				int numPeaksLeft = peakData.GetNumUnprocessedPeaks() ; 
				mint_percent_done = (100* (numTotalPeaks-numPeaksLeft))/ numTotalPeaks ; 
				if (mint_percent_done % 5 == 0)
				{
					mstr_status_mesg = System::String::Concat(S"Done with ", System::Convert::ToString(numTotalPeaks-numPeaksLeft), S" of ", System::Convert::ToString(numTotalPeaks), S" peaks.") ; 
				}
				if (currentPeak.mdbl_intensity < min_peptide_intensity)
					break ; 

				bool found_transform = false ;
				found_transform = mobj_transform->FindTransform(peakData, currentPeak, transformRecord, background_intensity) ; 
				if (found_transform && transformRecord.mshort_cs <= mobj_transform_parameters->get_MaxCharge())
				{
					vectTransformRecord.push_back(transformRecord) ;
				}
				found = peakData.GetNextPeak(minMZ, maxMZ, currentPeak) ; 
			}
			mint_percent_done = 100 ; 

			// Done with the transform. Lets copy them all to the given memory structure.
			int numTransformed = (int) vectTransformRecord.size() ; 
			mstr_status_mesg = System::String::Concat(S"Done with Mass Transform. Found ", System::Convert::ToString(numTransformed), S" features") ; 
			Console::WriteLine(mstr_status_mesg) ; 

			transformResults = new clsHornTransformResults* __gc [numTransformed] ; 
			for (int transformNum = 0 ; transformNum < numTransformed ; transformNum++)
			{
				transformResults[transformNum] = new clsHornTransformResults() ; 
				transformResults[transformNum]->Set(vectTransformRecord[transformNum]) ; 
			}
			mint_percent_done = 100 ; 
		}
	
  }
}