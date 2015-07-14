// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#pragma once
#using <mscorlib.dll>
#include <vector>
#include <fstream>
#include "clsPeak.h"
#include "clsHornTransformResults.h"
#include "clsHornTransformParameters.h"
#include "PeakProcessor/PeakData.h"
#include "HornTransform/MassTransform.h" 
#include "clsElementIsotopes.h" 
#include "HornTransformTheoreticalProfile/AtomicInformation.h" 
namespace GlypID
{
	namespace HornTransform
	{
		public __gc class clsHornTransform
		{
			int mint_percent_done ; 
			System::String *mstr_status_mesg ; 
			Engine::HornTransform::MassTransform __nogc *mobj_transform ; 
			clsHornTransformParameters* mobj_transform_parameters ; 
			void SetIsotopeFitOptions(System::String *str_averagine, System::String *str_tag, bool thrash_or_not, bool complete_fit) ; 
			void SetOptions(short max_charge, double max_mw, double max_fit, double min_s2n, double cc_mass, 
				double delete_threshold_intensity, double min_theoretical_intensity_for_score, short num_peaks_for_shoulder, 
				bool use_caching, bool o16_o18_media, bool check_against_charge_1) ; 
		public:

			__property clsHornTransformParameters* get_TransformParameters()
			{
				return mobj_transform_parameters ; 
			}
			__property void set_TransformParameters(clsHornTransformParameters* transform_parameters)
			{
				mobj_transform_parameters = dynamic_cast<clsHornTransformParameters *>(transform_parameters->Clone()) ; 
				SetOptions(mobj_transform_parameters->get_MaxCharge(), mobj_transform_parameters->get_MaxMW(), 
					mobj_transform_parameters->get_MaxFit(), mobj_transform_parameters->get_MinS2N(), mobj_transform_parameters->get_CCMass(),
					mobj_transform_parameters->get_DeleteIntensityThreshold(), mobj_transform_parameters->get_MinIntensityForScore(),
					mobj_transform_parameters->get_NumPeaksForShoulder(),
					mobj_transform_parameters->get_UseMercuryCaching(),  mobj_transform_parameters->get_O16O18Media(), 
					mobj_transform_parameters->get_CheckAllPatternsAgainstCharge1()) ; 
				SetIsotopeFitOptions(mobj_transform_parameters->get_AveragineFormula(), mobj_transform_parameters->get_TagFormula(),
					mobj_transform_parameters->get_ThrashOrNot(), mobj_transform_parameters->get_CompleteFit()) ; 
				const Engine::HornTransformTheoreticalProfile::AtomicInformation *ptr_atomic_info = mobj_transform_parameters->ElementIsotopeComposition->GetElementalIsotopeComposition() ; 
				if (mobj_transform != NULL)
				{
					mobj_transform->SetElementalIsotopeComposition(*ptr_atomic_info) ; 
					mobj_transform->SetIsotopeFitType((Engine::HornTransform::IsotopicFittingType) mobj_transform_parameters->get_IsotopeFitType()) ; 
				}
			}

			__property int get_PercentDone()
			{
				return mint_percent_done ; 
			}

			__property System::String* get_StatusMessage()
			{
				return mstr_status_mesg ; 
			}

			clsHornTransform(void);
			~clsHornTransform(void);
			void PerformTransform(float background_intensity, float min_peptide_intensity, float (&mzs) __gc [], float (&intensities) __gc [], 
				GlypID::Peaks::clsPeak* (&peaks) __gc [], GlypID::HornTransform::clsHornTransformResults* (&transformResults) __gc []) ; 
			bool FindPrecursorTransform(float background_intensity, float min_peptide_intensity, float (&parent_mzs) __gc[], float (&parent_intensities) __gc [], GlypID::Peaks::clsPeak* (&parent_peaks) __gc [], float parent_mz, GlypID::HornTransform::clsHornTransformResults* (&parentTransformResults) __gc []) ; 
			void AllocateValuesToTransform(float parent_mz, int abundance, short (&charges) __gc[], GlypID::HornTransform::clsHornTransformResults* (&parentTransformResults) __gc []) ; 


			
		};
	}
}