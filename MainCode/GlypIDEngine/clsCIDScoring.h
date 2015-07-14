// Written by Anoop Mayampurath, Indiana University
#pragma once
#using <mscorlib.dll>
#include <vector>
#include <fstream>
#include "clsPeak.h"
#include "clsHornTransformResults.h"
#include "clsHornTransformParameters.h"
#include "clsCIDScoringScanResults.h"
#include "clsScoringParameters.h" 
#include "clsSequence.h"
#include "SequenceManager/Sequence.h"
#include "PeakProcessor/PeakData.h"
#include "HornTransform/MassTransform.h" 
#include "MS2CIDScore/CIDScoring.h"
#include "clsElementIsotopes.h" 
#include "HornTransformTheoreticalProfile/AtomicInformation.h" 
#include "clsETDScoringScanResults.h"


namespace GlypID
{
	namespace CIDScoring
	{
		public __gc class clsCIDScoring
		{
			int mint_percent_done ; 
			int mint_msn_scan ; 
			int mint_parent_scan ; 
			int mint_msn_level ; 
			int mint_parent_level ; 
			double mdbl_parent_mz ; 	
			double mdbl_cc_mass ; 


			System::String *mstr_status_mesg ; 
			Engine::MS2CIDScoring::CIDScoring __nogc *mobj_cid_scoring ; 
			Engine::HornTransform::MassTransform __nogc *mobj_transform ; 			
			GlypID::HornTransform::clsHornTransformParameters* mobj_transform_parameters ; 
			GlypID::Scoring::clsScoringParameters* mobj_scoring_parameters ; 

			// Mass Transform Options
			void SetIsotopeFitOptions(System::String *str_averagine, System::String *str_tag, bool thrash_or_not, bool complete_fit) ; 
			void SetOptions(short max_charge, double max_mw, double max_fit, double min_s2n, double cc_mass, 
				double delete_threshold_intensity, double min_theoretical_intensity_for_score, short num_peaks_for_shoulder, 
				bool use_caching, bool o16_o18_media, bool check_against_charge_1) ; 
			
			
		public:

			__property GlypID::HornTransform::clsHornTransformParameters* get_TransformParameters()
			{
				return mobj_transform_parameters ; 
			}
			
			__property void set_TransformParameters(GlypID::HornTransform::clsHornTransformParameters* transform_parameters)
			{
				mobj_transform_parameters = dynamic_cast<GlypID::HornTransform::clsHornTransformParameters *>(transform_parameters->Clone()) ; 
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

			__property GlypID::Scoring::clsScoringParameters* get_ScoringParameters()
			{
				return mobj_scoring_parameters ; 
			}

			__property void set_ScoringParameters(GlypID::Scoring::clsScoringParameters *scoring_parameters)
			{				
				mobj_scoring_parameters = dynamic_cast<GlypID::Scoring::clsScoringParameters *>(scoring_parameters->Clone()) ; 				
			//	mobj_cid_scoring->Init() ; 
				//nothing to do here
			}

			__property int get_PercentDone()
			{
				return mint_percent_done ; 
			}

			__property System::String* get_StatusMessage()
			{
				return mstr_status_mesg ; 
			}

			clsCIDScoring(void);
			~clsCIDScoring(void);			
			bool ScoreCIDSpectra(GlypID::Peaks::clsPeak* (&hcd_peaks) __gc[], float (&mzs) __gc [], float (&intensities) __gc [], GlypID::HornTransform::clsHornTransformResults* (&precursor_records)  __gc[], GlypID::CIDScoring::clsCIDScoringScanResults* (&cid_scoring_results) __gc[])  ; 
			void SearchForY1Ion(GlypID::Sequence::clsSequence* (&peptides) __gc[],GlypID::Peaks::clsPeak* (&cid_peaks) __gc[], GlypID::CIDScoring::clsCIDScoringScanResults* (&cid_scoring_results) __gc[]) ; 
			void SetScanDetails(int msNScan, int parent_scan, int msNLevel, int parentLevel, double parentMz ); 			

			
			
			
		};
	}
}