// Written by Anoop Mayampurath, Indiana University
#pragma once
#using <mscorlib.dll>
#include <vector>
#include <fstream>
#include "clsPeak.h"
#include "clsHornTransformResults.h"
#include "clsHornTransformParameters.h"
#include "clsScoringParameters.h" 
#include "clsETDScoringScanResults.h"
#include "clsCIDScoring.h"
#include "clsCIDScoringScanResults.h"
#include "clsGlycan.h"
#include "PeakProcessor/PeakData.h"
#include "HornTransform/MassTransform.h" 
#include "MS2ETDScore/ETDScoring.h"
#include "GlycanCompositionManager/GlycanProcessor.h"
#include "clsElementIsotopes.h" 
#include "HornTransformTheoreticalProfile/AtomicInformation.h" 

namespace GlypID
{
	namespace ETDScoring
	{
		public __gc class clsETDScoring
		{
			int mint_msn_scan ; 
			int mint_parent_scan ; 
			int mint_msn_level ; 
			int mint_parent_level ; 
			double mdbl_parent_mz ; 

			Engine::MS2ETDScoring::ETDScoring __nogc *mobj_etd_scoring ; 
			Engine::HornTransform::MassTransform __nogc *mobj_transform ; 	
			MwtWinDll::MolecularWeightCalculator __gc *mobj_mwt; 
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
				
				// to do mobj_etd_scoring->SetOptions(mobj_scoring_parameters->get_MinNumPeaksToConsider(), mobj_scoring_parameters->get_MinHCDMz(), mobj_scoring_parameters->get_MaxHCDMz()) ; 
			}

			clsETDScoring(void) ; 
			~clsETDScoring(void) ; 

			void SetScanDetails(int msNScan, int parent_scan, int msNLevel, int parentLevel, double parentMz ) ; 

			double ScoreETDSpectra(GlypID::Peaks::clsPeak* (&etd_peaks) __gc[],  GlypID::ETDScoring::clsETDScoringScanResults* etd_scoring_result) ;

			double ScoreETDSpectraOLinked(GlypID::Peaks::clsPeak* (&etd_peaks) __gc[],  GlypID::ETDScoring::clsETDScoringScanResults* etd_scoring_result) ;

			double ScoreMS3CIDSpectra(GlypID::Peaks::clsPeak* (&etd_peaks) __gc[],  GlypID::ETDScoring::clsETDScoringScanResults* etd_scoring_result) ;

			void PreprocessETDSpectra(GlypID::Peaks::clsPeak* (etd_peaks) __gc[],  Engine::PeakProcessing::PeakData& etd_proc_peakData, double parent_mz ) ;

				

		};
	}
}