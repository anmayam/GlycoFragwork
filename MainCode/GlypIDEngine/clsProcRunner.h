 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)


#pragma once
#using <mscorlib.dll>
#include "clsPeakProcessor.h" 
#include "clsTransformResults.h" 
#include "clsHornTransformParameters.h" 
#include "clsHCDScoringResults.h"
#include "clsCIDHCDScoringResults.h"
#include "clsCIDETDScoringResults.h"
#include "clsRawData.h" 
#include "clsScoringParameters.h"
#include "clsCIDScoringResults.h" 


namespace GlypID
{
	public __value enum enmProcessState {IDLE = 0, RUNNING, COMPLETE, ERROR}; 
	
	public __gc class clsProcRunner
	{
		int mint_percent_done ; 
		int mint_current_scan ; 
		double mdbl_cc_mass ; 
		bool mbln_search_fasta ; 
		enmProcessState menm_state ; 
		System::String *mstr_file_name ; 	
		System::String *mstr_protein_list_name ; 
		System::String *mstr_glycan_composition_list_name ; 

		GlypID::Readers::FileType menm_file_type ;
		GlypID::Peaks::clsPeakProcessorParameters *mobj_peak_parameters ; 
		GlypID::HornTransform::clsHornTransformParameters *mobj_transform_parameters ; 		
		GlypID::Scoring::clsScoringParameters *mobj_scoring_parameters ; 
		GlypID::Results::clsTransformResults *mobj_transform_results ; 
		GlypID::Results::clsHCDScoringResults *mobj_hcd_scoring_results ; 
		GlypID::Results::clsCIDScoringResults *mobj_cid_scoring_results ; 
		GlypID::Results::clsCIDHCDCombinedScoringResults *mobj_cid_hcd_combined_scoring_results ;
		GlypID::Results::clsCIDETDCombinedScoringResults *mobj_cid_etd_combined_scoring_results ; 

		GlypID::Results::clsTransformResults* CreateTransformResults(System::String *file_name, GlypID::Readers::FileType file_type,
			 GlypID::Peaks::clsPeakProcessorParameters *peak_parameters, 
			GlypID::HornTransform::clsHornTransformParameters *transform_parameters);

		GlypID::Results::clsHCDScoringResults*  CreateHCDScoringResults(System::String *file_name, GlypID::Readers::FileType file_type, System::String *protein_file_name,
					System::String *glycan_composition_file_name, GlypID::Peaks::clsPeakProcessorParameters *peak_parameters,
				GlypID::HornTransform::clsHornTransformParameters *transform_parameters, GlypID::Scoring::clsScoringParameters *scoring_parameters, bool search_fasta) ; 

		GlypID::Results::clsCIDScoringResults* CreateCIDScoringResults(System::String *file_name, GlypID::Readers::FileType file_type, System::String *protein_file_name,
			System::String *glycan_composition_file_name , GlypID::Peaks::clsPeakProcessorParameters *peak_parameters,
				GlypID::HornTransform::clsHornTransformParameters *transform_parameters, GlypID::Scoring::clsScoringParameters *scoring_parameters, bool search_fasta) ; 

		
		GlypID::Results::clsCIDHCDCombinedScoringResults* CreateCIDHCDScoringResults(System::String *file_name, GlypID::Readers::FileType file_type, System::String *protein_file_name,
			System::String *glycan_composition_file_name , GlypID::Peaks::clsPeakProcessorParameters *peak_parameters, 
			GlypID::HornTransform::clsHornTransformParameters *transform_parameters, GlypID::Scoring::clsScoringParameters *scoring_parameters, bool search_fasta) ; 

		GlypID::Results::clsCIDETDCombinedScoringResults* CreateCIDETDScoringResults(System::String *file_name, GlypID::Readers::FileType file_type, 
		System::String *protein_file_name, System::String *glycan_composition_file_name , 
		GlypID::Peaks::clsPeakProcessorParameters *peak_parameters,
		GlypID::HornTransform::clsHornTransformParameters *transform_parameters, GlypID::Scoring::clsScoringParameters *scoring_parameters, bool search_fasta) ; 
	

		


	public:
		clsProcRunner(void);
		~clsProcRunner(void);
		void Reset()
		{
			mint_percent_done = 0 ; 
			mint_current_scan = 0 ; 
			menm_state = IDLE ; 
			mbln_search_fasta = false ; 
		}
		
		
		void CreateTransformResults() ; 
		void CreateHCDScoringResults() ; 
		void CreateCIDScoringResults() ; 
		void CreateCIDHCDScoringResults() ; 		
		void CreateCIDETDScoringResults() ; 

		
		__property int get_CurrentScanNum()
		{
			if (menm_state == IDLE)
				return  0 ;
			return  mint_current_scan ;
		}

		__property int get_PercentDone()
		{
			if (menm_state == IDLE)
				return  0 ;
			else if (menm_state == RUNNING)
				return  mint_percent_done ;
			else
				return 100 ; 
		}
		__property int get_ProcessState()
		{
			return menm_state ; 
		}
		__property System::String* get_FileName()
		{
			return mstr_file_name ; 
		}
		__property void set_FileName(System::String* file_name)
		{
			mstr_file_name = file_name ; 
		}	
		__property System::String* get_ProteinListFileName()
		{
			return mstr_protein_list_name ; 
		}
		__property void set_ProteinListFileName(System::String* file_name)
		{
			mstr_protein_list_name = file_name ; 
		}		
		__property System::String* get_GlycanCompositionFileName()
		{
			return mstr_glycan_composition_list_name ; 
		}
		__property void set_GlycanCompositionFileName(System::String* file_name)
		{
			mstr_glycan_composition_list_name = file_name ; 
		}		
		__property GlypID::Readers::FileType get_FileType()
		{
			return menm_file_type ; 
		}
		__property void set_FileType(GlypID::Readers::FileType file_type)
		{
			menm_file_type = file_type ; 
		}
		__property GlypID::Peaks::clsPeakProcessorParameters* get_PeakProcessorParameters()
		{
			return mobj_peak_parameters ; 
		}
		__property void set_PeakProcessorParameters(GlypID::Peaks::clsPeakProcessorParameters *peak_parameters)
		{
			mobj_peak_parameters = peak_parameters ; 
		}
		__property GlypID::HornTransform::clsHornTransformParameters* get_HornTransformParameters()
		{
			return mobj_transform_parameters ; 
		}
		__property void set_HornTransformParameters(GlypID::HornTransform::clsHornTransformParameters *transform_parameters)
		{
			mobj_transform_parameters = transform_parameters ; 
		}
		__property GlypID::Scoring::clsScoringParameters *get_ScoringParameters()
		{
			return mobj_scoring_parameters ; 
		}
		__property void set_ScoringParameters(GlypID::Scoring::clsScoringParameters *scoring_parameters)
		{
			mobj_scoring_parameters = scoring_parameters  ;
		}
		__property GlypID::Results::clsTransformResults* get_HornTransformResults()
		{
			return mobj_transform_results ; 
		}		
		__property void set_HCDScoringResults(GlypID::Results::clsHCDScoringResults *scoring)
		{
			mobj_hcd_scoring_results = scoring ; 
		}
		__property GlypID::Results::clsHCDScoringResults* get_HCDScoringResults() 
		{
			return mobj_hcd_scoring_results ; 
		}
		__property void set_CIDScoringResults(GlypID::Results::clsCIDScoringResults *scoring)
		{
			mobj_cid_scoring_results = scoring ; 
		}
		__property GlypID::Results::clsCIDScoringResults* get_CIDScoringResults() 
		{
			return mobj_cid_scoring_results ; 
		}
		__property void set_CIDHCDCombinedScoringResults(GlypID::Results::clsCIDHCDCombinedScoringResults *scoring)
		{
			mobj_cid_hcd_combined_scoring_results = scoring ; 
		}
		__property GlypID::Results::clsCIDHCDCombinedScoringResults* get_CIDHCDCombinedScoringResults() 
		{
			return mobj_cid_hcd_combined_scoring_results ; 
		}
		__property void set_CIDETDCombinedScoringResults(GlypID::Results::clsCIDETDCombinedScoringResults *scoring)
		{
			mobj_cid_etd_combined_scoring_results = scoring ; 
		}
		__property GlypID::Results::clsCIDETDCombinedScoringResults* get_CIDETDCombinedScoringResults() 
		{
			return mobj_cid_etd_combined_scoring_results ; 
		}
		__property void set_SearchFasta(bool value)
		{
			mbln_search_fasta = value ; 
		}
		__property bool get_SearchFasta()
		{
			return mbln_search_fasta ; 
		}
	};
}