#include "../PeakProcessor/PeakProcessor.h"
#include "../PeakProcessor/PeakData.h"
#include "./CIDInformationRecord.h"
#include "./Spectrum.h" 
#include <vector>



namespace Engine
{
	namespace MS2CIDScoring
	{
		class CIDScoring
		{
			private:						
				
				float p_value_cutoff;
				float  ms2_score_threshold; //the threshold above which, the ion is considered as glycopeptide
				Engine::MS2CIDScoring::Spectrum mobj_Spectrum ; 
				double mdbl_GlcNac_Mass ;  
				// charge carrier mass
				double mdbl_cc_mass; 
				
		
			public:				
				CIDScoring(); 
				~CIDScoring() ; 				
				double CalculateCIDScore(Engine::PeakProcessing::PeakData &pk_data, int charge) ; 
				bool CheckForOxoniumIons(Engine::PeakProcessing::PeakData &pk_data) ; 
				double DetermineY1IonThroughPeptideSearchingCID(Engine::PeakProcessing::PeakData &pk_data, Engine::MS2CIDScoring::CIDInformationRecord &cid_record) ; 
				double CalculateCIDScorePValue(int score) ; 
				void InitPTable(); 				



		}; 
	}
}
