// Written by Anoop Mayampurath, Indiana Univeristy
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0
#pragma once


#using <System.Xml.dll>
using namespace System::Xml ; 

namespace GlypID
{
	namespace Scoring
	{

		public __gc class clsScoringParameters: public System::ICloneable
		{
			static System::String *DEFAULT_ISOTOPE_FILE = S"isotope.xml" ; 
			bool mbln_use_scan_range ; 
			int mint_min_scan ; 
			int mint_max_scan ; 
			double mdbl_min_hcd_mz ; 
			double mdbl_max_hcd_mz ; 
			double mdbl_min_cid_mz ; 
			double mdbl_max_cid_mz ; 
			int mint_min_num_peaks_to_consider ; 
			int mint_min_path_length ; 
			double mdbl_ppm_tolerance ; 
			double mdbl_Da_tolerance ; 
			bool mbln_use_ppm ; 
			bool mbln_allocate_default_charges ; 
			bool mbln_filter_comb_list ; 
			bool mbln_include_best_hit ; 
			bool mbln_use_y1_for_id ; 
			
		public:
			clsScoringParameters(void);
			~clsScoringParameters(void);
			void SaveScoringParameters(System::Xml::XmlTextWriter *xwriter) ; 
			void LoadScoringParameters(System::Xml::XmlReader *rdr) ; 			

			virtual Object* Clone()
			{
				clsScoringParameters *new_params = new clsScoringParameters() ; 
				//Copy stuff
				new_params->set_MinHCDMz(this->get_MinHCDMz()) ; 
				new_params->set_MaxHCDMz(this->get_MaxHCDMz()) ; 				
				new_params->set_MinNumPeaksToConsider(this->get_MinNumPeaksToConsider()) ; 
				new_params->set_AllocateDefaultCharges(this->get_AllocateDefaultCharges()) ; 
				new_params->set_DaTolerance(this->get_DaTolerance()) ; 
				new_params->set_FilterCombList(this->get_FilterCombList()) ; 
				new_params->set_MinPathLength(this->get_MinPathLength()) ; 
				new_params->set_PPMTolerance(this->get_PPMTolerance()); 
				new_params->set_UsePPM(this->get_UsePPM()) ; 
				new_params->set_IncludeBestHit(this->get_IncludeBestHit()) ; 
				return new_params ; 
			}

			__property void set_UseY1ForIdentification (bool val)
			{
				mbln_use_y1_for_id = val;
			}

			__property bool get_UseY1ForIdentification()
			{
				return mbln_use_y1_for_id ; 
			}

			__property void set_IncludeBestHit (bool val)
			{
				mbln_include_best_hit = val; 
			}

			__property bool get_IncludeBestHit()
			{
				return mbln_include_best_hit ; 
			}

			__property void set_MinNumPeaksToConsider(int val)
			{
				mint_min_num_peaks_to_consider = val ; 
			}

			__property int get_MinNumPeaksToConsider()
			{
				return mint_min_num_peaks_to_consider; 
			}
			__property void set_MinPathLength(int val)
			{
				mint_min_path_length = val ; 
			}
			__property int get_MinPathLength()
			{
				return mint_min_path_length ; 
			}
			__property void set_PPMTolerance (double val)
			{
				mdbl_ppm_tolerance = val ; 
			}
			__property double get_PPMTolerance()
			{
				return mdbl_ppm_tolerance ; 
			}

			__property void set_DaTolerance (double val)
			{
				mdbl_Da_tolerance = val ; 
			}
			__property double get_DaTolerance()
			{
				return mdbl_Da_tolerance ; 
			}

			__property void set_UsePPM (bool val)
			{
				mbln_use_ppm = val ; 
			}
			__property bool get_UsePPM()
			{
				return mbln_use_ppm ; 
			}	
			__property void set_AllocateDefaultCharges (bool val)
			{
				mbln_allocate_default_charges = val ; 
			}
			__property bool get_AllocateDefaultCharges()
			{
				return mbln_allocate_default_charges; 
			}
			__property void set_FilterCombList(bool val)
			{
				mbln_filter_comb_list= val ; 
			}
			__property bool get_FilterCombList()
			{
				return mbln_filter_comb_list; 
			}

			__property void set_MinHCDMz(double mz)
			{
				mdbl_min_hcd_mz = mz ; 
			}
			__property void set_MaxHCDMz(double mz)
			{
				mdbl_max_hcd_mz = mz ; 
			}
			__property double get_MaxHCDMz()
			{
				return mdbl_max_hcd_mz ; 
			}
			__property double get_MinHCDMz()
			{
				return mdbl_min_hcd_mz ; 
			}
			__property void set_MinCIDMz(double mz)
			{
				mdbl_min_cid_mz = mz ; 
			}
			__property void set_MaxCIDMz(double mz)
			{
				mdbl_max_cid_mz = mz ; 
			}
			__property double get_MaxCIDMz()
			{
				return mdbl_max_cid_mz ; 
			}
			__property double get_MinCIDMz()
			{
				return mdbl_min_cid_mz ; 
			}
			
		};
	}
}

