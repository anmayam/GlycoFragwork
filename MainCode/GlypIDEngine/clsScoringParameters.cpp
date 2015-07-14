// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// All CID and HCD scoring parameters go here



#include ".\clsscoringparameters.h"
#using <mscorlib.dll>

#using <System.Xml.dll>
namespace GlypID
{
	namespace Scoring
	{
		clsScoringParameters::clsScoringParameters(void)
		{
			
			mint_min_num_peaks_to_consider = 1; 
			mint_min_path_length = 2 ; 
			mdbl_min_hcd_mz = 0 ; 
			mdbl_max_hcd_mz = 700 ; 
			mdbl_min_cid_mz = 100 ; 
			mdbl_max_cid_mz = 2000 ; 
			mdbl_ppm_tolerance = 20 ; 
			mdbl_Da_tolerance = 0.5 ; 
			mbln_use_ppm = true ; 
			mbln_allocate_default_charges = true ; 
			mbln_filter_comb_list = true ; 
			mbln_include_best_hit = true ; 
			mbln_use_y1_for_id = false ; 
		}

		clsScoringParameters::~clsScoringParameters(void)
		{
		}

		void clsScoringParameters::SaveScoringParameters(System::Xml::XmlTextWriter *xwriter)
		{
			xwriter->WriteStartElement(S"ScoringParameters") ; 
			xwriter->WriteWhitespace(S"\n\t\t") ; 

			xwriter->WriteElementString(S"MinHCDMz", System::Convert::ToString(this->mdbl_min_hcd_mz)); 
			xwriter->WriteWhitespace(S"\n\t\t") ; 
			xwriter->WriteElementString(S"MaxHCDMz", System::Convert::ToString(this->mdbl_max_hcd_mz)); 
			xwriter->WriteWhitespace(S"\n\t\t") ; 
			xwriter->WriteElementString(S"MinNumHCDPeaks", System::Convert::ToString(this->mint_min_num_peaks_to_consider)); 
			xwriter->WriteWhitespace(S"\n\t\t") ; 
			xwriter->WriteElementString(S"MinCIDMz", System::Convert::ToString(this->mdbl_min_cid_mz)); 
			xwriter->WriteWhitespace(S"\n\t\t") ; 
			xwriter->WriteElementString(S"MaxCIDMz", System::Convert::ToString(this->mdbl_max_cid_mz)); 
			xwriter->WriteWhitespace(S"\n\t\t") ;
			xwriter->WriteElementString(S"MinPathLength", System::Convert::ToString(this->mint_min_path_length)); 
			xwriter->WriteWhitespace(S"\n\t\t") ; 					
			xwriter->WriteElementString(S"UsePPMForSearchMassTolerance", System::Convert::ToString(this->mbln_use_ppm)); 
			xwriter->WriteWhitespace(S"\n\t\t") ; 		
			xwriter->WriteElementString(S"MinSearchTolInDa", System::Convert::ToString(this->mdbl_Da_tolerance)); 
			xwriter->WriteWhitespace(S"\n\t\t") ; 		
			xwriter->WriteElementString(S"MinSearchTolInPPM", System::Convert::ToString(this->mdbl_ppm_tolerance)); 
			xwriter->WriteWhitespace(S"\n\t\t") ; 	
			xwriter->WriteElementString(S"AllocateDefaultCharges", System::Convert::ToString(this->mbln_allocate_default_charges)); 
			xwriter->WriteWhitespace(S"\n\t\t") ; 	
			xwriter->WriteElementString(S"FilterCombinationList", System::Convert::ToString(this->mbln_filter_comb_list)); 
			xwriter->WriteWhitespace(S"\n\t\t") ; 	
			xwriter->WriteElementString(S"IncludeBestHit", System::Convert::ToString(this->mbln_include_best_hit)); 
			xwriter->WriteWhitespace(S"\n\t\t") ;
			xwriter->WriteEndElement();
			xwriter->WriteWhitespace(S"\n\t") ; 

		}

		void clsScoringParameters::LoadScoringParameters(XmlReader *rdr)
		{
			while (rdr->Read())
			{
				switch (rdr->NodeType)
				{
					case XmlNodeType::Element:
						if (rdr->Name->Equals(S"MinHCDMz"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for MinHCDMZ in parameter file") ; 
							}
							this->set_MinHCDMz(System::Convert::ToDouble(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"MaxHCDMz"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for MaxHCDMZ in parameter file") ; 
							}
							this->set_MaxHCDMz(System::Convert::ToDouble(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"MinNumHCDPeaks"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for MinNumHCDPeaks in parameter file") ; 
							}
							this->set_MinNumPeaksToConsider(System::Convert::ToInt32(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"MinCIDMz"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for MinCIDMz in parameter file") ; 
							}
							this->set_MinCIDMz(System::Convert::ToDouble(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"MaxCIDMz"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for MaxCIDMz in parameter file") ; 
							}
							this->set_MaxCIDMz(System::Convert::ToDouble(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"MinPathLength"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for MinPathLength in parameter file") ; 
							}
							this->set_MinPathLength(System::Convert::ToInt32(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"UsePPMForSearchMassTol"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for MinSearchTolInPPM in parameter file") ; 
							}
							this->set_UsePPM(System::Convert::ToBoolean(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"IncludeBestHit"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for Report All Matches in parameter file") ; 
							}
							this->set_IncludeBestHit(System::Convert::ToBoolean(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"MinSearchTolInDa"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for MinSearchTolInPPM in parameter file") ; 
							}
							this->set_DaTolerance(System::Convert::ToDouble(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"MinSearchTolInPPM"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for MinSearchTolInPPM in parameter file") ; 
							}
							this->set_PPMTolerance(System::Convert::ToDouble(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"AllocateDefaultCharges"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for AllocateDefaultCharges in parameter file") ; 
							}
							this->set_AllocateDefaultCharges(System::Convert::ToBoolean(rdr->Value)) ; 
						}
						else if (rdr->Name->Equals(S"FilterCombinationList"))
						{
							rdr->Read() ; 
							while(rdr->NodeType == XmlNodeType::Whitespace || rdr->NodeType == XmlNodeType::SignificantWhitespace)
							{
								rdr->Read() ; 
							}
							if (rdr->NodeType != XmlNodeType::Text)
							{
								throw new System::Exception (S"Missing information for FilterCombinationList in parameter file") ; 
							}
							this->set_FilterCombList(System::Convert::ToBoolean(rdr->Value)) ; 
						}
						break ; 
					case XmlNodeType::EndElement:
						if (rdr->Name->Equals(S"ScoringParameters"))
							return ;
						break ; 
					default:
						break ; 
				}
			}
		}
	}
}
