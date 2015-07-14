// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#pragma once
#using <mscorlib.dll>
#using <System.Xml.dll>

using namespace System;
#include "HornTransformTheoreticalProfile/AtomicInformation.h"


namespace GlypID
{
public __gc class clsElementIsotopes: public ICloneable

	{
	public:
		Engine::HornTransformTheoreticalProfile::AtomicInformation __nogc *mobjAtomicInfo ; 
		clsElementIsotopes() ; 
		

		virtual Object* Clone() ;
		virtual clsElementIsotopes* Assign(clsElementIsotopes* otherOne) ; 

		~clsElementIsotopes() ; 		
		void GetElementalIsotope(int index, Int32& atomicity, Int32& num_isotopes, 
			System::String* (&element_name), System::String* (&element_symbol),
			float __gc& average_mass, float __gc& mass_variance, float (&isotope_mass) __gc[], 
			float (&isotope_prob) __gc[]) ; 
		void UpdateElementalIsotope(int index, Int32& atomicity, Int32& isotope_num, 
			System::String* (&element_name), System::String* (&element_symbol),
			double  __gc& isotope_mass, double __gc &isotope_prob) ;
		int GetNumberOfElements();
		void SetElementalIsotopeComposition(Engine::HornTransformTheoreticalProfile::AtomicInformation 
			__nogc *atomic_info) ; 
		const Engine::HornTransformTheoreticalProfile::AtomicInformation* GetElementalIsotopeComposition() ; 
		void LoadElementIsotopes(System::Xml::XmlReader *rdr) ; 
		void SaveElementIsotopes(System::Xml::XmlTextWriter *xwriter) ; 
	};	
}
