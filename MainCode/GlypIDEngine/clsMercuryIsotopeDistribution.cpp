// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#pragma once
#include "clsMercuryIsotopeDistribution.h"
#include "GlypIDEngineUtils.h"

namespace GlypID
{
	clsMercuryIsotopeDistribution::clsMercuryIsotopeDistribution()
	{
		mMercuryIsotopeDistribution = new Engine::HornTransformTheoreticalProfile::MercuryIsotopeDistribution();
		set_MercurySize(8192);
		set_ChargeCarrierMass(1.00727638); 		
		set_Resolution(100000);
		set_ChargeState(1);
		
	}

	clsMercuryIsotopeDistribution::~clsMercuryIsotopeDistribution()
	{
		if (mMercuryIsotopeDistribution != NULL)
		{
			delete mMercuryIsotopeDistribution;
			mMercuryIsotopeDistribution = NULL ; 
		}
	}


	System::Drawing::PointF clsMercuryIsotopeDistribution::CalculateDistribution(System::Collections::Hashtable* elementCounts)[]
	{
		std::vector<double> vect_x;
		std::vector<double> vect_y;

		std::vector<double> vect_isotope_x;
		std::vector<double> vect_isotope_y;

		System::Drawing::PointF points __gc [];

		Engine::HornTransformTheoreticalProfile::MolecularFormula formula ; 
		GlypID::Utils::ConvertElementTableToFormula(mMercuryIsotopeDistribution->mobj_elemental_isotope_composition, elementCounts, formula);

		mMercuryIsotopeDistribution->CalculateDistribution(mChargeState, mResolution, formula, vect_x, vect_y, 0, vect_isotope_x, vect_isotope_y, false);

		int num_pts = (int) vect_x.size();
		points = new System::Drawing::PointF __gc [num_pts] ;

		for (int i = 0 ; i < num_pts ; i++)
		{
			points[i] = System::Drawing::PointF((float) vect_x[i], (float) vect_y[i]);
		}

		return points;
	};

	void clsMercuryIsotopeDistribution::CalculateMasses(System::Collections::Hashtable* elementCounts) 
	{		
		Engine::HornTransformTheoreticalProfile::MolecularFormula formula ; 
		GlypID::Utils::ConvertElementTableToFormula(mMercuryIsotopeDistribution->mobj_elemental_isotope_composition, elementCounts, formula);

		mMercuryIsotopeDistribution->CalculateMasses(formula);
	}
}