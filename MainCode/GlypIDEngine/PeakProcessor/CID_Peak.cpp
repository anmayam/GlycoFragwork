/*****************************************************************************\
 * CID Peak.cpp
 * 
 *
 * author: Yin Wu (wuyin@indiana.edu) Modified by Anoop
 *
\*****************************************************************************/


#include "CID_Peak.h"
#include "..\Utilities\system.h"

namespace Engine
{
	namespace PeakProcessing
	{		

		
		CID_Peak::~CID_Peak(void)
		{
			
		}

		void CID_Peak::init()
		{
			this->mz = 0;
			this->intensity = 0;
			this->original_intensity = 0;
			this->original_mass = 0;
			this->glycan_mass = 0;
			this->glycan_charge = 1;
			this->ms2_score = 0;
			this->b_ion_in_ms2 = 0;
			this->ms2_scanheader = 0;
			this->charge = 1;
			this->intensity_rank = 0;
			this->flag = 0;
			this->parent_isoform_id = 0;
			
		}	


		CID_Peak::CID_Peak(float mz, float intensity)
		{
			  this->init();
			  this->mz = mz;
			  this->intensity = intensity;
			  this->glycan_mass = 0;
			  this->glycan_charge = 1;
		}

		


		//used only if the charge state of the peak is known
		CID_Peak::CID_Peak(float mass, float intensity, int charge)
		{
			this->init();
			this->mz = mz;
			this->intensity = intensity;
			this->charge = charge;
			this->glycan_mass = 0;
			this->glycan_charge = charge;
			this->intensity_rank = 0;			
		}
		
		CID_Peak::CID_Peak(float mass, float intensity, float glycan_mass, int glycan_cs)
		{
			this->init();
			this->mz = mz;
			this->glycan_mass = glycan_mass;
			this->glycan_charge = glycan_cs;
			this->intensity = intensity;
			this->intensity_rank = 0;
		}


		CID_Peak::CID_Peak(void)
		{
			  this->init();
		}



		float CID_Peak::originalMass()
		{
			//Anoop : be careful when using this . Have not resolved whether glycan_mass is mass or mz

			//get the original mass of the peak in the spectrum (after the transformation)
			//it depends on the current mass, glycan_mass, 
			//and glycan_charge
			if ( glycan_charge > 0 )
					return (mz + glycan_mass - H2O + PROTON * glycan_charge)/glycan_charge;
			else
				return 0;
		}


		//an empty peak is like a NULL value of pointer
		bool CID_Peak::isEmpty()
		{
			if ( intensity == 0 && mz == 0 )
				return true;
			else
				return false;
		}

		//set a peak to be empty
		void CID_Peak::setEmpty() 
		{
			this->init();
			this->intensity = 0;
			this->mz= 0;
		}		
	}
}
	