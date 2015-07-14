//---------------------------------------------------------------------------

// Original code from Tiger's proteomics omics code
// Modified by Anoop for GlypID

#pragma hdrstop

#include "Sequence.h"
#include "../Utilities/system.h"

namespace Engine
{
	namespace SequenceManager
	{
		Sequence::Sequence(void)
		{
			site_position.clear();
			is_decoy = false ; 
			
		}
		void Sequence::setInfo(const string& AInfo)
		{
			this->info = AInfo;
			size_t ipos = AInfo.find(' ');
			if (ipos == string::npos)
			{
				ipos = AInfo.find('\t');
			}
			if (ipos == string::npos)
			{
				this->name = AInfo;
				this->description = "";
			}
			else
			{
				this->name = AInfo.substr(0,ipos);
				this->description = AInfo.substr(ipos, AInfo.length());
			}
		}
		
		
		float Sequence::calculateMass(bool need_accurate_mass)
		{
			//add modification when calculating peptide mass
			TAminoacids tas;
			tas['C'].resetMass(tas['C'].monoMass() + IODOACEDAMIDE_MONO, tas['C'].averageMass() + IODOACEDAMIDE_AVG); 
			
			string str = this->getSeqString();
			float pep_mass;
			if ( need_accurate_mass )
				pep_mass = (float) tas.monoPeptideMass(str);
			else
				pep_mass = (float) tas.averagePeptideMass(str);			
			return pep_mass ; 
		}

		void Sequence::setSeqString(string seq)
		{
			seqString = seq; 
		}

		void Sequence::setSeqName(string seqname)
		{
			name = seqname ; 
		}

		void Sequence::setNSitePosition(const string &site)
		{
			this->site_position = site ; 		
			
		}

		
	}
}
//---------------------------------------------------------------------------

//Anoop commenting this out as I have no clue what this does
//#pragma package(smart_init)
