//---------------------------------------------------------------------------


/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

// modified by Anoop for GlypID




#ifndef SequenceH
#define SequenceH

#include <string>
#include <map>
#include <vector>
#include <sstream>
#include "aminoacid.h"

using namespace std;

namespace Engine
{
	namespace SequenceManager
	{
		class Sequence 
		{
		public:
			Sequence (const string& AInfo, const string& ASeqString,const string &site): seqString(ASeqString)
			{
				setInfo(AInfo)		; 
				setNSitePosition(site) ; 
				seqMass = 0 ; 
			}
			Sequence(const string& AInfo, const string& ASeqString): seqString(ASeqString) 
			{
				setInfo(AInfo);
				seqMass = 0 ; 
			}
			Sequence(void);
					
			virtual ~Sequence() { }
  
			const string& getSeqString() const 
			{
				return seqString;
			}

			const string& getName() const 
			{
				return name;
			}

			const string& getInfo() const 
			{
				return info;
			}

			const string& getDescription() const
			{
				return description;
			}

			const string& getSitePosition() const
			{
				return site_position ; 
			}

			map<string, vector<string> >& getAnnotation()
			{
				return annotation;
			}
			void setSeqString(string seq);
			void setSeqName(string seqname) ;
			float calculateMass(bool need_accurate_mass); 

			float seqMass ; 

			void setInfo(const string& AInfo);
			void setNSitePosition(const string& site) ; 

			bool is_decoy ; 
			

		private:
			string seqString;
			string name;
			string info;
			string description;
			string site_position; 
			map<string, vector<string> > annotation;
			
			
			
		};
	}
}

//---------------------------------------------------------------------------
#endif
