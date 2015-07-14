//---------------------------------------------------------------------------

#ifndef RangeLocationH
#define RangeLocationH

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

#include <string>
#include <sstream>
#include "BioException.h"

using namespace std;


namespace Engine
{
	namespace SequenceManager
	{
		/**
		 * A simple implementation of Location that contains all points between
		 * getMin and getMax inclusive.
		 * <p>
		 * This will in practice be the most commonly used pure-java implementation.
		 *
		 * @author Matthew Pocock
		 */
		class RangeLocation 
		{
				public:
				  RangeLocation(int min, int max) 
				  {
						if(max < min)
						{
						  ostringstream out;
						  out << "max must exceed min: min=" << min << ", max=" << max;
						  throw BioException(out.str());
						}
						this->min = min;
						this->max = max;
				  }


				  int getMin() const 
				  {
					return min;
				  }

				  int getMax() const
				  {
					return max;
				  }

				  RangeLocation translate(int dist)
				  {
					return RangeLocation(min + dist, max + dist);
				  }

				private:
				  int min;
				  int max;
			};
	}
}
//---------------------------------------------------------------------------
#endif
