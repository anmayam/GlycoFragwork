//---------------------------------------------------------------------------
/// From Tiger's proteomics library
#ifndef DigestH
#define DigestH

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

/**
 * This class contains methods for calculating the results of proteolytic digestion
 * of a protein sequence
 *
 * <b> this class is not designed to be thread safe </b>
 *
 * @author Michael Jones
 * @author Mark Schreiber (refactoring, some documentation)
 */

#include <string>
#include "Protease.h"
#include "Sequence.h"
#include "RangeLocation.h"
#include "RangeLocationFilter.h"

using namespace std;

namespace Engine
{
	namespace SequenceManager
	{
		class Digest 
		{
			public:
				static string PEPTIDE_FEATURE_TYPE;
				
				Digest()
				{
					protease = NULL;
					sequence = NULL;
					maxMissedCleavages= 0;
					peptideQue = NULL;
					filter = NULL;
				}

				virtual ~Digest()
				{
					if(NULL != peptideQue)
					{
					  delete peptideQue;
					}
				}
				
				void setProtease(Protease* protease)
				{
					this->protease = protease;
				}

				void setSequence(Sequence* sequence)
				{
					this->sequence = sequence;
				}

				Sequence* getSequence()
				{
					return sequence;
				}

			  /**
			   * Sets the maximum number of partial digest products to be annotated.
			   * @param maxMissedCleavages the max number of partial digest products
			   */
			   void setMaxMissedCleavages(int maxMissedCleavages) 
			   {
					this->maxMissedCleavages = maxMissedCleavages;
			   }

			  void setRangeLocationFilter(IRangeLocationFilter* filter)
			  {
					this->filter = filter;
			  }

			  /** Adds peptides as features to the Sequence in this class. The feature will
			   * contain a small annotation specifying the protease with the key "protease".
			   * For Example:
			   * <PRE>
			   *
			   *         Sequence sequence = ...
			   *         Digest bioJavaDigest = new Digest();
			   *
			   *         bioJavaDigest.setMaxMissedCleavages(2);
			   *         bioJavaDigest.setProtease(ProteaseManager.getProteaseByName(Protease.ASP_N));
			   *         bioJavaDigest.setSequence(sequence);
			   *         bioJavaDigest.addDigestFeatures();
			   * </PRE>
			   * @throws BioException if the Protease or Sequence are null.
			   */
				  void addDigestFeatures();
		private:
			  Protease* protease;
			  Sequence* sequence;
			  int maxMissedCleavages;
			  vector<RangeLocation>* peptideQue;
			  IRangeLocationFilter* filter;

			  void addMissedCleavages();
			  void createPeptideFeature(const RangeLocation& loc);
			};

	}
}
//---------------------------------------------------------------------------
#endif
