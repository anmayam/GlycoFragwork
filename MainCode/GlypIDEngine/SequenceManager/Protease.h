//---------------------------------------------------------------------------

#ifndef ProteaseH
#define ProteaseH

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
#include <set>
#include <iostream>

/** The protease class stores parameters needed by Digest to digest a protein sequence.
* A custom protease can be created or one derived from the attributes set in the
* ProteaseManager.xml resource.
* @author Michael Jones
* @author Mark Schreiber (refactoring to ProteaseManager)
*/

// Modified by Anoop for GlypId

using namespace std;

namespace Engine
{
 
	namespace SequenceManager
	{
		class Protease 
		{
			public:
			  virtual ~Protease(void)
			  { };

			  static set<string> getProteaseList();

			  /**
			  * Retrieves a reference to the named Protease.
			  * (Internally calls ProteaseManager.getProteaseByName())
			  * @param proteaseName A protease name that is registered in the ProteaseManager (case sensitive)
			  * @return A Protease instance for the given protease name
			  */
			  static Protease& getProteaseByName(string proteaseName);

			  friend class ProteaseManager;
			  
			  friend ostream& operator << (ostream& stream, Protease& protease);

			  /*
			  public static final String TRYPSIN = ProteaseManager.TRYPSIN;
			  public static final String LYS_C = ProteaseManager.LYS_C;
			  public static final String ARG_C = ProteaseManager.ARG_C;
			  public static final String ASP_N = ProteaseManager.ASP_N;
			  public static final String GLU_C_BICARB = ProteaseManager.GLU_C_BICARB;
			  public static final String GLU_C_PHOS = ProteaseManager.GLU_C_PHOS;
			  public static final String CHYMOTRYP = ProteaseManager.CHYMOTRYP;
			  public static final String CNBr = ProteaseManager.CNBr;
			  */

			  /**
			  * The list of residues that the protease will cleave at.
			  * @return the residues as a SymbolList
			  */
			  const string getCleaveageResidues()	const	{
				return cleavageResidues;
			  }

			  /**
			  * Gets the name of this Protease
			  * @return the name as a String
			  */
			  const string getName() const {
				return name;
			  }

			  /**
			  * The list of residues that will prevent cleavage if they follow the cleavage
			  * residue.
			  */
			  const string getNotCleaveResidues() const {
				return notCleaveResidues;
				}

			  bool isEndoProtease()	{
				return endoProtease;
			  }

			private:
			  string cleavageResidues;
			  string notCleaveResidues;
			  bool endoProtease;
			  string name;

			protected:
			  Protease(string cleaveRes,
				bool endoProtease,
				string notCleaveRes,
				string name) {
				  this->cleavageResidues = cleaveRes;
				  this->endoProtease = endoProtease;
				  this->notCleaveResidues = notCleaveRes;
				  this->name = name;
			  }
		};

	}
}


#endif
