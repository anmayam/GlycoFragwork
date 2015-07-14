//---------------------------------------------------------------------------


#pragma hdrstop

#include "Digest.h"
#include "BioException.h"

namespace Engine
{
	namespace SequenceManager
	{

		string Digest::PEPTIDE_FEATURE_TYPE = "Peptide";

		void Digest::addDigestFeatures() 
		{
		  if(NULL == protease)
		  {
			  throw BioException("Protease is null, use Digest.setProtease()");
		  }
		  
		  if(NULL == sequence)
		  {
			  throw BioException("Sequence is null, use Digest.setSequence()");
		  }

		  string cleaveSites = protease->getCleaveageResidues();
		  string notCleave = protease->getNotCleaveResidues();

		  if(cleaveSites.empty() && notCleave.empty())
		  {
			throw BioException("Protease contains no cleave information");
		  }

		  delete peptideQue;
		
		  peptideQue = new vector<RangeLocation>();

		  bool endoProtease = protease->isEndoProtease();
		  size_t nTerm = 1;

		  const string& seq = sequence->getSeqString();
		  for (size_t j = 1; j <= seq.size(); j++) {
			char aa = seq[j-1];

			if(cleaveSites.find(aa, 0) != string::npos){
			  if (endoProtease) {
				bool cleave = true;
				if (j < seq.size())  {
				  char nextAA = seq[j];
				  if(notCleave.find(nextAA, 0) != string::npos){
					cleave = false;
				  }
				}

				if (cleave)  {
				  peptideQue->push_back(RangeLocation(nTerm, j));
				  nTerm = j + 1;
				}
			  } else {
				if (j > 1) {
				  peptideQue->push_back(RangeLocation(nTerm, j-1));
				  nTerm = j;
				}
			  }
			}
		  }

		  if (nTerm <= seq.size()) {
			peptideQue->push_back(RangeLocation(nTerm, seq.size()));
		  }

		  addMissedCleavages();

			//Now add the locations as Peptide freatures to the Sequence
		  if(filter != NULL){
			filter->setSequence(sequence);
			for(vector<RangeLocation>::iterator iter = peptideQue->begin(),
			  iterEnd = peptideQue->end(); iter != iterEnd; iter++) {
			  if(filter->accept(*iter)){
				createPeptideFeature(*iter);
			  }
			}
		  }
		  else {
			for(vector<RangeLocation>::iterator iter = peptideQue->begin(),
			  iterEnd = peptideQue->end(); iter != iterEnd; iter++) {
			  createPeptideFeature(*iter);
			}
		  }
		}

		void Digest::addMissedCleavages()
		{
		  if(maxMissedCleavages > 0){
			vector<RangeLocation> missedList;
			for(size_t i = 0;i < peptideQue->size();i++){
			  const RangeLocation& loc = (*peptideQue)[i];

			  for(int iMiss = 0;iMiss < maxMissedCleavages;iMiss++){
				size_t j = i + 1 + iMiss;
				if(j == peptideQue->size()){
				  break;
				}

				const RangeLocation& loc2 = (*peptideQue)[j];
				missedList.push_back(RangeLocation(loc.getMin(), loc2.getMax()));
			  }
			}

			peptideQue->insert(peptideQue->end(), missedList.begin(), missedList.end());
		  }
		}

		void Digest::createPeptideFeature(const RangeLocation& loc) 
		{
		  vector<string>& peptides = sequence->getAnnotation()[PEPTIDE_FEATURE_TYPE];
		  string peptide = sequence->getSeqString().substr(loc.getMin() - 1, loc.getMax() - loc.getMin() + 1);
		  peptides.push_back(peptide);
		}
	}
}

//---------------------------------------------------------------------------

#pragma package(smart_init)
