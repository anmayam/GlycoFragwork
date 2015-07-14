//---------------------------------------------------------------------------


#pragma hdrstop

#include "RangeLocationFilter.h"
#include "Aminoacid.h"


namespace Engine
{
	namespace SequenceManager
	{
		void NGlycanFilter::setSequence(Sequence* sequence) 
		{
			  delete isNglycan;
			  
			  isNglycan = new bool[sequence->getSeqString().size()];
			  
			  for(size_t i = 0;i < sequence->getSeqString().size();i++){
				isNglycan[i] = false;
			  }

			  for(size_t i = 0;i < sequence->getSeqString().size() - 2;i++){
				if(sequence->getSeqString()[i] == 'N'){
				  if('P' == sequence->getSeqString()[i+1] ) {
					 //'D' == sequence->getSeqString()[i+1]){
					continue;
				  }

				  if('S' == sequence->getSeqString()[i+2] ||
					 'T' == sequence->getSeqString()[i+2]){
					isNglycan[i] = true;					
				  }
				}
			  }
		}



		bool NGlycanFilter::accept(RangeLocation& rl)
		{
			  for(int index = rl.getMin() - 1; index <= rl.getMax() - 1;index ++)
			  {
				if(isNglycan[index])
				{
				  return true;
				}
			  }
			  return false;
		}

		bool PeptideWeightFilter::accept(RangeLocation& rl) 
		{
		  string peptide = sequence.substr(rl.getMin() - 1, rl.getMax() - rl.getMin() + 1);

		  double mass = aminoacids.averagePeptideMass(peptide);

		  return mass >= minWeight && mass <= maxWeight;
		}

	}
	
}

//---------------------------------------------------------------------------

#pragma package(smart_init)
