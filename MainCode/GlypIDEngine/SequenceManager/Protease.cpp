#include "Protease.h"
#include "ProteaseManager.h"

namespace Engine
{
	namespace SequenceManager
	{

		ostream& operator << (ostream& stream,Protease& protease) {
		  stream << protease.name;

		  if(protease.cleavageResidues.empty()){
			stream << "\t-";
		  }
		  else{
			stream << "\t" + protease.cleavageResidues;
		  }

		  if(protease.notCleaveResidues.empty()){
			stream << "\t-";
		  }
		  else{
			stream << "\t" << protease.notCleaveResidues;
		  }

		  if(protease.endoProtease){
			stream << "\t1";
		  }
		  else {
			stream << "\t0";
		  }

		  return stream;
		}

		set<string> Protease::getProteaseList() {
			return ProteaseManager::getNames();
		}

		Protease& Protease::getProteaseByName(string proteaseName)
		{
			return ProteaseManager::getProteaseByName(proteaseName);
		}		
	}
}
