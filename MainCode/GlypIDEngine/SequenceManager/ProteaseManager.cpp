#include "ProteaseManager.h"
namespace Engine
{
	namespace SequenceManager
	{

		map<string, Protease*> ProteaseManager::name2Protease;

		set<string>
		ProteaseManager::
		getNames() {
		  set<string> result;
		  for(map<string, Protease*>::const_iterator iter = name2Protease.begin();
			iter != name2Protease.end(); iter++){
			result.insert((*iter).first);
		  }
		  return result;
		}

		Protease&
		ProteaseManager::createProtease(
			string cleaveRes,
			bool endoProtease,
			string notCleaveRes,
			string name) {
		  Protease* p = new Protease(cleaveRes, endoProtease, notCleaveRes, name);
		  registerProtease(p);
		  return *p;
		}

		bool
		ProteaseManager::registered(string proteaseName){
		  return name2Protease.find(proteaseName) != name2Protease.end();
		}

		void
		ProteaseManager::registerProtease(Protease* prot) {
		  if(registered(prot->getName()))
		  {
			throw BioException(string("A Protease has already been registered with the name ")+ prot->getName());
		  }

		  name2Protease[prot->getName()] = prot;
		}

		Protease&
		ProteaseManager::getProteaseByName(string proteaseName) {
		  Protease* protease = name2Protease[proteaseName];
		  if(NULL == protease)
		  {
		    throw BioException("No protease has been registered by that name");
		  }
		  return *protease;
		}

	}
}
