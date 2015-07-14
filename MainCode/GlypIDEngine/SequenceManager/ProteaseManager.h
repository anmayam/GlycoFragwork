//---------------------------------------------------------------------------
#ifndef ProteaseManagerH
#define ProteaseManagerH

#include <set>
#include <map>
#include <string>
#include "Protease.h"
#include "BioException.h"
using namespace std;

namespace Engine
{
	namespace SequenceManager
	{
		class ProteaseManager
		{
			public:
			  static set<string> getNames();

			  static Protease& getProteaseByName(string proteaseName);

			  static Protease& createProtease(string cleaveRes, bool endoProtease, string notCleaveRes, string name);

			  static bool registered(string proteaseName);
			private:
			  ProteaseManager(void)
			  { };

			  virtual ~ProteaseManager(void)
			  { };

			  static map<string, Protease*> name2Protease;

			  static void registerProtease(Protease* prot);
			};

	}
}
#endif
