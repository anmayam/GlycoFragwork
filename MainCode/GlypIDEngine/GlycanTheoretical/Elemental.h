// Written by Yin Wu (wuyin@indiana.edu)
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// modified by Anoop for GlypID
// Written by Yin Wu (wuyin@indiana.edu)
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// modified by Anoop for GlypID


#ifndef ELEMENTAL_H
#define ELEMENTAL_H

#define EPSILON 1.0e-16
#define PROTON 1.00727638

#include <math.h>

#include <string>
#include <map>
#include <vector>
#include <iostream>
using namespace std;

namespace Engine
{
	namespace GlycanTheoretical
	{

		class Elemental 
		{
			private:
				  string name;
				  map<int, double> isotopes; //stable isotopes stored in pair <mass, faction>

				  void init();
				  void trimEnvelope(vector<double>& env);
				  string elementListToFormula(map<string, int>& elements);
			public:
			  Elemental(string& name);
			  string &getName();
			  bool isEqualTo(Elemental &ele);
			  map<int, double>& getIsotopes();
			  bool isEmpty();
			  
			  map<string, int> readMoleculeFormula(string& formula);
			  vector<double> getSingleElementEnvelope(Elemental& ele, int quanitity);
			  vector<double> getSingleElementEnvelope3(Elemental& ele, int quanitity);

			  vector<double> fuseIsotopicEnvelopes(const vector<double>& env1, const vector<double>& env2);

			  vector<double> getIsotopicEnvelope(string& formula);

			  string seqToFormula(string& seq);
			  string glycanToFormula(string& glycan);
			  string glycoPepToFormula(string& seq);
			  double formulaToMass(string& formula);

			  int highestIsotopeOfFormula(int k, string &formula);
		};
	}
}


#endif
