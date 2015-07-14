// Written by Yin Wu (wuyin@indiana.edu)
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// modified by Anoop for GlypID

#include "Elemental.h"
#include "..\Utilities\String_Utils.h"
#include "..\Utilities\util.h"
#include "..\Utilities\system.h"
#include <assert.h>
#include <cctype>
using namespace std;

namespace Engine
{
	namespace GlycanTheoretical
	{		
		Elemental::Elemental(string& name) 
		{
		  this->init();
		  if ( name == "C" || name == "c" ) 
		  {
			this->name = "C";
			this->isotopes[12] = 0.988922;
			this->isotopes[13] = 0.011078;
		  }
		  else if ( name == "H" || name == "h" )
		  {
			this->name = "C";
			this->isotopes[1] = 0.99984426;
			this->isotopes[2] = 0.00015574;
		  }
		  else if ( name == "O" || name == "o" ) 
		  {
			this->name = "O";

			/*
			this->isotopes[16] = 0.997628; //
			this->isotopes[18] = 0.002; //
			*/

			/*
			this->isotopes[16] = 0.9977628; //0.997628; //
			this->isotopes[18] = 0.0022372; //0.002; //
			*/

			this->isotopes[16] = 0.99757;
			this->isotopes[17] = 0.00038;
			this->isotopes[18] = 0.00205;
		  }
		  else if ( name == "N" || name == "n" ) 
		  {
			this->name = "N";
			this->isotopes[14] = 0.996337;
			this->isotopes[15] = 0.003663;
		  }
		  else if ( name == "S" || name == "s" ) 
		  {
			this->name = "S";
			this->isotopes[32] = 0.95018; //0.954015; //averaged with S33 to sum fraction to 1
			this->isotopes[34] = 0.04215; //0.045985;
		  }
		}

		void Elemental::init() 
		{
			  this->name = "";
			  this->isotopes.clear();
		}

		string& Elemental::getName() 
		{ 
			return name;
		}
	
		bool Elemental::isEmpty() 
		{
			return name == "";
		}		



		bool Elemental::isEqualTo(Elemental& ele) 
		{
			if ( name.compare("") != 0 && name.compare(ele.getName()) == 0 )
				return true;
			else 
				return false;
		}

		map<int, double>& Elemental::getIsotopes() 
		{
			return isotopes;
		}


		////transform a glyco peptide peptide + glycan to its elemental formula
		string Elemental::glycoPepToFormula(string& seq) 
		{ 
			string glyco_pep = Engine::Utilities::StringUtils::to_upper_copy(seq);
			string peptide, glycan;

			//split the seq into peptide and glycan
			string::size_type loc = seq.find("(", 0);

			if ( loc == string::npos ) 
			{
				peptide = seq;
				glycan = "";	
			}
			else 
			{	
				peptide = string(seq, 0, loc);
				glycan = string(seq, loc, seq.length() - loc);
			}

			//convert both strings into formulas
			string peptide_formula = seqToFormula(peptide);
			string glycan_formula = glycanToFormula(glycan);
			map<string, int> pep_f = readMoleculeFormula(peptide_formula);
			map<string, int> glycan_f = readMoleculeFormula(glycan_formula);

			//combine the two formulas together
			map<string, int>::iterator iter1;
			for (iter1 = pep_f.begin(); iter1 != pep_f.end(); iter1++) 
			{
				if ( glycan_f.find(iter1->first) !=	 glycan_f.end() ) 
				  glycan_f[iter1->first] += iter1->second;
				else
				  glycan_f[iter1->first] = iter1->second;
			  }

			  //subtract an H2O
			  glycan_f["H"] -= 2;
			  glycan_f["O"] -= 1;

			  string result = elementListToFormula(glycan_f);
			  return result;
		}

		////transform a glycan composition string to its elemental formula
		string Elemental::glycanToFormula(string& glycan) 
		{
		  map<string, string> sugar_to_formula;
		  map<string, string>::iterator iter1;

		  sugar_to_formula["HEX"] = "C6H10O5";
		  sugar_to_formula["HNAC"] = "C8H13N1O5";
		  sugar_to_formula["DEHX"] = "C6H10O4";
		  sugar_to_formula["SACID"] = "C11H17N1O8";

		  map<string, int> mono;
		  map<string, int>::iterator iter2;

		  //find each monosaccharide
		  for (iter1 = sugar_to_formula.begin();  iter1 != sugar_to_formula.end(); iter1++) 
		  {
			string::size_type loc = glycan.find(iter1->first, 0);
			//int loc = (int) glycan.find(iter1->first, 0);

			//read the quantity of the monosaccharide
			if ( loc != string::npos ) \
			{
			  loc = loc + iter1->first.size(); //ignore the name
			  int quantity = 0;
			  while ( loc < glycan.size() && glycan[loc] >= '0' &&  glycan[loc] <= '9' ) 
			  {
				char cp[2];
				cp[0] = glycan[loc]; cp[1] = '\0';
				quantity = quantity * 10 + atoi(cp);
				loc++;
			  }
			  //iter1->second = quantity;
			  mono[iter1->first] = quantity;
			}
		 }


		  //convert the monosaccharide composition to elemental composition
		  map<string, int> f;
		  for (iter2 = mono.begin(); iter2 != mono.end(); iter2++) 
		  {
			int sugar_quantity = iter2->second;
			string elemental_formula = sugar_to_formula[iter2->first];
			map<string, int> temp_f = readMoleculeFormula(elemental_formula);
			map<string, int>::iterator iter3;

			for (iter3 = temp_f.begin(); iter3 != temp_f.end(); iter3++) 
			{
			  if ( f.find(iter3->first) != f.end() )
					f[iter3->first] += iter3->second * sugar_quantity;
			  else
					f[iter3->first] = iter3->second * sugar_quantity;
			}
		  }

		  //add an H2O and and a few protons
		  f["H"] += 2; 
		  f["O"] += 1;

		  string result = elementListToFormula(f);

		  /*
		  string result = "";
		  char line[256];
		  for (iter2 = f.begin(); 
			   iter2 != f.end(); iter2++) {
			sprintf(line, "%s%d", iter2->first.c_str(), iter2->second);
			result.append(string(line));
		  }
		  */

		  return result;  
		}

		//convert the elements to a string
		string Elemental::elementListToFormula(map<string, int>& elements) 
		{
			  map<string, int>::iterator iter2;
			  string result = "";
			  char line[256];
			  for (iter2 = elements.begin(); 
				   iter2 != elements.end(); iter2++) {
				sprintf(line, "%s%d", iter2->first.c_str(), iter2->second);
				result.append(string(line));
			  }

			  return result;
		}


		string Elemental::seqToFormula(string& seq) 
		{
		  map<string, string> aa_to_formula;
		  map<string, string>::iterator iter1;

		  map<string, int> f;
		  map<string, int>::iterator iter2;

		  aa_to_formula["G"] = "C2H3N1O1";
		  aa_to_formula["A"] = "C3H5N1O1";
		  aa_to_formula["S"] = "C3H5N1O2";
		  aa_to_formula["P"] = "C5H7N1O1";
		  aa_to_formula["V"] = "C5H9N1O1";
		  aa_to_formula["T"] = "C4H7N1O2";
		  aa_to_formula["C"] = "C3H5N1O1S1";
		  aa_to_formula["L"] = "C6H11N1O1";
		  aa_to_formula["I"] = "C6H11N1O1";
		  aa_to_formula["N"] = "C4H6N2O2";
		  aa_to_formula["D"] = "C4H5N1O3";
		  aa_to_formula["Q"] = "C5H8N2O2";
		  aa_to_formula["K"] = "C6H12N2O1";
		  aa_to_formula["E"] = "C5H7N1O3";
		  aa_to_formula["M"] = "C5H9N1O1S1";
		  aa_to_formula["H"] = "C6H7N3O1";
		  aa_to_formula["F"] = "C9H9N1O1";
		  aa_to_formula["R"] = "C6H12N4O1";
		  aa_to_formula["Y"] = "C9H9N1O2";
		  aa_to_formula["W"] = "C11H10N2O1";

		  //translate each aa into formula, and
		  //sum up the elements
		  for (int i=0; i< (int) seq.length(); i++) {
			string aa = string(seq, i, 1);
			iter1 = aa_to_formula.find(aa);
			if ( iter1 != aa_to_formula.end() ){
			  string aa_formula = iter1->second;
			  map<string, int> temp_f = readMoleculeFormula(aa_formula);
			  for (iter2 = temp_f.begin(); 
			   iter2 != temp_f.end(); iter2++) {
			f[iter2->first] += iter2->second;
			  }
			}
		  }

		  //add any modification here???
		  for (int i=0; i< (int) seq.length(); i++) {
			string aa = string(seq, i, 1);
			//add carboxymide
			if ( aa.compare("C") == 0 ) {
			  f["H"] += 1;
			  f["O"] += 2;
			  f["C"] += 2;
			}
		  }

		  //add an H2O and and a few protons
		  f["H"] += 2; 
		  f["O"] += 1;

		  string result = elementListToFormula(f);  

		  //convert the elements to a string
		  /*
		  string result = "";
		  char line[256];
		  for (iter2 = f.begin(); 
			   iter2 != f.end(); iter2++) {
			sprintf(line, "%s%d", iter2->first.c_str(), iter2->second);
			result.append(string(line));
		  }
		  */

		  return result;
		}


		//this function has not been verified for correctness
		//yet. It is not right, most probably ???
		//so far, it has never been used any more.
		double Elemental::formulaToMass(string& formula) 
		{
		  map<string, int> f;
		  map<string, int>::iterator iter1;

		  //double unit_H = 1.0078;
		  f = readMoleculeFormula(formula);

		  double mass = 0;
		  for (iter1 = f.begin(); 
			   iter1 != f.end(); iter1++) {
			if ( iter1->first == "C" )
			  mass += iter1->second * 12;
			else if ( iter1->first == "H" )
			  mass += iter1->second * 1;
			else if ( iter1->first == "O" )
			  mass += iter1->second * 16;
			else if ( iter1->first == "N" )
			  mass += iter1->second * 14;
			else if ( iter1->first == "S" )
			  mass += iter1->second * 32;
		  }

			return mass * PROTON;
		}

		//returns the element name and its corresponding quantity
		map<string, int> Elemental::readMoleculeFormula(string& formula)
		{
		  string element_name = "";
		  int element_quantity = 0;
		  map<string, int> result;

		  for (int i=0; i< (int) formula.length(); i++) 
		  {
			char c = formula[i];

			if ( c >= 'A' && c <= 'z' ) 
			{
			  //ready to read next element
			  element_name = string(formula, i, 1);
			  element_quantity = 0;

			  for (int j=i+1; j < (int) formula.length(); j++) 
			  {
				char c2 = formula[j];
				if (! (c2 >= '0' && c2 <= '9') )
				  break;
				else 
				{
				  char cp[2];
				  cp[0] = c2; cp[1] = '\0';
				  element_quantity = element_quantity * 10 + atoi(cp);
				}
			  }
			  
			  if ( element_quantity != 0 )
			  {
				Elemental ele = Elemental(element_name);
				if (! ele.isEmpty() ) 
				  result[ele.getName()] = element_quantity;
				else
				  cout << "wrong element = " << element_name << endl;
			  }
			}
		  }

		  return result;
		}


		//trim the envelope by removing the very small terms.
		void Elemental::trimEnvelope(vector<double>& env) 
		{
			std::vector <double> env1 ; 
			  int trim_index = -1;

			  for (int i= (int) env.size()-1; i>=0; i--) {
				if ( env[i] < EPSILON )
				  trim_index = i;
				else
				  break;
			  }

			  if ( trim_index >= 0) {
				vector<double>::iterator iter1 = env.begin() + trim_index;
				env1.assign(env.begin(), iter1);
				env.clear()  ;
				env.assign(env1.begin(), env1.end()) ; 
			  }
		}


			//returns the series of relative intensity starting from
			//the first mono isotope. the intensity is normalized 
			//regarding to the total intensity
		vector<double> Elemental::getIsotopicEnvelope(string &formula) 
		{
			  map<string, int> the_formula = readMoleculeFormula(formula);
			  map<string, int>::iterator iter1;
			  vector<double> result;
			  vector<double> single_ele_envelope;

			  //cout << "point 5" << formula << endl;
			  for (iter1 = the_formula.begin(); 
				   iter1 != the_formula.end(); iter1 ++) {
				string ele_name = iter1->first; 
				Elemental ele = Elemental(ele_name);
				int quantity = iter1->second;

				//cout << "point 6:" << ele_name << endl;
				single_ele_envelope = getSingleElementEnvelope(ele, quantity);
					  //cout << "point 8" << endl;
				result = fuseIsotopicEnvelopes(result, single_ele_envelope);
					  //cout << "point 7" << endl;
			  }

			  return result;
		}	


		//return the kth highest isotopic peak of the
		//given formula. k must be between 1 and 5.
		//e.g. when k is 2, the index of the 2nd highest isotopic peak
		//is returned. NOTE: all indices starts from 0
		// In case of error , -1 is returned.
		int Elemental::highestIsotopeOfFormula(int k, string &formula)
		{
		  //cout << "formula here : " << formula << endl;
		  Elemental ele = Elemental(string("C"));

		  if ( k <=0 || k > 5 )
			return -1;
		  else {
			vector<double> envelope = ele.getIsotopicEnvelope(formula);
			//cout << "point 3" << endl;
			if ( ((int) envelope.size()) < k )
			  return -1;
			else {
			 // cout << "point 1" << endl;
			  vector<double> tmp;
			  tmp.assign(envelope.begin(), envelope.end());
			  sort(tmp.begin(), tmp.end(), greater<double>());

			  double kth_highest_intensity = tmp[k-1];
		      
			  //cout << "point 2" << endl;

			  //look for the index of the kth highest
			  int result = -1;
			  for (int i=0; i< (int) envelope.size(); i++) {
			  if ( envelope[i] == kth_highest_intensity ) {
				result = i;
				break;
			  }
			  }

			  return result;
			}
		  }
		}


		//compute the isotopic envelope for a molecule which contains
		//only one element. this result will be used for fusing multiple
		//single element envelopes
		vector<double> Elemental::getSingleElementEnvelope3(Elemental& ele, int quantity) 
		{
		  vector<double> result;
		  map<int, double> temp_isotope = ele.getIsotopes();
		  map<int, double>::iterator iter1;
		  int min_isotope_mass = 0;
		  int max_isotope_mass = 0;

		  //find the minimum and maximum isotopic mass
		  for (iter1 = temp_isotope.begin(); 
			   iter1 != temp_isotope.end(); iter1 ++ ) {
			if ( min_isotope_mass == 0 ||
			 min_isotope_mass > iter1->first )
			  min_isotope_mass = iter1->first;

			if ( max_isotope_mass < iter1->first )
			  max_isotope_mass = iter1->first;    
		  }

		  //calculate the evenlope. 
		  //right now, deal with only three isotope factions p, q and r
		  //find the minimum isotopic mass
		  int mass[3];
		  double fract[3];
		  int c = 0;
		  for (iter1 = temp_isotope.begin(); 
			   iter1 != temp_isotope.end(); iter1 ++ ) {
			if ( c < 3 ) {
			  mass[c] = iter1->first;
			  fract[c] = iter1->second;
			}
			c++;
		  }

		  result = vector<double>(max_isotope_mass * quantity, 0.0);
		  for (int i=0; i<= quantity; i++) {
			  double fraction1 = Engine::Utilities::NChooseK(quantity, i) *
			  pow((double) fract[1], (double) i);
		    
			for (int j=0; j<= quantity -i; j++ ) {
			  double fraction2;

			  if ( quantity - i == 0 )
			fraction2 = 1;
			  else fraction2 = Engine::Utilities::NChooseK(quantity-i, j)*
			pow((double) fract[2], (double) j) *
			pow((double) fract[0], (double) quantity-i-j);
		      
			  double fraction = fraction1 * fraction2;
			  int m = i * mass[1] + j * mass[2] + (quantity-i-j)*mass[0];
			  int isotope_index = m - quantity * min_isotope_mass;
		      
			  if ( fraction > 0 )
			result[isotope_index] += fraction;
			}
		  }

		  trimEnvelope(result);
		  return result;
		}


		//compute the isotopic envelope for a molecule which contains
		//only one element. this result will be used for fusing multiple
		//single element envelopes
		vector<double> Elemental::getSingleElementEnvelope(Elemental& ele, int quantity)
		{
		  vector<double> result;
		  map<int, double> temp_isotope = ele.getIsotopes();
		  map<int, double>::iterator iter1;
		  int min_isotope_mass = 0;
		  int max_isotope_mass = 0;

		  if ( temp_isotope.size() <= 2 ) {
			//find the minimum and maximum isotopic mass
			for (iter1 = temp_isotope.begin(); 
			 iter1 != temp_isotope.end(); iter1 ++ ) {
			  if ( min_isotope_mass == 0 ||
			   min_isotope_mass > iter1->first )
			min_isotope_mass = iter1->first;

			  if ( max_isotope_mass < iter1->first )
			max_isotope_mass = iter1->first;    
			}

			//calculate the evenlope. 
			//right now, deal with only two isotope factions p and q !!!
			//find the minimum isotopic mass
			int mass[2];
			double fract[2];
			int c = 0;
			for (iter1 = temp_isotope.begin(); 
			 iter1 != temp_isotope.end(); iter1 ++ ) {
			  if ( c < 2 ) {
			mass[c] = iter1->first;
			fract[c] = iter1->second;
			  }
			  c++;
			}

			result = vector<double>(max_isotope_mass * quantity, 0.0);

			for (int i=0; i<= quantity; i++) {
				double fraction = Engine::Utilities::NChooseK(quantity, i)*
			pow((double) fract[0], (double) i) *
			pow((double) fract[1], (double) quantity-i);
			  int m = i * mass[0] + (quantity-i)*mass[1];
			  int isotope_index = m - quantity * min_isotope_mass;

			  if ( fraction > 0 )
			result[isotope_index] = fraction;
			}
		  }
		  else
			result = getSingleElementEnvelope3(ele, quantity);

		  trimEnvelope(result);
		  return result;
		}

		//combine multiple isotopic envelopes 
		vector<double> Elemental::fuseIsotopicEnvelopes(const vector<double>& env1,  const vector<double>& env2)
		{
		  vector<double> result;

		  if ( env1.size() <= 0 )
			result.assign(env2.begin(), env2.end());
		  else if ( env2.size() <= 0 )
			result.assign(env1.begin(), env1.end());
		  else {
			result = vector<double>((int) (env1.size() + env2.size()-1), 0);

			for (int i=0; i< (int) env1.size(); i++) {
			  for (int j=0; j< (int) env2.size(); j++) {
			int isotope_index = i+j;
			double delta_fraction = env1[i] * env2[j];
			result[isotope_index] += delta_fraction;
			  }
			}
			//result.assign(env2.begin(), env2.end());
		  }

		  trimEnvelope(result); 
		  return result;
		}  

	}
}