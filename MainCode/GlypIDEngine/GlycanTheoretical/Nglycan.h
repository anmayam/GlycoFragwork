// Written by Yin Wu (wuyin@indiana.edu)
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0


// modified by Anoop for GlypID
// Anoop Comment: Should thinking about breaking this out into diff files to avoid a cluster-blah


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <fstream>
#include <string.h>
#include <assert.h>
using namespace std;

namespace Engine
{
	namespace GlycanTheoretical
	{
		enum { MONOSACCHARIDE_HEXNAC = 0, MONOSACCHARIDE_HEX,
			  MONOSACCHARIDE_DEHEX, MONOSACCHARIDE_NEUAC, 
			  MONOSACCHARIDE_NEUGC, MONOSACCHARIDE_MOD_SULFATE, MONOSACCHARIDE_TOTAL, MONOSACCHARIDE_UNKNOWN };

		
		
		enum {COMPLEX_ANTENNA = 50, HIGH_MANNOSE_ANTENNA };

		enum   {COMPLEX_GLYCAN = 100, HIGH_MANNOSE_GLYCAN, HYBRID_GLYCAN, ALL_NGLYCAN};

		enum GLYCAN_MONOSACCHARIDE {GLYCAN_HAS_THIS_MONOSACCHARIDE_POSSIBLE = 150, 
		GLYCAN_HAS_THIS_MONOSACCHARIDE_YES, GLYCAN_HAS_THIS_MONOSACCHARIDE_NO, 
		GLYCAN_HAS_THIS_MONOSACCHARIDE_UNKNOWN};
	
		
		
	/*	map<int, float> monosaccharide_mass_accurate;
		map<int, float> monosaccharide_mass_average;
		map<int, string> monosaccharide_names;*/

		//class NGlycan;

		class GlycanComposition 
		{
			map<int, float> monosaccharide_mass_accurate;
			map<int, float> monosaccharide_mass_average;
			map<int, string> monosaccharide_names;
			int glycan_has_hexose, glycan_has_fucose, glycan_has_neuac, glycan_has_hexnac, glycan_has_neugc;
			int glycan_has_sulfate ; 
		


			public:

			  void initNGlycanConstants();
			  //key is monosaccharide, value is the count of it.
			  map<int, int> monosaccharide_count;  
			  int glycan_type;
			  bool consider_sulfated ; 
			  bool is_decoy ; 

			  float glycan_avg_mass; 
			  float glycan_accurate_mass ; 
			  

			  GlycanComposition();
			  void clone(GlycanComposition &new_gc);
			  void print();
			  GlycanComposition combine(GlycanComposition &gc);
			  float accurateMass();
			  bool hasNeuAc() ; 
			  bool hasDeHex();
			  /*int numHex() ; 
			  int numHexNAc() ; 
			  int numDeHex() ; 
			  int numNeuAc() ; */
			  
			  float averageMass();
			  int getID();
			  int totalNumMonosacchrides();
			  bool validate(GlycanComposition &min_gc,  GlycanComposition &max_gc,	int min_size, int max_size, int glycan_type);
			  string getCompositionString(); //returns the string form
			  string getConvertableString();
			  void parseCompositionString(string gc_comp) ; 
		};

		//class Bond 
		//{
		//	public:
		//	  int parent_side_bond, child_side_bond; 
		//};

		//class Monosaccharide
		//{
		//	public:
		//	  int mono_type;
		//	  Bond incoming_bond;
		//	  Monosaccharide *parent;
		//	  vector<Monosaccharide*> children;

		//	  Monosaccharide(int type, Monosaccharide *parent);
		//	  Monosaccharide();
		//	  void addChild(Monosaccharide *child);

		//	private:
		//	  static map<int, string> monosaccharide_names;
		//	  static map<int, float> monosaccharide_avg_mass;
		//	  static map<int, float> monosaccharide_mono_mass;

		//	  static void init();
		//};

		//class NGTree 
		//{
		//	private:
		//		  void deleteNGTree(Monosaccharide *root);
		//		  void getComposition(Monosaccharide* root, GlycanComposition &gc);
		//	public:
		//		  Monosaccharide *root;
		//		  bool with_fucose;

		//		  NGTree();
		//		  ~NGTree();
		//		  NGTree* clone();
		//		  Monosaccharide* cloneNGTree(Monosaccharide* root);
		//		  void copyInto(NGTree *destination);
		//		  void setRoot(Monosaccharide* root);
		//		  GlycanComposition getComposition();
		//		  int getNGTreeSize(Monosaccharide *root);
		//		  int size();
		//		  static void allPossibleNGlycans(std::vector<NGlycan*> &sub_trees);			
		//};


		//class NGlycan
		//{
		//	public:
		//	  NGTree pentamer;
		//	  vector<NGTree> antennas;

		//	  void print();
		//	  GlycanComposition getComposition();
		//	  void getIsoformCompositions(vector<GlycanComposition> &isoform_comp, GlycanComposition min_composition, GlycanComposition max_composition);
		//	 // bool validateGlycan(GlycanComposition &min_gc, GlycanComposition &max_gc,int max_antenna_fucose, int max_pentamer_fucose,int min_size, int max_size);
		//	  int getGlycanType();
		//	  int size();
		//};


		//class PentamerTemplate 
		//{
		//	public:
		//	  NGTree ng_template;

		//	  PentamerTemplate(bool has_fucose);
		//	  void getAllSubStructures(vector<NGTree*> &sub_trees);
		//	  void combineSubtreesOfChildren(vector<NGTree*> &result,  vector< vector<NGTree*> > &subtrees_of_all_children, int index);
		//	  void getAllSub(vector<NGTree*> &sub_trees, Monosaccharide *root);
		//};
	
		//class AntennaTemplate 
		//{
		//	private:
		//	  Monosaccharide* makeComplexAntenna(bool has_fucose);
		//	  Monosaccharide* makeHMAntenna(bool has_fucose);

		//	  void getAllSub(vector<NGTree*> &sub_trees, Monosaccharide *root);

		//	public:
		//	  NGTree ng_template;
		//	  vector<Monosaccharide> nodes;
		//	  bool with_fucose;

		//	  AntennaTemplate(int type, bool has_fucose);
		//	  void getAllSubStructures(vector<NGTree*> &sub_trees);
		//	
		//};


		//extern void printNGTree(Monosaccharide *root);
		//extern int strToMonosaccharideCardinality(string str);

		bool checkCardinality(int n, int cardinality);
	}
}

