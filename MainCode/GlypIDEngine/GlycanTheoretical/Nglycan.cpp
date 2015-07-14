//// Written by Yin Wu (wuyin@indiana.edu)
//// -------------------------------------------------------------------------------
//// 
//// Licensed under the Apache License, Version 2.0; you may not use this file except
//// in compliance with the License.  You may obtain a copy of the License at 
//// http://www.apache.org/licenses/LICENSE-2.0
//
//// modified by Anoop for GlypID
//
#include "Nglycan.h"
#include "../Utilities/system.h" 
#include "../Utilities/GStrTok.h"
//
namespace Engine
{
	namespace GlycanTheoretical
	{
//

//
		void GlycanComposition::initNGlycanConstants() 
		{
			monosaccharide_mass_accurate[MONOSACCHARIDE_HEXNAC] = 203.0794;
			monosaccharide_mass_accurate[MONOSACCHARIDE_HEX] = 162.0528;
			monosaccharide_mass_accurate[MONOSACCHARIDE_NEUAC] = 291.0954;
			monosaccharide_mass_accurate[MONOSACCHARIDE_DEHEX] = 146.0579;
			monosaccharide_mass_accurate[MONOSACCHARIDE_NEUGC] = 307.0903;
			monosaccharide_mass_accurate[MONOSACCHARIDE_MOD_SULFATE]= 80.0632 ; 

			monosaccharide_mass_average[MONOSACCHARIDE_HEXNAC] = 203.195;
			monosaccharide_mass_average[MONOSACCHARIDE_HEX] = 162.1424;
			monosaccharide_mass_average[MONOSACCHARIDE_NEUAC] = 291.2579;
			monosaccharide_mass_average[MONOSACCHARIDE_DEHEX] = 146.1463;
			monosaccharide_mass_average[MONOSACCHARIDE_NEUGC] = 307.2573;
			monosaccharide_mass_average[MONOSACCHARIDE_MOD_SULFATE] = 80.0632 ; 

			monosaccharide_names[MONOSACCHARIDE_HEXNAC] = "HexNAc";
			monosaccharide_names[MONOSACCHARIDE_HEX] = "Hex";
			monosaccharide_names[MONOSACCHARIDE_NEUAC] = "NeuAc";
			monosaccharide_names[MONOSACCHARIDE_NEUGC] = "NeuGc";
			monosaccharide_names[MONOSACCHARIDE_DEHEX] = "DeHex";
			monosaccharide_names[MONOSACCHARIDE_MOD_SULFATE] = "SO3" ; 

			consider_sulfated = false ; 
			glycan_accurate_mass = 0 ; 
			glycan_avg_mass = 0 ; 
		}

		//void printNGTree(Monosaccharide *root) 
		//{
		//  if ( root == NULL )
		//	return;
		//  else {
		//	if ( root->parent == NULL || 
		//	 root->children.size() > 0 ) //it is the root of the whole tree
		//	  cout << "(";

		//	cout << root->mono_type << ",";

		//	if ( root->children.size() > 0 ) {
		//	  cout << "(";
		//	  for (int i=0; i<(int) root->children.size(); i++)
		//	printNGTree(root->children[i]);
		//	  cout << ")";
		//	}

		//	if ( root->parent == NULL || 
		//	 root->children.size() > 0 ) 
		//	  cout << ")";
		//  }
		//}

		//void Monosaccharide::init() 
		//{
		//  /*
		//  monosaccharide_names.clear();

		//  monosaccharide_names.insert(make_pair((int) MONOSACCHARIDE_HEXNAC, "HexNAc"));
		//  monosaccharide_names.insert(make_pair((int) MONOSACCHARIDE_HEX, "Hex"));
		//  monosaccharide_names.insert(make_pair((int) MONOSACCHARIDE_DEHEX, "DeHex"));
		//  monosaccharide_names.insert(make_pair((int) MONOSACCHARIDE_NEUAC, "NeuAc"));
		//  */
		//}

		//Monosaccharide::Monosaccharide(int type, Monosaccharide *parent)
		//{
		//  /*
		//  if ( monosaccharide_names.size() <= 0 )
		//	init();
		//  */

		//  this->mono_type = type;
		//  this->parent = parent;
		//}


		//Monosaccharide::Monosaccharide() {
		//  /*
		//  if ( monosaccharide_names.size() <= 0 )
		//	init();
		//  */

		//  this->mono_type = MONOSACCHARIDE_UNKNOWN;
		//  this->parent = NULL;
		//}


		//void Monosaccharide::addChild(Monosaccharide *child) 
		//{

		//}


		//PentamerTemplate::PentamerTemplate(bool has_fucose) 
		//{

		//  //build a full pentamer structure
		//  Monosaccharide *temp_root, *temp_mono1, *temp_mono2;

		//  //add first HexNac
		//  temp_mono1 = new Monosaccharide(MONOSACCHARIDE_HEXNAC, NULL);
		//  temp_root = temp_mono1;

		//  if ( has_fucose ) {
		//	//add first Fucose
		//	/*
		//	temp_mono2 = new Monosaccharide(MONOSACCHARIDE_DEHEX, temp_mono1);
		//	temp_mono1->children.push_back(temp_mono2);
		//	*/
		//  }

		//  this->ng_template.with_fucose = has_fucose;

		//  //add second HexNac
		//  temp_mono2 = new Monosaccharide(MONOSACCHARIDE_HEXNAC, temp_mono1);
		//  temp_mono1->children.push_back(temp_mono2);
		//  temp_mono1 = temp_mono2;

		//  //add First Mannose
		//  temp_mono2 = new Monosaccharide(MONOSACCHARIDE_HEX, temp_mono1);
		//  temp_mono1->children.push_back(temp_mono2);
		//  temp_mono1 = temp_mono2;
		//  
		//  //add two child mannose
		//  temp_mono2 = new Monosaccharide(MONOSACCHARIDE_HEX, temp_mono1);
		//  temp_mono1->children.push_back(temp_mono2);
		//  temp_mono2 = new Monosaccharide(MONOSACCHARIDE_HEX, temp_mono1);
		//  temp_mono1->children.push_back(temp_mono2);

		//  this->ng_template.root = temp_root;
		//  //printNGTree(temp_root);
		//}



		//AntennaTemplate::AntennaTemplate(int type, bool has_fucose)aa 
		//{
		//  switch ( type ) {
		//  case COMPLEX_ANTENNA:
		//	this->ng_template.root = makeComplexAntenna(has_fucose);
		//	break;
		//  case HIGH_MANNOSE_ANTENNA:
		//	this->ng_template.root = makeHMAntenna(has_fucose);
		//	break;
		//  default:
		//	assert ( false );
		//  }
		//}


		//NGTree::NGTree()
		//{
		//  this->root = NULL;
		//  this->with_fucose = false;
		//}


		//void NGTree::setRoot(Monosaccharide* root)
		//{
		//}

		//void NGTree::deleteNGTree(Monosaccharide* root) 
		//{

		//  if ( root != NULL ) {
		//	for ( int i=0; i< (int) root->children.size(); i++ ) {
		//	  Monosaccharide *mono = root->children[i];
		//	  delete( mono );
		//	}
		//	delete( root );
		//  }
		//}

		//NGTree::~NGTree() {

		//  if ( root != NULL ) {
		//	deleteNGTree( root );
		//  }
		//}


		////==================DEBUG ZONE===========================
		//void AntennaTemplate::getAllSubStructures(vector<NGTree*> &sub_trees)
		//{
		//  //getAllSub(sub_trees, this->ng_template.root);
		//  
		//  NGTree *new_sub_tree = new NGTree();

		//  sub_trees.push_back(new_sub_tree); //empty tree is a default subtree

		//  for ( int i=0; i < (int) this->nodes.size(); i++) {
		//	new_sub_tree = new NGTree();
		//	new_sub_tree->with_fucose = this->with_fucose;
		//	new_sub_tree->root = new Monosaccharide();
		//	(* (new_sub_tree->root)) = this->nodes[0];
		//	new_sub_tree->root->parent = NULL;
		//	Monosaccharide *tail = new_sub_tree->root;
		//	Monosaccharide *temp;

		//	for (int j=1; j <= i; j++) {
		//	  temp = new Monosaccharide();
		//	  (*temp) = this->nodes[j];
		//	  temp->parent = tail;
		//	  tail->children.push_back( temp );
		//	  tail = temp;
		//	}

		//	sub_trees.push_back(new_sub_tree);
		//  }


		//  //add this part for NEUGC.
		//  //consider changing it soon!!!????
		//  //replace the NeuAc with NeuGc and rebuild the tree.
		//  if ( this->nodes.back().mono_type == MONOSACCHARIDE_NEUAC ) 
		//  {
		//	Monosaccharide tmp_mono = this->nodes.back();
		//	Monosaccharide neugc = 
		//	  Monosaccharide(MONOSACCHARIDE_NEUGC, NULL);
		//      
		//	this->nodes[this->nodes.size()-1] = neugc;

		//	new_sub_tree = new NGTree();
		//	new_sub_tree->with_fucose = this->with_fucose;
		//	new_sub_tree->root = new Monosaccharide();
		//	(* (new_sub_tree->root)) = this->nodes[0];
		//	new_sub_tree->root->parent = NULL;
		//	Monosaccharide *tail = new_sub_tree->root;
		//	Monosaccharide *temp;

		//	for (int j=1; 
		//	 j < (int) this->nodes.size(); j++) {
		//	  temp = new Monosaccharide();
		//	  (*temp) = this->nodes[j];
		//	  temp->parent = tail;
		//	  tail->children.push_back( temp );
		//	  tail = temp;
		//	}

		//	sub_trees.push_back(new_sub_tree);
		//	this->nodes[this->nodes.size()-1] = tmp_mono;
		//  }
		//  
		//}


		//void AntennaTemplate::getAllSub(vector<NGTree*> &sub_trees,
		//				Monosaccharide *root) 
		//{
		//  /*
		//  NGTree *new_sub_tree;

		//  if ( this == NULL || this->nodes.size() <= 0 )
		//	return;
		//  else {
		//	new_sub_tree = new NGTree();
		//	new_sub_tree->root = new Monosaccharide(this->mono_type, 
		//					   this->parent);
		//	(*(new_sub_tree->root)) = (*this); //copy the content over.
		//	Monosaccharide *temp_mono = this->parent;

		//	//reversly add the chain
		//	while ( parent_mono != NULL ) {
		//	  new_sub_tree->root->parent = temp_mono;
		//	  temp_mono = new_sub_tree->root;
		//	}
		//  }
		//  */
		//}

		//void PentamerTemplate::getAllSub(vector<NGTree*> &sub_trees,
		//				 Monosaccharide *root) 
		//{
		//  vector< vector<NGTree*> > subtrees_of_all_children;
		//  vector<NGTree*> result;

		//  //get all the subtrees for all of the children
		//  for (int i=0; i< (int) root->children.size(); i++) {
		//	if ( root->children[i] != NULL ) {
		//	  vector<NGTree*> subtrees_of_child;
		//	  subtrees_of_all_children.push_back(subtrees_of_child);
		//	  getAllSub(subtrees_of_all_children[i], root->children[i]);    
		//	}
		//  }

		//  //use a recursive function to combine all possible subtrees 
		//  //of all the children. save them in the result
		//  /*
		//	//recursive implementation
		//  if ( subtrees_of_all_children.size() > 0 ) {
		//	combineSubtreesOfChildren(result, subtrees_of_all_children, 0);
		//    
		//	for (int i=0; i< (int) result.size(); i++)
		//	  result[i]->setRoot(root);
		//  }
		//  */

		//  /*
		//	//for loop implementation
		//  if ( subtrees_of_all_children.size() > 0 ) {
		//	NGTree sub_t;

		//	//first child
		//	for (int i=0; i < (int) subtrees_of_all_children[0].size(); i++ ) {
		//	  NGTree first_child = subtrees_of_all_children[0][i]->clone();

		//	  if ( subtrees_of_all_children.size() > 1 ) {
		//	//second child
		//	for (int j=0; j < (int) subtrees_of_all_children[1].size(); j++ ) {
		//	  NGTree *second_child = subtrees_of_all_children[1][i];
		//	}
		//	  }
		//	} 
		//  }
		//  */


		//  //clear memory
		//  for (int i=0; i< (int) subtrees_of_all_children.size(); i++) 
		//  {
		//	vector<NGTree*> subtrees_of_child = subtrees_of_all_children[i];

		//	for (int j=0; j< (int) subtrees_of_child.size(); j++) {
		//	  delete( subtrees_of_child[j] );
		//	}

		//	subtrees_of_all_children[i].clear();
		//  }
		//  subtrees_of_all_children.clear();

		//  sub_trees.assign(result.begin(), result.end());
		//  result.clear();
		//  return;
		//}



		//*
		//NGTree::NGTree() {
		//  root = NULL;
		//}
		//*/

		//NGTree* NGTree::clone()
		//{
		//  NGTree *result = NULL;
		//  return result;
		//}

		//void PentamerTemplate::combineSubtreesOfChildren(vector<NGTree*> &result, 
		//			  vector< vector<NGTree*> > &subtrees_of_all_children, 
		//			  int index) 
		//{


		//}


		//void PentamerTemplate::getAllSubStructures(vector<NGTree*> &sub_trees) 
		//{
		//  getAllSub(sub_trees, this->ng_template.root);
		//}



		//
		//Monosaccharide* AntennaTemplate::makeComplexAntenna(bool has_fucose) {
		//  //build a full pentamer structure
		//  Monosaccharide *temp_root, *temp_mono1, *temp_mono2;

		//  //add first HexNac
		//  temp_mono1 = new Monosaccharide(MONOSACCHARIDE_HEXNAC, NULL);
		//  temp_root = temp_mono1;

		//  if ( has_fucose ) {
		//	//add first Fucose
		//	temp_mono2 = new Monosaccharide(MONOSACCHARIDE_DEHEX, temp_mono1);
		//	temp_mono1->children.push_back(temp_mono2);
		//  }

		//  //add First Gal
		//  temp_mono2 = new Monosaccharide(MONOSACCHARIDE_HEX, temp_mono1);
		//  temp_mono1->children.push_back(temp_mono2);
		//  temp_mono1 = temp_mono2;
		//  
		//  //add First NeuAc
		//  temp_mono2 = new Monosaccharide(MONOSACCHARIDE_NEUAC, temp_mono1);
		//  temp_mono1->children.push_back(temp_mono2);

		//  printNGTree(temp_root);
		//  return temp_root;
		//}


		//Monosaccharide* AntennaTemplate::makeHMAntenna(bool has_fucose) {
		//  //build a full pentamer structure
		//  Monosaccharide *temp_root, *temp_mono1, *temp_mono2;

		//  //add First Mannose
		//  temp_mono1 = new Monosaccharide(MONOSACCHARIDE_HEX, NULL);
		//  temp_root = temp_mono1;

		//  if ( has_fucose ) {
		//	//add first Fucose
		//	temp_mono2 = new Monosaccharide(MONOSACCHARIDE_DEHEX, temp_mono1);
		//	temp_mono1->children.push_back(temp_mono2);
		//  }

		//  //add Second Mannose
		//  temp_mono2 = new Monosaccharide(MONOSACCHARIDE_HEX, temp_mono1);
		//  temp_mono1->children.push_back(temp_mono2);
		//  temp_mono1 = temp_mono2;
		//  
		//  //add Third Mannose
		//  temp_mono2 = new Monosaccharide(MONOSACCHARIDE_HEX, temp_mono1);
		//  temp_mono1->children.push_back(temp_mono2);

		//  printNGTree(temp_root);
		//  return temp_root;
		//}
		//*/


		//Monosaccharide* AntennaTemplate::makeComplexAntenna(bool has_fucose)
		//{

		//  this->nodes.clear();

		//  this->with_fucose = has_fucose;

		//  this->nodes.push_back( Monosaccharide(MONOSACCHARIDE_HEXNAC, NULL) );
		//  this->nodes.push_back( Monosaccharide(MONOSACCHARIDE_HEX, NULL) );
		//  this->nodes.push_back( Monosaccharide(MONOSACCHARIDE_NEUAC, NULL) );

		//  return NULL;
		//}


		//Monosaccharide* AntennaTemplate::makeHMAntenna(bool has_fucose) 
		//{
		//  this->nodes.clear();

		//  this->with_fucose = has_fucose;

		//  this->nodes.push_back( Monosaccharide(MONOSACCHARIDE_HEX, NULL) );
		//  this->nodes.push_back( Monosaccharide(MONOSACCHARIDE_HEX, NULL) );
		//  this->nodes.push_back( Monosaccharide(MONOSACCHARIDE_HEX, NULL) );
		//  this->nodes.push_back( Monosaccharide(MONOSACCHARIDE_HEX, NULL) );

		//  return NULL;
		//}


		//void NGTree::allPossibleNGlycans(vector<NGlycan*> &sub_trees) 
		//{
		//  PentamerTemplate pentamer_with_fucose = PentamerTemplate(true);
		//  PentamerTemplate pentamer_no_fucose = PentamerTemplate(false);

		//  AntennaTemplate complex_antenna_with_fucose = AntennaTemplate(COMPLEX_ANTENNA,
		//							true);
		//  AntennaTemplate complex_antenna_no_fucose = AntennaTemplate(COMPLEX_ANTENNA,
		//							  false);  
		//  AntennaTemplate high_mannose_antenna_with_fucose = 
		//	AntennaTemplate(HIGH_MANNOSE_ANTENNA, true);
		//  AntennaTemplate high_mannose_antenna_no_fucose = 
		//	AntennaTemplate(HIGH_MANNOSE_ANTENNA, false);

		//  vector<NGTree*> all_fucosed_complex_antennas, all_no_fucose_complex_antennas, 
		//	all_fucosed_high_mannose_antennas, all_no_fucose_high_mannose_antennas,
		//	all_antennas;
		//  vector<PentamerTemplate*> all_pentamers;

		//  NGlycan *ng;

		//  //get all possible pentamers
		//  all_pentamers.push_back( &pentamer_with_fucose );
		//  all_pentamers.push_back( &pentamer_no_fucose );

		//  //add high-mannose type trees
		//  //1. get all the structure and sub-structures of antennas
		//  complex_antenna_with_fucose.getAllSubStructures(all_fucosed_complex_antennas);
		//  complex_antenna_no_fucose.getAllSubStructures(all_no_fucose_complex_antennas);
		//  high_mannose_antenna_with_fucose.getAllSubStructures(all_fucosed_high_mannose_antennas);
		//  high_mannose_antenna_no_fucose.getAllSubStructures(all_no_fucose_high_mannose_antennas);


		//  //2. get all the substrcutre of the antennas
		//  if ( glycopeptide_matching_glycan_type == HYBRID_GLYCAN ||
		//	glycopeptide_matching_glycan_type == ALL_NGLYCAN ||
		//	   glycopeptide_matching_glycan_type == COMPLEX_GLYCAN ) {
		//	all_antennas.insert(all_antennas.end(), all_fucosed_complex_antennas.begin(),
		//			all_fucosed_complex_antennas.end());
		//	all_antennas.insert(all_antennas.end(), all_no_fucose_complex_antennas.begin(),
		//			all_no_fucose_complex_antennas.end());
		//  }
		//  
		//  if ( glycopeptide_matching_glycan_type == HYBRID_GLYCAN ||
		//	glycopeptide_matching_glycan_type == ALL_NGLYCAN ||
		//	   glycopeptide_matching_glycan_type == HIGH_MANNOSE_GLYCAN ) {
		//	all_antennas.insert(all_antennas.end(), all_fucosed_high_mannose_antennas.begin(),
		//			all_fucosed_high_mannose_antennas.end());
		//	all_antennas.insert(all_antennas.end(), all_no_fucose_high_mannose_antennas.begin(),
		//			all_no_fucose_high_mannose_antennas.end());
		//  }


		//  if ( glycopeptide_matching_glycan_type != HYBRID_GLYCAN &&
		//	   glycopeptide_matching_glycan_type != COMPLEX_GLYCAN &&
		//	   glycopeptide_matching_glycan_type != HIGH_MANNOSE_GLYCAN &&
		//	   glycopeptide_matching_glycan_type != ALL_NGLYCAN ) {
		//	cout << "Unknown glycan type: " << glycopeptide_matching_glycan_type << endl;
		//	assert(0);
		//  }
		//  
		//  //3.now combine all of them to produce all possible complex trees
		//  //a complete glycan tree can have one pentamer and up to 4 antennas
		//  //the structure contains both complex, high-mannose and hybrid
		//  for (int i=0; i< (int) all_pentamers.size(); i++) {
		//	for (int j=0; j< (int) all_antennas.size(); j++) {
		//	  for (int k=0; k< (int) all_antennas.size(); k++) {
		//	for (int m=0; m< (int) all_antennas.size(); m++) {
		//	  for (int n=0; n< (int) all_antennas.size(); n++) {
		//		ng = new NGlycan();
		//		ng->antennas.push_back(NGTree());
		//		ng->antennas.push_back(NGTree());
		//		ng->antennas.push_back(NGTree());
		//		ng->antennas.push_back(NGTree());
		//		all_pentamers[i]->ng_template.copyInto(&(ng->pentamer));
		//		all_antennas[j]->copyInto(&(ng->antennas[0]));
		//		all_antennas[k]->copyInto(&(ng->antennas[1]));
		//		all_antennas[m]->copyInto(&(ng->antennas[2]));
		//		all_antennas[n]->copyInto(&(ng->antennas[3]));
		//		sub_trees.push_back(ng);
		//	  }
		//	}
		//	  }
		//	}
		//  }
		//  
		//  all_antennas.clear();
		//}


		////duplicate a stand-alone copy of self into 
		////the destination
		//void NGTree::copyInto(NGTree *destination)
		//{
		//  assert ( destination != NULL );

		//  destination->with_fucose = this->with_fucose;

		//  if ( this->root == NULL )
		//	destination->root = NULL;
		//  else
		//	destination->root = cloneNGTree(this->root);
		//}


		//Monosaccharide* NGTree::cloneNGTree(Monosaccharide* root)
		//{
		//  assert( root != NULL );
		//  Monosaccharide* result;

		//  result = new Monosaccharide();

		//  //copy the content 
		//  result->mono_type = root->mono_type;
		//  result->incoming_bond = root->incoming_bond;
		//  result->parent = root->parent;

		//  result->children.clear();
		//  
		//  for (int i=0; i< (int) root->children.size(); i++) {
		//	if ( root->children[i] != NULL ) {
		//	  Monosaccharide* temp = 
		//	cloneNGTree(root->children[i]);
		//	  temp->parent = result;
		//	  result->children.push_back(temp);
		//	}
		//  }

		//  return result;
		//}


		//void NGlycan::print() 
		//{
		//  cout << "==============start of tree" << endl;
		//  printNGTree(pentamer.root);
		//  if ( pentamer.with_fucose )
		//	cout << " fucose";

		//  cout << endl;

		//  for (int i=0; i< (int) antennas.size(); i++) {
		//	printNGTree( antennas[i].root );
		//	if ( antennas[i].with_fucose )
		//	  cout << " fucose";
		//	cout << endl;
		//  }
		//  cout << "==============end of tree" << endl;
		//}


		GlycanComposition::GlycanComposition()
		{
		  monosaccharide_count.clear();
		  monosaccharide_count.insert(pair<int, int>(MONOSACCHARIDE_HEXNAC, 0));
		  monosaccharide_count.insert(pair<int, int>(MONOSACCHARIDE_HEX, 0));
		  monosaccharide_count.insert(pair<int, int>(MONOSACCHARIDE_DEHEX, 0));
		  monosaccharide_count.insert(pair<int, int>(MONOSACCHARIDE_NEUAC, 0));
		  monosaccharide_count.insert(pair<int, int>(MONOSACCHARIDE_NEUGC, 0));
		  monosaccharide_count.insert(pair<int, int>(MONOSACCHARIDE_MOD_SULFATE, 0)); 
		  monosaccharide_count.insert(pair<int, int>(MONOSACCHARIDE_UNKNOWN, 0));
		  consider_sulfated = false ; 
		  glycan_type = ALL_NGLYCAN;
		  glycan_accurate_mass = 0 ; 
		  glycan_avg_mass = 0 ; 
		  is_decoy = false ; 
		}

		bool GlycanComposition::validate(GlycanComposition &min_gc,GlycanComposition &max_gc,int min_size, int max_size, int glycan_type) 
		{
				
		  int s = this->totalNumMonosacchrides();

		  if ( s < min_size || s > max_size )
			return false;


		  //check the cardinality specification of the glycan
		  if (! (checkCardinality(this->monosaccharide_count[MONOSACCHARIDE_HEX], 
			glycan_has_hexose) &&
			checkCardinality(this->monosaccharide_count[MONOSACCHARIDE_DEHEX], 
			glycan_has_fucose) &&
			checkCardinality(this->monosaccharide_count[MONOSACCHARIDE_HEXNAC], 
			glycan_has_hexnac) &&
			checkCardinality(this->monosaccharide_count[MONOSACCHARIDE_NEUAC], 
			glycan_has_neuac) && 
			checkCardinality(this->monosaccharide_count[MONOSACCHARIDE_NEUGC], 
			glycan_has_neugc) &&
			checkCardinality(this->monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE], 
			glycan_has_sulfate)))
			return false;

		  if ( this->totalNumMonosacchrides() <
			   min_gc.monosaccharide_count[MONOSACCHARIDE_TOTAL] ||
			   this->totalNumMonosacchrides() >
			   max_gc.monosaccharide_count[MONOSACCHARIDE_TOTAL] ||
			   this->monosaccharide_count[MONOSACCHARIDE_HEX] <
			   min_gc.monosaccharide_count[MONOSACCHARIDE_HEX] ||
			   this->monosaccharide_count[MONOSACCHARIDE_DEHEX] <
			   min_gc.monosaccharide_count[MONOSACCHARIDE_DEHEX] ||
			   this->monosaccharide_count[MONOSACCHARIDE_HEXNAC] <
			   min_gc.monosaccharide_count[MONOSACCHARIDE_HEXNAC] ||
			   this->monosaccharide_count[MONOSACCHARIDE_NEUAC] <
			   min_gc.monosaccharide_count[MONOSACCHARIDE_NEUAC] ||
			   this->monosaccharide_count[MONOSACCHARIDE_HEX] >
			   max_gc.monosaccharide_count[MONOSACCHARIDE_HEX] ||
			   this->monosaccharide_count[MONOSACCHARIDE_DEHEX] >
			   max_gc.monosaccharide_count[MONOSACCHARIDE_DEHEX] ||
			   this->monosaccharide_count[MONOSACCHARIDE_HEXNAC] >
			   max_gc.monosaccharide_count[MONOSACCHARIDE_HEXNAC] ||
			   this->monosaccharide_count[MONOSACCHARIDE_NEUAC] >
			   max_gc.monosaccharide_count[MONOSACCHARIDE_NEUAC] ||
			   this->monosaccharide_count[MONOSACCHARIDE_NEUGC] >
			   max_gc.monosaccharide_count[MONOSACCHARIDE_NEUGC] ) {
			return false;
		  }

		  
		  if ( ! (glycan_type == ALL_NGLYCAN || 
			  this->glycan_type == glycan_type) ) 
		  {
			return false;
		  }

		  return true;
		}


		void GlycanComposition::print() 
		{
		  if ( this->glycan_type == COMPLEX_GLYCAN )
			cout << "Complex ";
		  else if ( this->glycan_type == HIGH_MANNOSE_GLYCAN )
			cout << "Highman ";
		  else if ( this->glycan_type == HYBRID_GLYCAN )
			cout << "Hybrid ";

		  if (!consider_sulfated)
		  {


		  cout << "Composition = { HEXNAC= " 
			   << monosaccharide_count[MONOSACCHARIDE_HEXNAC] << ", HEX= "
			   << monosaccharide_count[MONOSACCHARIDE_HEX] << ", DEHEX= "
			   << monosaccharide_count[MONOSACCHARIDE_DEHEX] << ", NEUAC= "
			   << monosaccharide_count[MONOSACCHARIDE_NEUAC] << ", NEUGC= "
			   << monosaccharide_count[MONOSACCHARIDE_NEUGC] << ", UNKNOWN= "
			   << monosaccharide_count[MONOSACCHARIDE_UNKNOWN] << "}" 
			   << " ,accurate_mass = " << this->accurateMass() 
			   << " ,average_mass = " << this->averageMass() << endl;
		  }
		  else
		  {
			   //to do
		  }
		}

		//get a unique ID for each composition
		int GlycanComposition::getID() 
		{
		  int n = 10;
		  int id = 0;
		  
		  id += monosaccharide_count[MONOSACCHARIDE_HEXNAC];
		  id += monosaccharide_count[MONOSACCHARIDE_HEX] * n;
		  id += monosaccharide_count[MONOSACCHARIDE_DEHEX] * n * n;
		  id += monosaccharide_count[MONOSACCHARIDE_NEUAC] * n * n * n;
		  id += monosaccharide_count[MONOSACCHARIDE_NEUGC] * n * n * n * n;
		  id += monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE] * n * n * n * n * n; 
		  id += glycan_type * n * n * n * n * n * n;

		  return id;
		}

		//return the total number of monosaccharides in the glycan composition
		int GlycanComposition::totalNumMonosacchrides() 
		{
		  int total = 0;
		  
		  total+= monosaccharide_count[MONOSACCHARIDE_HEXNAC];
		  total += monosaccharide_count[MONOSACCHARIDE_HEX];
		  total += monosaccharide_count[MONOSACCHARIDE_DEHEX];
		  total += monosaccharide_count[MONOSACCHARIDE_NEUAC];
		  total += monosaccharide_count[MONOSACCHARIDE_NEUGC];

		  assert( total > 0 );
		  return total;
		}

	


		//returns the string form of the composition
		string GlycanComposition::getCompositionString() 
		{
		  //string result = string("");
		  char tmp[MAX_LINE];
		
		  if (!consider_sulfated)
		  {
		  sprintf(tmp, "HexNAc * %d + Hex * %d + DeHex * %d + NeuAc * %d + NeuGc * %d", 
			  monosaccharide_count[MONOSACCHARIDE_HEXNAC],
			  monosaccharide_count[MONOSACCHARIDE_HEX],
			  monosaccharide_count[MONOSACCHARIDE_DEHEX],
			  monosaccharide_count[MONOSACCHARIDE_NEUAC],
			  monosaccharide_count[MONOSACCHARIDE_NEUGC]);
		  }
		  else
		  {
			  sprintf(tmp, "HexNAc * %d + Hex * %d + DeHex * %d + NeuAc * %d + SO3 * %d + NeuGc * %d ", 
			  monosaccharide_count[MONOSACCHARIDE_HEXNAC],
			  monosaccharide_count[MONOSACCHARIDE_HEX],
			  monosaccharide_count[MONOSACCHARIDE_DEHEX],
			  monosaccharide_count[MONOSACCHARIDE_NEUAC],
			  monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE],
			  monosaccharide_count[MONOSACCHARIDE_NEUGC]);
		  }
		  return string(tmp);
		}

		//Anoop:parses a composition string and sets the count of each monosaccharide
		void GlycanComposition::parseCompositionString(string gc_comp)
		{
			size_t found1 = -1;
			size_t found2;
			size_t found3; 
			size_t t_found ; 
			size_t decoy_found ; 
			string sugar;
			string mono ; 
			string mono_count ; 

			found1 = gc_comp.find("Decoy", 0) ; 
			if (found1 == string::npos)
			{
				found1 = -1 ; 
				is_decoy = false ; 
			}
			else
			{
				is_decoy = true ; 
			}
				

			found2 = gc_comp.find("+", found1+1) ; 
			while (found2 != string::npos)
			{			
				sugar  = gc_comp.substr(found1+1, found2-found1-1);
				found3 = sugar.find("*") ; 
				mono = sugar.substr(0, found3-1) ; 			
				mono_count = sugar.substr(found3+2, sugar.length()-found3-3) ; 

				if (strcmp(mono.c_str(), "HexNAc") == 0)
				{
					monosaccharide_count[MONOSACCHARIDE_HEXNAC] = atoi(mono_count.c_str()) ; 
				}
				if (strcmp(mono.c_str(), "Hex") == 0)
				{
					monosaccharide_count[MONOSACCHARIDE_HEX] = atoi(mono_count.c_str()) ; 
				}
				if (strcmp(mono.c_str(), "DeHex") == 0)
				{
					monosaccharide_count[MONOSACCHARIDE_DEHEX] = atoi(mono_count.c_str()) ; 
				}
				if (strcmp(mono.c_str(), "NeuAc") == 0)
				{
					monosaccharide_count[MONOSACCHARIDE_NEUAC] = atoi(mono_count.c_str()) ; 
				}		
				if (strcmp(mono.c_str(), "SO3") == 0)
				{
					monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE] = atoi(mono_count.c_str()); 
					consider_sulfated = true ; 
				}
				
				found1 = found2+1 ;
				found2 = gc_comp.find("+", found1) ; 				
			}

			sugar = gc_comp.substr(found1+1, gc_comp.length()-found1-1); 
			found3 = sugar.find("*") ; 
			mono = sugar.substr(0, found3-1) ; 			
			mono_count = sugar.substr(sugar.length()-1,1) ; 
			if (strcmp(mono.c_str(), "NeuGc") == 0)
			{
				monosaccharide_count[MONOSACCHARIDE_NEUGC] = atoi(mono_count.c_str()) ; 
			}
				
		}


		//returns the string form of the composition
		string GlycanComposition::getConvertableString()
		{
		  //string result = string("");
		  char tmp[MAX_LINE];

		  sprintf(tmp, "(HNAC%dHEX%dDEHX%dSACID%d)", 
			  monosaccharide_count[MONOSACCHARIDE_HEXNAC],
			  monosaccharide_count[MONOSACCHARIDE_HEX],
			  monosaccharide_count[MONOSACCHARIDE_DEHEX],
			  monosaccharide_count[MONOSACCHARIDE_NEUAC]);
		  return string(tmp);
		}

		

		bool GlycanComposition::hasNeuAc()
		{
			if (this->monosaccharide_count[MONOSACCHARIDE_NEUAC] == 0)
				this->glycan_has_neuac = 0; 
			else
				this->glycan_has_neuac = 1; 
			return (this->glycan_has_neuac) ; 
		}

		bool GlycanComposition::hasDeHex()
		{
			if (this->monosaccharide_count[MONOSACCHARIDE_DEHEX] ==0)
				this->glycan_has_fucose =0 ; 
			else
				this->glycan_has_fucose = 1; 
			return (this->glycan_has_fucose) ; 
	
		}

		


		float GlycanComposition::accurateMass()
		{
		  float mass = 0;

		  mass += this->monosaccharide_count[MONOSACCHARIDE_HEXNAC] *
			monosaccharide_mass_accurate[MONOSACCHARIDE_HEXNAC];

		  mass += this->monosaccharide_count[MONOSACCHARIDE_HEX] *
			monosaccharide_mass_accurate[MONOSACCHARIDE_HEX];

		  mass += this->monosaccharide_count[MONOSACCHARIDE_DEHEX] *
			monosaccharide_mass_accurate[MONOSACCHARIDE_DEHEX];

		  mass += this->monosaccharide_count[MONOSACCHARIDE_NEUAC] *
			monosaccharide_mass_accurate[MONOSACCHARIDE_NEUAC];

		  mass += this->monosaccharide_count[MONOSACCHARIDE_NEUGC] *
			monosaccharide_mass_accurate[MONOSACCHARIDE_NEUGC];
		  
		  if(consider_sulfated)
		  {
			  mass += this->monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE] *
				  monosaccharide_mass_accurate[MONOSACCHARIDE_MOD_SULFATE];
		  }

		  mass += H2O;

		
		  return mass;

		 
		}


		float GlycanComposition::averageMass() 
		{
		  float mass = 0;

		  mass += this->monosaccharide_count[MONOSACCHARIDE_HEXNAC] *
			monosaccharide_mass_average[MONOSACCHARIDE_HEXNAC];

		  mass += this->monosaccharide_count[MONOSACCHARIDE_HEX] *
			monosaccharide_mass_average[MONOSACCHARIDE_HEX];

		  mass += this->monosaccharide_count[MONOSACCHARIDE_DEHEX] *
			monosaccharide_mass_average[MONOSACCHARIDE_DEHEX];

		  mass += this->monosaccharide_count[MONOSACCHARIDE_NEUAC] *
			monosaccharide_mass_average[MONOSACCHARIDE_NEUAC];

		  mass += this->monosaccharide_count[MONOSACCHARIDE_NEUGC] *
			monosaccharide_mass_average[MONOSACCHARIDE_NEUGC];

		  if(consider_sulfated)
		  {
			  mass += this->monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE] *
				  monosaccharide_mass_average[MONOSACCHARIDE_MOD_SULFATE];
		  }

		  mass += H2O;
		  
		  return mass;
		}


		GlycanComposition GlycanComposition::combine(GlycanComposition &gc) 
		{
		  GlycanComposition new_gc;

		  new_gc.monosaccharide_count[MONOSACCHARIDE_HEXNAC]= 
			this->monosaccharide_count[MONOSACCHARIDE_HEXNAC] +
			gc.monosaccharide_count[MONOSACCHARIDE_HEXNAC];

		  new_gc.monosaccharide_count[MONOSACCHARIDE_HEX]= 
			this->monosaccharide_count[MONOSACCHARIDE_HEX] +
			gc.monosaccharide_count[MONOSACCHARIDE_HEX];

		  new_gc.monosaccharide_count[MONOSACCHARIDE_DEHEX]= 
			this->monosaccharide_count[MONOSACCHARIDE_DEHEX] +
			gc.monosaccharide_count[MONOSACCHARIDE_DEHEX];

		  new_gc.monosaccharide_count[MONOSACCHARIDE_NEUAC]= 
			this->monosaccharide_count[MONOSACCHARIDE_NEUAC] +
			gc.monosaccharide_count[MONOSACCHARIDE_NEUAC];

		  new_gc.monosaccharide_count[MONOSACCHARIDE_NEUGC]= 
			this->monosaccharide_count[MONOSACCHARIDE_NEUGC] +
			gc.monosaccharide_count[MONOSACCHARIDE_NEUGC];

		  new_gc.monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE] =
			  this->monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE] + 
			  gc.monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE] ; 

		  new_gc.monosaccharide_count[MONOSACCHARIDE_UNKNOWN]= 
			this->monosaccharide_count[MONOSACCHARIDE_UNKNOWN] +
			gc.monosaccharide_count[MONOSACCHARIDE_UNKNOWN];

		  new_gc.glycan_accurate_mass = new_gc.accurateMass() ; 
		  new_gc.glycan_avg_mass = new_gc.averageMass() ; 
		  return new_gc;
		}

		//some monosaccharides are interchangable in terms of functions,
		//such as NeuGc and NeuAc. this function generates all possible glycan
		//compositions with Gc and Ac. and check if they are valid under user-defined,
		//max and min composition restriction.
		//save the results in the isoform_comp
		//void NGlycan::getIsoformCompositions(vector<GlycanComposition> &isoform_comp, GlycanComposition min_composition, GlycanComposition max_composition)
		//{
		//	//GlycanComposition gc;

		// // gc = pentamer.getComposition().combine(gc);

		// // for (int i=0; i< (int) antennas.size(); i++) {
		// //   gc = antennas[i].getComposition().combine(gc);
		// // }

		// // gc.glycan_type = this->getGlycanType();

		//	GlycanComposition gc = this->getComposition();

		//  //replace each NeuAc by NeuGc, and come up with an different glycan composition.
		//  int num_neuac = gc.monosaccharide_count[MONOSACCHARIDE_NEUAC];
		//  isoform_comp.clear();
		//  for (int i=0; i <= gc.monosaccharide_count[MONOSACCHARIDE_NEUAC];
		//	  i++) 
		//  {
		//	  GlycanComposition isoform_composition;
		//	  gc.clone(isoform_composition);
		//	  isoform_composition.monosaccharide_count[MONOSACCHARIDE_NEUGC] = i;
		//	  isoform_composition.monosaccharide_count[MONOSACCHARIDE_NEUAC] = 
		//		  gc.monosaccharide_count[MONOSACCHARIDE_NEUAC] - i;
		//	  isoform_comp.push_back(isoform_composition);
		//  }
		//}

		//clone the content of current glycan composition into the
		//new gc.
		void GlycanComposition::clone(GlycanComposition &new_gc)
		{
			new_gc.monosaccharide_count[MONOSACCHARIDE_HEXNAC]= 
			this->monosaccharide_count[MONOSACCHARIDE_HEXNAC];

		  new_gc.monosaccharide_count[MONOSACCHARIDE_HEX]= 
			this->monosaccharide_count[MONOSACCHARIDE_HEX];

		  new_gc.monosaccharide_count[MONOSACCHARIDE_DEHEX]= 
			this->monosaccharide_count[MONOSACCHARIDE_DEHEX];

		  new_gc.monosaccharide_count[MONOSACCHARIDE_NEUAC]= 
			this->monosaccharide_count[MONOSACCHARIDE_NEUAC];

		  new_gc.monosaccharide_count[MONOSACCHARIDE_NEUGC]= 
			this->monosaccharide_count[MONOSACCHARIDE_NEUGC];

		  new_gc.monosaccharide_count[MONOSACCHARIDE_UNKNOWN]= 
			this->monosaccharide_count[MONOSACCHARIDE_UNKNOWN];
		  
		  new_gc.monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE] = 
			  this->monosaccharide_count[MONOSACCHARIDE_MOD_SULFATE] ; 			

		  new_gc.glycan_type = this->glycan_type;
		  new_gc.glycan_accurate_mass = this->glycan_accurate_mass ; 
		  new_gc.glycan_avg_mass = this->glycan_avg_mass ; 
		}

		//GlycanComposition NGlycan::getComposition()
		//{
		//  GlycanComposition gc;

		//  gc = pentamer.getComposition().combine(gc);

		//  for (int i=0; i< (int) antennas.size(); i++) {
		//	gc = antennas[i].getComposition().combine(gc);
		//  }

		//  gc.glycan_type = this->getGlycanType();
		//  return gc;
		//}

		//int NGlycan::getGlycanType()
		//{
		//	return COMPLEX_GLYCAN;
		//}


		//GlycanComposition NGTree::getComposition() 
		//{
		//  GlycanComposition gc;

		//  getComposition(this->root, gc);

		//  if ( this->with_fucose && this->root != NULL )
		//	gc.monosaccharide_count[MONOSACCHARIDE_DEHEX]++;

		//  return gc;
		//}



		//void NGTree::getComposition(Monosaccharide* root, GlycanComposition &gc)
		//{
		//  if ( root == NULL )
		//	return;
		//  else 
		//  {
		//		for (int i=0; i< (int) root->children.size(); i++) 
		//		{
		//			getComposition(root->children[i], gc);
		//		}

		//		assert ( root->mono_type == MONOSACCHARIDE_DEHEX ||
		//		 root->mono_type == MONOSACCHARIDE_HEXNAC ||
		//		 root->mono_type == MONOSACCHARIDE_HEX ||
		//		 root->mono_type == MONOSACCHARIDE_NEUAC ||
		//		 root->mono_type == MONOSACCHARIDE_NEUGC ||
		//		 root->mono_type == MONOSACCHARIDE_UNKNOWN );
		//		gc.monosaccharide_count[root->mono_type]++;
		//	}
		//}


		////check if the current glycan is valid or not
		//bool NGlycan::validateGlycan(GlycanComposition &min_gc, GlycanComposition &max_gc, int max_antenna_fucose, int max_pentamer_fucose, int min_size, int max_size) 
		//{
		//  int s = this->size();

		//  if ( s < min_size || s > max_size )
		//	return false;

		//  //check the composition constrain
		//  GlycanComposition gc = this->getComposition();

		//  //check the fucose on the pentamer
		//  GlycanComposition pentamer_composition = 
		//	this->pentamer.getComposition();

		//  if ( pentamer_composition.monosaccharide_count[MONOSACCHARIDE_DEHEX] >
		//	   max_pentamer_fucose ) {
		//	return false;
		//  }

		//  //check number of fuocsed antennas
		//  int num_fucosed_antenna = 0;
		//  for ( int i=0; i< (int) antennas.size(); i++ ) {
		//	if ( antennas[i].size() > 0 && antennas[i].with_fucose )
		//	  num_fucosed_antenna++;
		//  }

		//  if ( num_fucosed_antenna > max_antenna_fucose )
		//	return false;

		//  /*
		//  //check the cardinality specification of the glycan
		//  if (! (checkCardinality(gc.monosaccharide_count[MONOSACCHARIDE_HEX], 
		//	glycan_has_hexose) &&
		//	checkCardinality(gc.monosaccharide_count[MONOSACCHARIDE_DEHEX], 
		//	glycan_has_fucose) &&
		//	checkCardinality(gc.monosaccharide_count[MONOSACCHARIDE_HEXNAC], 
		//	glycan_has_hexnac) &&
		//	checkCardinality(gc.monosaccharide_count[MONOSACCHARIDE_NEUAC], 
		//	glycan_has_neuac) && 
		//	checkCardinality(gc.monosaccharide_count[MONOSACCHARIDE_NEUGC], 
		//	glycan_has_neugc)) )
		//	return false;

		//  if ( gc.totalNumMonosacchrides() <
		//	   min_gc.monosaccharide_count[MONOSACCHARIDE_TOTAL] ||
		//	   gc.totalNumMonosacchrides() >
		//	   max_gc.monosaccharide_count[MONOSACCHARIDE_TOTAL] ||
		//	   gc.monosaccharide_count[MONOSACCHARIDE_HEX] <
		//	   min_gc.monosaccharide_count[MONOSACCHARIDE_HEX] ||
		//	   gc.monosaccharide_count[MONOSACCHARIDE_DEHEX] <
		//	   min_gc.monosaccharide_count[MONOSACCHARIDE_DEHEX] ||
		//	   gc.monosaccharide_count[MONOSACCHARIDE_HEXNAC] <
		//	   min_gc.monosaccharide_count[MONOSACCHARIDE_HEXNAC] ||
		//	   gc.monosaccharide_count[MONOSACCHARIDE_NEUAC] <
		//	   min_gc.monosaccharide_count[MONOSACCHARIDE_NEUAC] ||
		//	   gc.monosaccharide_count[MONOSACCHARIDE_HEX] >
		//	   max_gc.monosaccharide_count[MONOSACCHARIDE_HEX] ||
		//	   gc.monosaccharide_count[MONOSACCHARIDE_DEHEX] >
		//	   max_gc.monosaccharide_count[MONOSACCHARIDE_DEHEX] ||
		//	   gc.monosaccharide_count[MONOSACCHARIDE_HEXNAC] >
		//	   max_gc.monosaccharide_count[MONOSACCHARIDE_HEXNAC] ||
		//	   gc.monosaccharide_count[MONOSACCHARIDE_NEUAC] >
		//	   max_gc.monosaccharide_count[MONOSACCHARIDE_NEUAC] ||
		//	   gc.monosaccharide_count[MONOSACCHARIDE_NEUGC] >
		//	   max_gc.monosaccharide_count[MONOSACCHARIDE_NEUGC] ) {
		//	return false;
		//  }
		//  */

		//  return true;
		//}


		////return the number of monosaccharides in the glycan
		//int NGlycan::size() 
		//{
		//  int s = pentamer.size();

		//  for ( int i=0; i< (int) antennas.size(); i++ ) {
		//	s += antennas[i].size();
		//  }

		//  return s;
		//}


		////return the number of monosaccharides in the glycan
		//int NGTree::size() 
		//{
		//  int s = 0;

		//  if ( this->root != NULL && this->with_fucose )
		//	s = 1; //count the fucose

		//  s += getNGTreeSize(this->root);
		//  return s;
		//}

		//int NGTree::getNGTreeSize(Monosaccharide *root) {
		//  int s = 0;

		//  if ( root == NULL )
		//	s = 0;
		//  else {
		//	s = 1;

		//	for (int i=0; i< (int) root->children.size(); i++) {
		//	  s += getNGTreeSize(root->children[i]);
		//	}
		//  }

		//  return s;
		//}

		////this method is used in parsing the command line input.
		////Depending on the input parameter, it returns whether a
		////glycan should have some certain number of monosaccharides
		////or not. it returns 4 values: POSSBILE, YES, NO and UNKNOWN
		//int strToMonosaccharideCardinality(string str) 
		//{
		//	int cardinality = GLYCAN_HAS_THIS_MONOSACCHARIDE_UNKNOWN;

		//  if (str.find("possible", 0) != string::npos ) 
		//	cardinality = GLYCAN_HAS_THIS_MONOSACCHARIDE_POSSIBLE;
		//  else if (str.find("yes", 0) != string::npos ) 
		//	cardinality = GLYCAN_HAS_THIS_MONOSACCHARIDE_YES;
		//  else if (str.find("no", 0) != string::npos ) 
		//	cardinality = GLYCAN_HAS_THIS_MONOSACCHARIDE_NO;
		//  else {
		//	cardinality = GLYCAN_HAS_THIS_MONOSACCHARIDE_UNKNOWN;
		//  }

		//  return cardinality;
		//}

		////check if the given number n obeys the given cardinality
		bool checkCardinality(int n, int cardinality)
		{
		  if ( n < 0 )
			return false;
		  else if ( n == 0 && 
					cardinality == GLYCAN_HAS_THIS_MONOSACCHARIDE_YES )
			return false;
		  else if ( n > 0 && 
					cardinality == GLYCAN_HAS_THIS_MONOSACCHARIDE_NO )
			return false;		
		  else if ( cardinality != GLYCAN_HAS_THIS_MONOSACCHARIDE_NO &&
			cardinality != GLYCAN_HAS_THIS_MONOSACCHARIDE_YES &&
			cardinality != GLYCAN_HAS_THIS_MONOSACCHARIDE_POSSIBLE ||
			cardinality == GLYCAN_HAS_THIS_MONOSACCHARIDE_UNKNOWN )
			return false;
		  else return true;
		}
	}
}