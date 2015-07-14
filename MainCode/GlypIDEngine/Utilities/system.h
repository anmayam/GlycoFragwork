///*****************************************************************************\
// * system.h
// *
// * Are the global variables and functions are declared in this file.
// *
// * author: Yin Wu (wuyin@indiana.edu)
// *
//\*****************************************************************************/
//
///*************Operating system options ********/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
//#include <strings.h>
#include <assert.h>
#include <math.h>
using namespace std;

//#ifdef MWT
//#import "..\bin\MwtWindll.dll" 

//#include <qstringlist.h>

//#include "./include/omics/core/Sequence.h" 
//
//
////using namespace  proteomics;
//
//#include "debug.h"
//#include "constants.h"
//#include "Peak.h"
//#include "Distribution.h"
//#include "Spectrum.h"
//#include "Dataset.h"
//#include "GStrTok.h"
//#include "NglycanSpectrum.h"
//#include "..\GlycanTheoretical\Nglycan.h"
//
class NglycanSpectrum;
class NglycanPeak;
//
//
//
///***********************************************/
///*********** MACROS ****************************/
///***********************************************/
///******* FILE I/O ********/
#define MAX_LINE 8192 
//enum {RUN_HEADER_SECTION=0, SCAN_HEADER_SECTION, SCAN_DATA_SECTION} ;


//******** peptides *******/
enum { MANNOSE_ONLY = 50, ALL_MONOSACCHARIDES, 
       SELF_DEFINED_156, SELF_DEFINED_99, SELF_DEFINED_176,
       SELF_DEFINED_103, SELF_DEFINED_93, SELF_DEFINED_72,
       SELF_DEFINED_56, SELF_DEFINED_33, SELF_DEFINED_00, 
       SELF_DEFINED_01, SELF_DEFINED_02};
//

enum {ASCENDING, DESCENDING} ; 

//
//********* system commands *******/
#define TURNED_ON 1
#define TURNED_OFF 0
#define INDEX_START 1
#define TRUE 1
#define FALSE 0
//
enum {PREDICT= 200, SAMPLE, MATCH_PEPTIDE, 
      EXTRACT_SPECTRA, EXTRACT_MASTER_INC_LIST,
      SEARCH_IONLIST, EVALUATE_MS2 };
//
//enum {UNKNOWN_DATA_TYPE=250, FT_DATA, ION_TRAP_DATA};
enum {INC_LIST_SORT_BY_ION_INTENSITY = 300, 
      INC_LIST_SORT_BY_GLYCO_CLUSTER_SIZE };
enum {INC_LIST_FULL_DYNAMIC = 350, INC_LIST_SEMI_DYNAMIC, INC_LIST_DYNAMIC_INCREASE };
enum {INC_LIST_SEGMENT_EVENLY = 400, INC_LIST_SEGMENT_BY_TIME_PERIOD};

enum {ISOFORM_MATCH_BY_FIRST_MONO=450, ISOFORM_MATCH_BY_HIGHEST_MONO,
    ISOFORM_MATCH_BY_FIRST_AND_HIGHEST_MONO, 
    ISOFORM_MATCH_BY_FIRST_OR_HIGHEST_MONO};


#define FT_MZ_ERROR 0.015
#define FT_PPM_ERROR 5
#define MZ_ERROR 0.7
#define H2O 18.01056
#define INTRA_SPECTRUM_ION_TRAP_MZ_ERROR 0.4

#define IODOACEDAMIDE_MONO 57.02146 
#define IODOACEDAMIDE_AVG 57.0520 

#define EPSILON 0.000001
#define FLAG_ON 1
#define FLAG_OFF 0
#define PI 3.1415926



//
///***********************************************/
///***********Global Functions/Variables *********/
///***********************************************/
//#ifndef MAIN_CPP
//#define EXTERN extern
//EXTERN bool glycan_library_file_name_provided;
//EXTERN char in_filename[MAX_LINE], 
//  ion_list_filename[MAX_LINE],
//  out_filename[MAX_LINE], 
//  candidate_protein_filename[MAX_LINE],
//  score_distribution_filename[MAX_LINE],
//  exp_inc_list_filename[MAX_LINE],
//  ionlist_output_filename[MAX_LINE], 
//  local_directory[MAX_LINE], 
//  glycan_library_file_name[MAX_LINE];
//EXTERN string inclusion_list_output_directory;
//EXTERN int exp_inc_list_size, inclusion_list_sort_by, 
//  inc_list_segment_method, inc_list_priority_threshold, 
//  data_file_id, inclusion_list_mode, inc_list_num_segments,
//  glycan_cluster_screening;
//EXTERN int glycopeptide_matching_glycan_type, 
//  glycopeptide_matching_max_fucose_on_pentamer, 
//  glycopeptide_matching_max_fucose_on_antena;
//EXTERN int isoform_match_method;
//EXTERN float  ms2_score_threshold; //the threshold above which, the ion is considered as glycopeptide
//EXTERN float intensity_cap;
//EXTERN int ft_ppm_error;
//EXTERN int max_tryptic_mis_cleavage;
//EXTERN float inc_list_segment_time_period, full_elution_time;
//EXTERN NglycanSpectrum ns;
//EXTERN GlycanComposition max_glycan_composition,
//  min_glycan_composition;
//EXTERN int glycan_has_hexose, glycan_has_fucose, 
//glycan_has_neuac, glycan_has_hexnac, glycan_has_neugc;
//EXTERN Distribution overall_score_dist;
//EXTERN Distribution ms2_score_dist;
//EXTERN Dataset data;
//EXTERN Dataset ms2_data;
//EXTERN Dataset ms3_data;
//EXTERN Dataset degly_data;
//EXTERN Dataset degly_ms2_data;
//EXTERN vector<Environment> environment_list;
//EXTERN vector<Sequence> candidate_protein_set;
//EXTERN float MH_masses_mannose_only[];
//EXTERN int num_MH_mannose_only;
//EXTERN float MH_masses_all_mono[];
//EXTERN int num_MH_all_mono;
//EXTERN float monosaccharide_masses[];
//EXTERN int num_monosaccharides;
//EXTERN float diagnostic_oxonium_ions[];			       
//EXTERN int num_oxonium_ions;
//EXTERN float pk_value;
//EXTERN float arg2;
//EXTERN float arg3;
//EXTERN float arg4;
//EXTERN float arg5;
//EXTERN float arg6;
//EXTERN float arg7;
//EXTERN float arg8;
//EXTERN float arg9;
//EXTERN float arg10;
//EXTERN float arg11;
//EXTERN float p_value_cutoff;
//EXTERN float peptide_elute_time_shift;
//EXTERN bool enforce_B_ion_check; /* If TRUE when matching an MS2 with the precursor ion in the MS1, a B-ion
//				  peak (backbone + GlcNAc) must be found in MS2. Otherwise, the MS2 will not
//				  be attched to the precursor.*/
//EXTERN bool is_score_distribution_file_available;
//EXTERN int data_type;
//EXTERN bool has_candidate_protein_set;
//EXTERN float ms1_mz_err, ms2_mz_err, 
//  ms3_mz_err, sfrag_mz_err;
//EXTERN float ms2_iso_width;  /*The width of MS1 isoform, when gathering MS2 
//			       ions. e.g. if the designated precursor of MS2 
//			       is mp, then the ions in range (mp-width, mp+width) 
//			       are gathered for fragmentation. */
//EXTERN float B_ion_intensity_ratio; /*the ratio between the intensity of the B-ion 
//				      and that of the most intense ion. Used to 
//				     verify that a discovery ion is strong enough
//				     to be a B-ion */
//EXTERN int command;
//
//EXTERN float self_defined156[];
//EXTERN int num_self_defined156;
//EXTERN float self_defined99[];
//EXTERN int num_self_defined99;
//EXTERN float self_defined176[];
//EXTERN int num_self_defined176;
//EXTERN float self_defined103[];
//EXTERN int num_self_defined103;
//EXTERN float self_defined93[];
//EXTERN int num_self_defined93;
//EXTERN float self_defined72[];
//EXTERN int num_self_defined72;
//EXTERN float self_defined56[];
//EXTERN int num_self_defined56;
//EXTERN float self_defined33[];
//EXTERN int num_self_defined33;
//EXTERN float self_defined_base[];
//EXTERN int num_self_defined_base;
//EXTERN map<int, GlycanComposition> glycan_compositions_in_search;
//EXTERN vector<GlycanComposition> external_glycan_compositions_list;
//EXTERN vector<GlycanComposition> computed_glycan_compositions_list;
//EXTERN vector<GlycanComposition> combined_glycan_compositions_list;
//EXTERN void init();
//#undef EXTERN
//#endif
//
//#ifdef IO_CPP
//#define EXTERN
//#else 
//#define EXTERN extern
//#endif
//EXTERN int max_line_len;
//EXTERN bool writeDTA(ofstream* f, Spectrum* spectrum, bool with_meta_data);
//EXTERN bool writeEIC(ofstream* f, 
//		     NglycanSpectrum* spectrum, 
//		     float threshold);
//EXTERN bool writeEIC4Excel(ofstream* f, 
//			   NglycanSpectrum* spectrum,
//			   float threshold);
//EXTERN void writeNglycanSpectra(ofstream *out_file,
//				vector<NglycanSpectrum> &nglycan_specs,
//				float threshold);
//EXTERN void writeNglycanSpectra4Excel(ofstream *out_file,
//				      vector<NglycanSpectrum> &nglycan_specs,
//				      float threshold);
//EXTERN void extract_spectra(char *input_filename);
//EXTERN bool writeEIC4ExcelHeader(ofstream* f);
//EXTERN void DtaNoComments(ifstream *in_f, ofstream *out_f);
//EXTERN Spectrum* readDTA(ifstream *f);
//EXTERN void reformatGlyFile(char* filename);
//EXTERN ifstream* openInputFile(const char* filename);
//EXTERN ofstream* openOutputFile(const char* filename);
//EXTERN void readCfg(ifstream *f);
//EXTERN int readIncList(char *file_name, 
//		       vector<Isoform> &inc_list,
//		       map<int, string> &data_file_id_list);
//EXTERN void writeIonList(char *file_name, 
//			 vector<Isoform> &inc_list,
//			 map<int, string> &data_file_id_list);
//EXTERN void writeExpIncList(char *file_name, 
//			    vector<Isoform> &inc_list,
//			    int inc_list_size,
//			    float mz_tolerance);
//EXTERN void writeExpIncListWithIonID(char *file_name, 
//		     vector<Isoform> &ion_list);
//EXTERN void readIncListWithIonID(char *exp_inc_list_with_ion_id_file, 
//      vector<int> &inclist_ion_id);
//
//EXTERN void readFastaFile(char *fasta_filename,
//			  vector<Sequence> &seq_list);
//EXTERN Spectrum* readEIC(ofstream* f);
//EXTERN Spectrum* readFTMS(ifstream *f);
//EXTERN void readScoreDistributionLog(ifstream *f, 
//				     vector<Environment>& env_list);
//EXTERN void loadGlyDatasets( const string& filename );
//EXTERN void loadGlycanList(vector<GlycanComposition> &glycan_compositions_list);
//EXTERN void computeGlycanList(vector<GlycanComposition> &glycan_compositions_list);
//EXTERN void combineGlycanList(vector<GlycanComposition> &gc_list1, 
//                              vector<GlycanComposition> &gc_list2,
//                              vector<GlycanComposition> &gc_list_result);
//EXTERN void printIsoformIntensityDist(vector<Isoform> &iso_list);
//EXTERN void writeScreeningInclist(char *file_name, vector<Isoform> &ion_list,
//                           int data_file_id, float mz_tolerance); 
//void writeDatasetLogFile(Dataset &data, 
//			 Dataset &ms2_data,
//			 vector<NglycanPeak>
//			 &refined_clusters);
//
//void writeGlycanClusterLogFile(char *file_name, 
//			       vector<NglycanPeak>
//			       &glycan_clusters,
//			       int data_file_id);
//void writeMS2LogFile(char *file_name, 
//		     Dataset &ms1_dataset,
//		     Dataset &ms2_dataset,
//		     int data_file_id);
//void writeSelectedGlycanClusterLogFile(char *file_name, 
//				       vector<NglycanPeak>
//				       &glycan_clusters,
//				       int data_file_id);
//void checkInclistCoverRate(vector<Isoform> &old_inc_list, 
//                           Dataset &ms2_data);
//#undef EXTERN
//
//
//#ifdef UTIL_CPP
//#define EXTERN
//#else 
//#define EXTERN extern
//#endif
//EXTERN float approxFindWithCharge(vector<float> &list, float value, 
//                           int charge, float error_tolerance);
//EXTERN float approxFind(vector<float> &list, float value);
//EXTERN float max(float v1, float v2);
//EXTERN bool approxEqual(float v1, float v2);
//EXTERN bool peakSmaller(Peak p1, Peak p2);
//EXTERN bool peakGreater(Peak p1, Peak p2);
//EXTERN int flip(int i);
////EXTERN int abs(int v);
//EXTERN void sort_array(float *src, float *dest, 
//		       int length, int order);
//EXTERN void SpectrumToPixelList(Spectrum *s, 
//				vector<GPixel> &plist);
//EXTERN int rangedRandom(int low, int high);
//EXTERN void removeDuplicates(vector<int> &t_list);
//EXTERN void removeDuplicates(vector<float> &t_list);
//EXTERN float idNumberToCombination(int id1, int id2);
//EXTERN int getID1FromCombination(float comb);
//EXTERN int getID2FromCombination(float comb);
//
////a function that computs the combination
////of N choose k
//EXTERN long NChooseK(int n, int k);
//EXTERN string number2Str(double num, 
//			 const char *format);
//#undef EXTERN
//
//
//#ifdef UI_CPP
//#define EXTERN
//#else 
//#define EXTERN extern
//#endif
//
//
//#undef EXTERN
//
//#endif 
