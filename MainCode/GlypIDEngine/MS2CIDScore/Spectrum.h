///*****************************************************************************\
// * Spectrum.h
// *
// * File reader/writer
// *
// * author: Yin Wu (wuyin@indiana.edu)
// *
//\*****************************************************************************/
//
//#ifndef SPECTRUM_H
//#define SPECTRUM_H
//
//#include <stdlib.h>
//#include <iostream>
//#include <list>
//#include <iterator>
//#include <algorithm>
//#include <string>
//#include <math.h>
//#include <assert.h>
//using namespace std;
//
//#include "Peak.h"
//#include "Isoform.h"
//
//class Isoform;
//class Peak;
//

#include <algorithm>
#include <assert.h>
#include <vector>
#include "../PeakProcessor/PeakProcessor.h"
#include "../PeakProcessor/PeakData.h"
#include "../PeakProcessor/CID_Peak.h"
#include "Distribution.h" 

using namespace std;

namespace Engine
{
	namespace MS2CIDScoring
	{
		class Spectrum 
		{
		private:
			static const int num_monosaccharides = 4; //length of the array above
			static float monosaccharide_masses[num_monosaccharides];
			static const int num_oxonium_ions = 5; 
			static float oxonium_ions[num_oxonium_ions] ; 

			double mdbl_min_mz ; 
			double mdbl_max_mz ; 
			

			int mint_start_pos;
			int mint_end_pos ; 

			std::vector<Engine::PeakProcessing::CID_Peak> mobj_score_path ; 
			std::vector<Engine::PeakProcessing::CID_Peak> peaks ;
			int *l_table ;

			Engine::PeakProcessing::SORT_MODE  m_sort_mode ; 
			Engine::PeakProcessing::SORT_ORDER m_sort_order ;
			Engine::MS2CIDScoring::Distribution mobj_ms2_distribution ; 


		public:

			Spectrum(void) ; 
			~Spectrum(void) ; 

			float getMs2Score();

			void setSortOrder(Engine::PeakProcessing::SORT_ORDER order) ; 
			void setSortMode(Engine::PeakProcessing::SORT_MODE mode) ; 

			float getMs2ScoreByBackboneAndGlycan(float backbone_mass, float glycan_mass, int glycan_charge);

			//pre-compute the score for the ms2 spectrum.
			float computeMs2Score(vector<float> &mass_dif, Engine::PeakProcessing::PeakData &pk_data, int charge) ; 

			// Check if oxonium ion is present
			bool containsOxoniumion(Engine::PeakProcessing::PeakData &peak_data) ; 
			
		

			//The recurrence function that does the actual
			//enumeration for allMassDif().
			//Note that: The init value of "mass" should be 0.
			void enumerateMassDif(vector<float> &list, 	float mass, int mass_idx, int allowed_missing);

			//gives a list of all possible featured mass differences
			//betwee the glycan chains, depending on how many missing
			//peaks allowed. this list is derived from the
			//pre-defined (or pre-loeaded) possible masses
			void allMassDif(vector<float> &m_list, int allowed_missing);

			void initTable(int *tbl, int size) ; 
			
			// Initialize p value table
			void sampleMS2Distribution() ; 

			// Function to calculate p value
			float MS2ScoreToPValue(int score); 


			float approxFindWithCharge(vector<float> &list, float value, int charge, float error_tolerance) ; 
			void updateLinkTable(int *tbl, int size,int pk1, int pk2) ; 
			int getLink(int *tbl, int size, int pk1, int pk2) ; 
			void setLink(int *tbl, int size, int pk1, int pk2, int val) ; 

			//Converts PeakData peaks into yin's structure
			void ConvertPeaksToCIDPeaks(Engine::PeakProcessing::PeakData &peak_data); 
			bool addPeak(Engine::PeakProcessing::CID_Peak &p) ; 
			//bool peakSmaller(Engine::PeakProcessing::CID_Peak p1, Engine::PeakProcessing::CID_Peak p2) ;
			//bool peakGreater(Engine::PeakProcessing::CID_Peak p1, Engine::PeakProcessing::CID_Peak p2) ; 
			void topNPeaks(int n)  ; 
			void rankPeaksByIntensity() ; 
			void sortPeaks(Engine::PeakProcessing::SORT_MODE mode, Engine::PeakProcessing::SORT_ORDER order) ;
			

		};
	}
}


//namespace Engine
//{
//	namespace MS2CIDScore
//	{
//
//		enum {MS_UNKNOWN = 10, MS_ONE, MS_TWO, MS_THREE, MS_SFRAG};
//
//		class Spectrum {
//		  
//		public:
//		  float precursor_mass, precursor_charge;
//		  vector<Peak> peaks;
//		  Peak most_intense_peak;
//		  float flag;
//
//		  vector<Isoform> isolist; //list of extracted isoforms. 
//								   //only for high resolution data.
//
//		  //used in the glycan experiment only
//		  int turn_num, scanheader_num, file_id;
//		  int parent_scan; //the scan header of the parent scan (the parent of
//						   //ms3 is ms2, the parent of ms2 is ms1
//		  int ms_type; //MS_ONE, MS_TWO, MS_THREE, or IN SOURCE FRAGMENTATION
//		  float start_time; //in unit of minute
//		  float scores[MAX_CHARGE];  //there is a score for 
//									 //each charge level
//		  vector<Peak> score_path[MAX_CHARGE];
//
//		  //used by MS2 only
//		  float glycan_set_start[MAX_CHARGE]; //match to scores[]
//		  float total_intensity;
//		  int *l_table;
//		  int start_pos;
//		  bool oxonium_ion, is_valid_ms2;
//		  float oxonion_ion;
//
//		  Spectrum(float mass, float charge); 
//		  Spectrum(int turn, int scanheader, 
//			   float mass, float charge); 
//		  Spectrum();
//		  ~Spectrum();
//		  void print();
//		  bool addPeak(float mass, float intensity);
//		  bool addPeak(float mass, float intensity, int charge);
//		  bool addPeak(Peak &p);
//		  void normalize();
//		  Spectrum* clone();
//		  void sortPeaks(int mode, int order);
//		  void rankPeaksByIntensity();
//		  void topNPeaks(int n);
//		  void trimMS2(int n, float err);
//		  float totalIntensity();
//		  float totalScore();
//		  float totalScore(float backbone_mass,
//				   float glycan_mass, 
//				   int glycan_charge);
//		  float totalScore2(float backbone_mass,
//					float glycan_mass, 
//					int glycan_charge);
//
//		  void findCentroids();
//		  void bottomNPercent(float ratio);
//		  void grassLevel(float grass);  
//		  float getBinIntensity(float mz, float bin_size);
//		  float score(vector<float> &mass_dif, 
//				  int charge);
//		  void printLinks(int start);
//		  void refine(float window);
//		  void refineByGreatestPeak(float window);
//		  void refineConvolution();
//		  void capIntensity(float cap);
//		  void massInterval(float start, float end);
//		  Peak* pepmassMatchMS2(float error, 
//					float pepmass,
//					int charge);
//		  void backup();
//		  void restore();
//		  float peakVariance();
//		  bool findOxoniumIon();
//		  bool validate(int top_n_peaks);
//		  const vector<Peak>& getBackup() const;
//		  void printTable(int *tbl, int size);
//		  void peak_picking();
//		  //void deIsotop(float err);
//		  void findPeaks(vector<Peak> &candidate_peaks,
//				 vector<Peak> &detected_peaks,
//				 float err);
//		  void freeMemory();
//
//		  //test if the current spectrum has the peak of mz
//		  //in window (-err, +err)
//		  bool hasPeak(float mz, float err); 
//		  bool isIntensityRanked();
//		  int peakRank(float mz, float err);  
//		  int peakRank(float mz, float err,
//				   float mz_lower_bound,
//				   float mz_upper_bound);
//
//		  float hasBIon( float backbone );
//		  float getSupportedGlycopepBackbone();
//		  void setSupportedGlycopepBackbone(float backbone);
//		  void trimBaseLineForIsoform(float base_line_intensity);
//		  
//		  static void allMassDif(vector<float> &list, 
//					 int allowed_missing);
//		  static Spectrum* combineSpectrums(vector<Spectrum*> &specs);
//
//		private:
//		  float supported_glycopeptide_backbone;
//		  bool is_intensity_ranked; 
//		  vector<Peak> peaks_backup; //a backup of the peaks.
//									 //needed when complex operation 
//									 //is performed. used by calling 
//									 //backup and restore.
//
//		  void initTable(int *tbl, int size);
//
//		  int getLink(int *tbl, int size,
//				  int pk1, int pk2) {
//			assert ( pk1 < size && pk2 < size );
//			if ( pk1 <= pk2 )
//			  return tbl[(pk2*(pk2+1))/2+pk1];
//			else 
//			  return tbl[(pk1*(pk1+1))/2+pk2];    
//		  };
//
//
//		  void setLink(int *tbl, int size,
//				   int pk1, int pk2, int val) {
//			assert ( pk1 < size && pk2 < size );
//			if ( pk1 <= pk2 )
//			  tbl[(pk2*(pk2+1))/2+pk1] = val;
//			else 
//			  tbl[(pk1*(pk1+1))/2+pk2] = val;    
//		  };
//
//		  void updateLinkTable(int *tbl, int size, 
//					  int pk1, int pk2);
//		  float getBinIntensity(int index, 
//					float bin_size);
//		  static void enumerateMassDif(vector<float> &list,
//						   float mass, int mass_idx,
//						   int allowed_missing);
//		};
//
//
//		class Link {
//		public:
//		  int start, end;
//		  float length;
//
//		  Link(int s, int e, float l) {
//			start = s;
//			end = e;
//			length = l;
//		  };
//
//		};
//
//		class LinkTable {
//		public:
//		  void clear() {
//			l_table.clear();
//		  };
//
//		  void update(int s, int e, float length);
//		  float getPath(int s, int e);
//
//		private:
//		  vector<Link> l_table;
//		};
//	}
//}
//
//#endif 
//
