// Written by Yin Wu (wuyin@indiana.edu)
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// modified by Anoop for GlypID

#ifndef ISOFORM_H
#define ISOFORM_H

#include "Peak.h"
#include "Spectrum.h"
#include "debug.h"

namespace Engine
{
	namespace IsoformProcessor
	{

		//the maximum number of peaks that an envelope can have
		#define MAX_ENV_SIZE 35
		#define MAX_ISOFORM_WIDTH 20.0
		#define FIRST_MONO_TO_HIGHEST_MONO_RATIO_BOUND 0.1
		#define ISO_ID_START 1

		class Peak;
		class Spectrum;

		class Isoform 
		{
			private: 
			  //a integer used in the glycan cluster sreening method
			  //indicating with the current ion was ever in the screening
			  //inclusion list or not. 
			  //the value of the screening_history is the File ID of 
			  //the raw file in which the ion is in the screening 
			  //inclusion list. If this screening_history equals the 
			  //next file ID, that means the Ion should be put in the
			  //next screening inclusion list. and the screening_flag
			  //should be TRUE.
			  bool glycan_cluster_screening_flag;
			  int glycan_cluster_screening_history;

			  vector<float> glycan_cluster_screening_dependency;
			  float glycan_cluster_screening_is_approved_by_ion;

		public:
			  int id; //id starts at 1
			  int charge;
			  float priority, priority_2nd, priority_3rd;

			  float ms2_score;
			  int start_scan, end_scan, size, last_mono_index;
			  float start_time, end_time; //the time interval in which the isoform appears
			  float total_intensity; //sum of intensity of the envelopes in different MS1 scans.
			  float top_intensity; //the highest aggregated intensity of the envelope among envelopes in different MS1 scans.
			  float single_peak_top_intensity; //the highest single peak intensity of the envelope among envelopes in different MS1 scans.
			  float first_mono; //the m/z of the first mono peak
			  float flag; //used by programmer.
			  string glycan_cluster_str;

			  Peak envelope[MAX_ENV_SIZE];
			  vector<Isoform> siblings; //Isoforms that belong to the sample ion
			  vector<Spectrum*> attached_sub_ms;
			  vector<float> ion_merging_history; //record which ions are merged with this ion 
												 //during inter-experiment data comparison. a 
												 //float number with ion_id.file_id
												 //represents an ion
			  vector<int> glycan_cluster_info;

			  float key_sub_ms_header;
			  vector<string> annotations;
			  int file_id_in_which_glycopeptide_confirmed;
			  int file_id_in_which_ion_is_found;

			  void init(); //init the class variables
			  Isoform();
			  Isoform(int id, int charge, float first_mono,
				  float start_time, float end_time,
				  int start_scan, int end_scan);
			  Isoform(int id, int charge, float first_mono,
				  float start_time, float end_time,
				  int start_scan, int end_scan,
				  vector<Peak>& plist);
			  ~Isoform();

			  //recalculate the mz of first mono and
			  //total intensity.
			  //according to the peaks in the envelope 
			  //(take their averages). it does not check
			  //for any possible error or mis-match.
			  //the ratio_bound is used to eliminate 
			  //false first mono
			  void update(float first_mono_to_highest_mono_ratio_bound);

			  //used to track which screening ion does the current ion
			  //depend on. if the current ion depends on some other
			  //ion, then it is pending to be approved to the inclusion
			  //list. if it does not depend on any ion, it is either
			  //in the inclusion list, or it is rejected already.
			  //the values are the ion ID of the depended ions.
			  void updateScreeningDependency(vector<float> &screen_ion_id_list);
			  void addScreeningDependency(vector<float> &screen_ion_id_list);
			  vector<float>& getGlycanClusterScreeningDependency();
			  void clearGlycanClusterScreeningDependency();

			  void updateFileID(int id);

			  //check if current isoform overlaps with iso
			  bool overlap(Isoform& iso, float err);
			  bool overlap(Isoform& iso, int ppm_err);
			  bool overlap(Isoform &iso, int ppm_err, float elution_time_shift);

			  bool isCompleteIsoform();
			  bool isIdentifiedAsGlycopeptide();

			  //verify if, in the glycan cluster screening method, it 
			  //has been approved to be put in inclusion list.
			  bool isApprovedInScreening();
			  //set the ion to be "been approved"
			  void setApprovedInScreening(int ion_id, 
				  int data_file_id);
			  float getApprovedInScreeningBy();
			  int getGlycanClusteringScreeningHistory();
			  void setGlycanClusteringScreeningHistory(int file_id);


			  //combine() and finalize() are used together
			  //to combine sibling Isoforms from multiple scans.
			  void combine(Isoform& iso);
			  void finalize();

			  void print();

			  void addAnnotation(string &annotation);

			  //test if the precursor mz is for this isoform or not
			  bool acceptSubMS(Spectrum* ms2);

			  //test if the isoform is covered by the ion acquiziont
			  //window of the MS2
			  bool isCoveredBySubMS(Spectrum* ms2, float mz_tolerance);

			  //test if an peak (m/z) is part of the current
			  //isoform. given (-err, +err) window
			  bool isPartOf(float mz, float err);

			  //convert the current isoform to the peak whose
			  // m/z is the first mono.
			  Peak toFirstMonoPeak();
			  Peak toMissedFirstMonoPeak();

			  //returns the m/z of the last mono peak.
			  float getLastMono();

			  //returns the index (starting from 1) 
			  //of the last mono peak.
			  int getLastMonoIndex();
			  
			  //returns the m/z of the most intense peak.
			  float getMostIntenseMono();

			  //returns the intensity of the most intense peak
			  float getIntensityOfMostIntenseMono();

			  //add a new sub MS(2) to the isoform. the added
			  //sub ms are sorted by their ms2 score.
			  void addSubMS(Spectrum *ms);

			  //find the best glyco-MS2
			  Spectrum* getKeySubMS();

			  //returns the m/z of the ion, which can be
			  //put into the inclusion list for MS/MS scan
			  //for the envolope
			  float getIncListIon();

			  vector<int> getMsMsPerExperimentHistory();

			  //check if the current ion is in the screening
			  //ion list for raw of the the data_file_id.
			  //e.g. if ion is put into the screening
			  //ion list, and is subject to MS/MS in the experiment
			  //which generate file i. then
			  // ion.isInScreenIonList(i) should return true.
			  bool isInScreenIonList(int data_file_id);
			  void setInScreenIonList(int data_file_id);

			  void setSinglePeakTopIntensity();

			  //clone a copy of self in a new block of memory.
			  //since Isoform uses an arrary of points of class
			  //Peak. this is the only safe way to copy an Isoform
			  //to another.
			  /*
			  Isoform* clone() {
				Isoform* iso = new iso();
				(*iso) = (*this);
				for (int i=0; i < MAX_ENV_SiZE; i++) {
				  if ( this->envelop[i] != NULL ) {
				iso->envelope[i] = new Peak(0, 0);
				(*(iso->envelope[i])) = (*(this->envelope[i]));
				  }
				}
			  };
			  */

			  //deIsotop an MS spectrum. this spectrum 
			  //must be high resuolution scan.
			  //return a vector of Isoforms
			  static void deIsotop(vector<Peak>& plist, float start_time, 
						   float end_time, int start_scan, 
						   int end_scan, int ppm, int max_missing,
						   vector<Isoform>& isolist);

		
		};
	}
}


#endif 
