// Written by Yin Wu (wuyin@indiana.edu)
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// modified by Anoop for GlypID



#ifndef ISOFORM_CPP
#define ISOFORM_CPP
//
//#include <stdlib.h>
//#include <iostream>
//#include <vector>
//#include <map>
//#include <iterator>
//#include <math.h>
//#include <assert.h>
//using namespace std;
//
//#include "Isoform.h"
//#include "system.h"
//
//
//namespace Engine
//{
//	namespace IsoformProcessor
//	{
//
//		//init the class variables
//		void Isoform::init() 
//		{
//			  this->file_id_in_which_glycopeptide_confirmed = 0;
//			  this->file_id_in_which_ion_is_found = 0;
//			  this->key_sub_ms_header = 0;
//			  this->total_intensity = this->top_intensity = this->single_peak_top_intensity = 0;
//			  this->glycan_cluster_info.clear();
//			  this->glycan_cluster_screening_flag = false;
//			  this->setGlycanClusteringScreeningHistory(0);
//			  this->glycan_cluster_screening_dependency.clear();
//			  this->glycan_cluster_screening_is_approved_by_ion = 0;
//
//			  id = 0;
//			  charge = 0;
//			  start_time = end_time = 0;
//			  total_intensity = 0;
//			  first_mono = 0;
//			  start_scan = end_scan = 0;
//			  size = 0;
//			  last_mono_index = 0;
//			  ms2_score = 0;
//			  priority = priority_2nd = priority_3rd = 0;
//
//			  for (int i=0; i < MAX_ENV_SIZE; i++) 
//				envelope[i].setEmpty();
//
//			  this->siblings.clear();
//			  this->ion_merging_history.clear();
//		}
//
//
//		Isoform::Isoform() 
//		{
//		  init();
//		  id = 0;
//		  charge = 0;
//		  start_time = end_time = 0;
//		  total_intensity = 0;
//		  first_mono = 0;
//		  start_scan = end_scan = 0;
//		  size = 0;
//		  last_mono_index = 0;
//		  this->priority = this->priority_2nd = this->priority_3rd = 0;
//
//		  for (int i=0; i < MAX_ENV_SIZE; i++) 
//			envelope[i].setEmpty();
//		}
//
//		Isoform::Isoform(int id, int charge, float first_mono, float start_time, float end_time, int start_scan, int end_scan) 
//		{
//		  init();
//		  this->id = id;
//		  this->start_time = start_time;
//		  this->end_time = end_time;
//		  this->charge = charge;
//		  this->first_mono = first_mono;
//		  this->total_intensity = 0;
//		  this->start_scan = start_scan;
//		  this->end_scan = end_scan;
//		  this->size = 0;
//		  this->last_mono_index = 0;
//		  this->priority = this->priority_2nd = this->priority_3rd = 0;
//
//		  for (int i=0; i < MAX_ENV_SIZE; i++) 
//			envelope[i].setEmpty();
//		}
//
//		//construct an Isoform from a list of peaks. the
//		//The first mono must be specified. this constructor
//		//does not check for correctness of the input.
//		Isoform::Isoform(int id, int charge, float first_mono, float start_time, float end_time, int start_scan, int end_scan, vector<Peak>& plist)
//		{
//			  //cout << " point 1 " << endl;
//			  init();
//			  this->id = id;
//			  this->start_time = start_time;
//			  this->end_time = end_time;
//			  this->charge = charge;
//			  this->first_mono = first_mono;
//			  this->total_intensity = 0;
//			  this->start_scan = start_scan;
//			  this->end_scan = end_scan;
//			  this->size = 0;
//			  this->last_mono_index = 0;
//			  this->priority = 0;
//
//			  for (int i=0; i < MAX_ENV_SIZE; i++) 
//				envelope[i].setEmpty();
//
//			  assert( charge != 0 );
//			  for (int i=0; i< (int) plist.size(); i++) 
//			  {
//				  int tmp = (int) roundf((plist[i].mass - first_mono)*charge);
//				  //cout << "when ion is " << plist[i].mass << " tmp is " << tmp << endl;
//				  //cout << " pt 1 = " << (plist[i].mass - first_mono)*charge << endl;
//				  //cout << " pt 1 = " << roundf((plist[i].mass - first_mono)*charge) << endl;
//
//					//make a copy of plist[i]
//					//duplicated peaks in plist are ignored.
//					if ( envelope[i].isEmpty() && tmp < MAX_ENV_SIZE )
//					{
//					  envelope[tmp] = plist[i];
//						this->size++;
//						if ( tmp > this->last_mono_index )
//							this->last_mono_index = tmp;
//					}
//				}
//
//			  //cout << " point 2 " << endl;  
//			  //some peaks maybe missing from the pklist (depending on the missing tolerance)
//			  //add an peak with 0 intensity to indicate the missing peak.
//			  //NOTE: the intensity must be 0 or a very small value, because intensity will be
//			  //used as a weight to adjust first-mono later.
//			  bool found_last_mono = false;
//			  for (int i= MAX_ENV_SIZE-1; i >= 0; i--) 
//			  {
//				//find the last non-empty mono peak and mark it
//				//cout << " i " << i << endl; 
//				if ( ! envelope[i].isEmpty() )
//				  found_last_mono = true;
//
//				 // cout << " point 3 " << endl; 
//				if ( found_last_mono )
//				{
//					//cout << " mono peak " << envelope[i].mass << " at charge " << charge << endl; 
//					  if ( envelope[i].isEmpty() ) 
//						envelope[i].mass = first_mono + ((float) i)/charge;
//				}
//			  }
//
//			  //cout << " point 4 " << endl;  
//		}
//
//		Isoform::~Isoform()
//		{
//			  siblings.clear();
//			  attached_sub_ms.clear();
//			  ion_merging_history.clear();
//		}
//
//
//		bool Isoform::isCompleteIsoform() 
//		{
//			return this->size >= 4;
//		}
//
//
//		//this function returns the file IDs of 
//		//experiments in which the ion (isoform)
//		//is scanned for MS/MS. for example, if the ion
//		//scanned for MS/MS in both experiment 1 and
//		//experiment 2. it will return a list of 
//		//integers 1 and 2.
//		vector<int> Isoform::getMsMsPerExperimentHistory() 
//		{
//			map<int, int> fileID_hash;
//			map<int, int>::iterator iter;
//			vector<int> result;
//
//			for (int i=0; i< (int) this->attached_sub_ms.size(); i++) {
//				Spectrum *s = this->attached_sub_ms[i];
//
//				assert( s != NULL && s->file_id >= 0 );
//
//				if ( s != NULL && s->file_id >= 0 )
//					fileID_hash[s->file_id] = s->scanheader_num;
//			}
//
//			result.clear();
//			for (iter = fileID_hash.begin(); iter != fileID_hash.end(); 
//				iter++ ) {
//					result.push_back(iter->first);
//			}
//
//			return result;
//		}
//
//		//recalculate the mz of first mono and
//		//total intensity.
//		//according to the peaks in the envelope 
//		//(take their averages). it does not check
//		//for any possible error or mis-match.
//		void Isoform::update(float first_mono_to_highest_mono_ratio_bound) 
//		{
//		  float fm = 0; //mz value of first mono
//		  float inten = 0; //total intensity
//		  float single_peak_top_inten = 0;
//		  int idx;
//
//		  assert( charge > 0 );
//		  float gap = 1.0/charge;
//
//		  //find the highest mono
//		  inten = 0;
//		  for (int i=0; i< MAX_ENV_SIZE ; i++ ) {
//			if ( envelope[i].intensity > inten )
//			  inten = envelope[i].intensity;
//		  }
//
//		  if ( inten <= 0 ) {
//			init();
//			return;
//		  }
//
//		  //determine first mono
//		  idx = -1;
//		  for (int i=0; i< MAX_ENV_SIZE ; i++ ) {
//			//find the first non-empty peak
//			if ( ! envelope[i].isEmpty() ) {
//
//			  //any of the three conditions
//			  //1. current peak is the last peak OR
//			  //2. next peak is intensity zero OR
//			  //3. current/highest ratio is greater than some ratio
//			  if ( i == MAX_ENV_SIZE -1 ||
//			   envelope[i+1].intensity <= 0 ||
//			   envelope[i].intensity/inten >= 
//			   first_mono_to_highest_mono_ratio_bound ) {
//			idx = i;
//			break;
//			  }
//			}
//		  }
//
//		  if ( idx < 0 ) { //first mono not found
//			init();
//			return;
//		  }
//		  else {
//			//shift the whole envelope.
//			if ( idx != 0 ) {
//			  for (int i=0; i< MAX_ENV_SIZE ; i++ ) {
//			if ( idx + i < MAX_ENV_SIZE )
//			  envelope[i] = envelope[idx+i];
//			else
//			  envelope[i].setEmpty();
//			  }
//			}
//		  }
//
//		  //update information
//		  single_peak_top_inten = 0;
//		  inten = 0;
//		  this->size = 0;
//		  this->last_mono_index = 0;
//		  for (int i=0; i< MAX_ENV_SIZE ; i++ ) {
//			if ( (! envelope[i].isEmpty()) && envelope[i].intensity != 0) {
//			  fm = (fm * inten + (envelope[i].mass - i*gap)*envelope[i].intensity)/
//			(inten + envelope[i].intensity);
//			  inten += envelope[i].intensity;
//			  this->size++;
//			  this->last_mono_index = i;
//
//			  if ( envelope[i].intensity > 
//				  single_peak_top_inten )
//				  single_peak_top_inten = 
//				  envelope[i].intensity;
//			}
//		  }
//
//		  first_mono = fm;
//		  total_intensity = inten;
//		  top_intensity = inten;
//		}
//
//
//		//print the isoform
//		void Isoform::print() 
//		{
//			cout << "Isoform " << id << ": " << start_time 
//					<< " ~ " << end_time << " mins" << endl;
//			cout << "scans =" << start_scan << " ~ " 
//					<< end_scan << endl;
//			cout << "1st mono =" << first_mono << endl;
//			 cout << "charge =" << charge << endl;
//
//			for (int i=0; i<MAX_ENV_SIZE; i++) 
//			{
//				if ( ! envelope[i].isEmpty() )
//					  cout << i << ": " << envelope[i].mass << ", " << envelope[i].intensity << endl;
//			}
//		}
//
//		void Isoform::addAnnotation(string &annotation) 
//		{
//			for (int i=0; i< (int) annotations.size(); i++) 
//			{
//				string &ann = annotations[i];
//
//				if ( ann.compare(annotation) == 0 ) 
//				{
//				  //the annotation already exisits. return without
//				  //adding duplcation
//				  return;
//				}
//			}
//
//			this->annotations.push_back(annotation);
//		}
//
//
//		//deIsotop an MS spectrum. this spectrum 
//		//must be high resuolution scan. the result
//		//is saved in isolist. the matching err is
//		//defined as (-err, +err). the isoform stops
//		//when there are two tandem missing peaks.
//		void Isoform::deIsotop(vector<Peak>& plist, float start_time,  float end_time, int start_scan, int end_scan, int ppm, int max_missing, vector<Isoform>& isolist)
//		{
//			  vector<Peak> peaks, envelope2;
//			  vector<Peak>::iterator iter2, iter3;
//			  float err = 0;
//			 
//			  assert( max_missing > 0 );
//
//			  float first_mono, inten;
//			  int iso_id = 1;
//			  int charge_state = 0;
//			  Isoform iso;
//
//			  peaks.assign(plist.begin(), plist.end());
//
//			  //init the flag of all peaks.
//			  //in this function. the flag will be used to store the isoform ID
//			  for (iter2 = peaks.begin(); 
//				   iter2 != peaks.end(); iter2++) {
//				(*iter2).flag = 0;
//				(*iter2).parent_isoform_id = 0;
//			  }
//
//			  isolist.clear();
//			  for (iter2 = peaks.begin(); 
//				   iter2 != peaks.end(); iter2++) {
//
//				//ignore the peaks whose charge can not be assigned.
//				//ignore the peaks who has been assigned to some envelope
//				//ignore the peaks whose intensity is zero (which is abnormal)
//				if ( (*iter2).charge == 0 || (*iter2).flag > 0 
//				 || (*iter2).intensity <= 0) 
//				  continue;
//			    
//				//start a new isoform
//				iso_id += 1;
//				first_mono = (*iter2).mass;
//				err = PPM2Absolute(ppm, first_mono);
//				inten = 0;
//				charge_state = (*iter2).charge;
//
//				//search for the subsequence mono peaks belonging to the 
//				//current isoform
//				int mono_count = 1;
//				int tandem_missing_peak_count = 0;
//				float gap = 1.0/charge_state;
//				envelope2.clear();
//
//				for (iter3 = iter2; iter3 != peaks.end(); iter3++) {
//				  //check for missing peak
//				  if ( (*iter3).mass >  first_mono + (mono_count-1)*gap + err ) {
//				tandem_missing_peak_count++;
//				if ( tandem_missing_peak_count > max_missing )
//				  break;
//				else {
//				  mono_count += 1;
//				  if ( iter3 != iter2 )
//					iter3--; //roll back the iterator by one 
//				  continue;
//				}  
//				  }
//				  else if ( (*iter3).flag == 0 && (*iter3).charge == charge_state && 
//					(*iter3).mass >= first_mono + (mono_count - 1)*gap - err &&
//					(*iter3).mass <= first_mono + (mono_count - 1)*gap + err) {
//				//dynamically adjust the mz of first mono
//				first_mono = (first_mono * inten + 
//						  ((*iter3).mass - (mono_count-1)*gap) 
//						  * (*iter3).intensity)/
//				  (inten + (*iter3).intensity);
//				inten += (*iter3).intensity;
//
//				//cout << " ion " << (*iter3).mass << " is added" << endl;
//
//				envelope2.push_back((*iter3));
//				(*iter3).parent_isoform_id = iso_id;
//				(*iter3).flag = (float) iso_id; //mark the peaks as "having been visited"
//
//				mono_count += 1; //next mono
//				tandem_missing_peak_count = 0;	
//				  }
//				}
//
//				//cout << "adding ion " << first_mono << endl;
//
//				iso = Isoform(iso_id, charge_state, first_mono,
//					  start_time, end_time, start_scan, 
//					  end_scan, envelope2);    
//				iso.setSinglePeakTopIntensity();
//
//				//cout << first_mono << " is added " << endl;
//
//				iso.update(FIRST_MONO_TO_HIGHEST_MONO_RATIO_BOUND);
//
//				/*
//				if ( start_scan >= 2773 && start_scan <= 2800 ) {
//				  cout << "here: " << iso.first_mono << ", " << iso.charge << endl;
//				}
//				*/
//
//				isolist.push_back(iso);
//			  }
//
//			peaks.clear();
//		}
//
//
//		void Isoform::setSinglePeakTopIntensity() 
//		{
//			this->single_peak_top_intensity = 0;
//
//			for (int i=0; i< MAX_ENV_SIZE; i++) 
//			{
//				if ( this->envelope[i].isEmpty() )
//					continue;
//
//				if ( this->single_peak_top_intensity < this->envelope[i].intensity )
//					this->single_peak_top_intensity = this->envelope[i].intensity;
//			}
//		}
//
//
//
//		//check if current isoform overlaps with iso
//		//two isoform overlaps only if at least 3 peaks
//		//of their envolops overlap, and their charge
//		//states are the same.
//		/*
//		bool Isoform::overlap(Isoform& iso, float err) {
//		  float dif, small_mono;
//		  float tmp1, tmp2;
//		  int overlap_start_index, overlap_end_index;
//
//		  //cout << "pppp 1 " << endl;
//
//		  if ( this->charge != iso.charge )
//			return false;
//
//		  //cout << "pppp 2 " << endl;
//
//		  if ( this->first_mono > iso.first_mono ) {
//			small_mono = iso.first_mono;
//			dif = this->first_mono - small_mono;
//			tmp1 = roundf(dif*charge);
//			overlap_start_index = ((int) tmp1);
//			overlap_end_index = iso.getLastMonoIndex();
//		  }
//		  else {
//			small_mono = this->first_mono;
//			dif = iso.first_mono - small_mono;
//			tmp1 = roundf(dif*charge);
//			overlap_start_index = ((int) tmp1);
//			overlap_end_index = this->getLastMonoIndex();
//		  }
//
//		  //cout << "lap starts " << overlap_start_index 
//		  //   << " lap ends " << overlap_end_index << endl;
//
//		  //two isoforms can not differ by 
//		  //more than two mono peaks
//		  if ( ((int) tmp1) > 3  )
//			return false;
//
//		  //cout << "pppp 3 " << endl;
//
//		  tmp2 = tmp1/charge;
//		 
//		  //if the two isoforms can not be caliburated
//		  if ( dif > tmp2 + err ||
//			   dif < tmp2 - err )
//			return false;
//
//		  //cout << "pppp 4 " << overlap_end_index - overlap_start_index 
//		  //   << endl;
//
//		  //at least 3 peaks overlaping.
//		  int num_overlapping_peaks = overlap_end_index - 
//			overlap_start_index + 1;
//		  if ( num_overlapping_peaks < 2 )
//			return false;
//
//		  //cout << "pppp 5 " << endl;
//		  return true;
//		}
//		*/
//
//		bool Isoform::overlap(Isoform &ion, int ppm_err, float elution_time_shift) 
//		{
//			assert( ppm_err > 0 && elution_time_shift > 0 );
//
//			//check elution time overlap
//			if ( this->start_time >
//				ion.end_time + elution_time_shift ||
//				this->end_time <
//				ion.start_time - elution_time_shift )
//				return false;
//			else 
//				return this->overlap(ion, ppm_err);
//		}
//
//		bool Isoform::overlap(Isoform& iso, int ppm_err) 
//		{
//		  float absolute_err = PPM2Absolute(ppm_err,
//							this->first_mono);
//		  assert( absolute_err >= 0 );
//		  return overlap(iso, absolute_err);
//		}
//
//		bool Isoform::overlap(Isoform& iso, float err)
//		{
//			  float dif, small_mono;
//			  float tmp1, tmp2;
//			  int overlap_start_index, overlap_end_index;
//			  bool is_overlaped = true;
//			  int non_overlap_reason_code = 0;
//
//			  if ( this->charge != iso.charge ) {
//				is_overlaped = false;
//				non_overlap_reason_code = 1;
//			  }
//
//			  if ( is_overlaped ) {
//				if ( this->first_mono > iso.first_mono ) {
//				  small_mono = iso.first_mono;
//				  dif = this->first_mono - small_mono;
//				tmp1 = roundf(dif*charge);
//				  overlap_start_index = ((int) tmp1) - 1;
//				  overlap_end_index = iso.getLastMonoIndex();
//				}
//				else {
//				  small_mono = this->first_mono;
//				  dif = iso.first_mono - small_mono;
//				 // tmp1 = (float) ((int) (dif*charge));
//				tmp1 = roundf(dif*charge);
//				  overlap_start_index = ((int) tmp1) - 1;
//				  overlap_end_index = this->getLastMonoIndex();
//				}
//			  
//				//two isoforms can not differ by 
//				//more than two mono peaks
//				//  if ( ((int) tmp1) > 2  ) {
//				if ( ((int) tmp1) > 10  ) {
//				  /*
//				if ( mim1 > 1548 && mim1 < 1549 &&
//				mim2 > 1548 && mim2 < 1549 ) {
//				cout << "stop1 here " << endl;
//				iso.print(); this->print();
//				}*/
//				  is_overlaped = false;
//				  non_overlap_reason_code = 2;
//				}
//			    
//				tmp2 = tmp1/charge;
//			    
//				//if the two isoforms can not be caliburated
//				if ( is_overlaped ) { 
//				  if ( dif > tmp2 + err ||
//				   dif < tmp2 - err ) {
//				is_overlaped = false;
//				non_overlap_reason_code = 3;
//				  }
//				}
//			    
//				if ( is_overlaped ) {
//				  //at least 3 peaks overlaping.
//				  int num_overlapping_peaks = overlap_end_index - 
//				overlap_start_index;
//				  if ( num_overlapping_peaks < 3 ) {
//				/*
//				  if ( mim1 > 1548 && mim1 < 1549 &&
//				  mim2 > 1548 && mim2 < 1549 ) {
//				  cout << "stop3 here " << endl;
//				  iso.print(); this->print();
//				  }*/
//				is_overlaped = false;
//				non_overlap_reason_code = 4;
//				  }
//				}
//			}
//
//		#ifdef DEBUG_ION_MERGE
//		  if ( ( this->id == DEBUG_ION_MERGE_ION_ONE &&
//			 iso.id == DEBUG_ION_MERGE_ION_TWO ) ||
//			  ( this->id == DEBUG_ION_MERGE_ION_TWO &&
//			 iso.id == DEBUG_ION_MERGE_ION_ONE ) ) {
//			  iso.print(); this->print();
//			  if ( is_overlaped ) 
//			cout << "overlap test is true " << endl;
//			  else
//			cout << "overlap test is false with code " 
//				 << non_overlap_reason_code << endl;
//			}
//		#endif
//
//		return is_overlaped;
//	}
//
//
//		//combine() and finalize() are used together
//		//to combine sibling Isoforms from multiple scans.
//		void Isoform::combine(Isoform& iso) 
//		{
//				siblings.push_back(iso);
//		}	
//
//
//		//combine() and finalize() are used together.
//		//after a set of sibling isoforms are combined.
//		//finalize() extract the most complete isoform from them.
//		//finalize() uses the first-mono of the most intense envelope
//		//as the final first mono.
//		void Isoform::finalize() 
//		{
//		  float min_first_mono;
//		  float intensity_sum[MAX_ENV_SIZE];
//		  int idx, min_start_scan, max_end_scan;
//		  float min_start_time, max_end_time;
//
//		  if ( siblings.size() <= 0 )
//			return; //there is no sibling
//
//		  //1. find the most intense isoform among siblings
//		  //2. find the smallest first mono
//		  int most_intense_iso_index = 0;
//		  min_first_mono = this->first_mono;
//		  min_start_time = this->start_time;
//		  max_end_time = this->end_time;
//		  min_start_scan = this->start_scan;
//		  max_end_scan = this->end_scan;
//
//		  for (int i=0; i < (int) siblings.size(); i++) {
//			if ( siblings[most_intense_iso_index].total_intensity <
//			 siblings[i].total_intensity )
//			  most_intense_iso_index = i;
//
//			if ( min_first_mono > siblings[i].first_mono )
//			  min_first_mono = siblings[i].first_mono;
//		    
//			if ( min_start_time > siblings[i].start_time ) {
//			  min_start_time = siblings[i].start_time;
//			  min_start_scan = siblings[i].start_scan;
//			}
//
//			if ( max_end_time < siblings[i].end_time ) {
//			  max_end_time = siblings[i].end_time;
//			  max_end_scan = siblings[i].end_scan;
//			}
//		  }
//
//		  //align the envelopes by their mono peaks. and
//		  //sum up the intensity at each mono peak
//		  memset(intensity_sum, 0, sizeof(float)*MAX_ENV_SIZE);
//		  siblings.insert(siblings.end(), (*this)); //include self
//		  for (int i=0; i < (int) siblings.size(); i++) 
//		  {
//			Isoform &iso = siblings[i];
//			for (int j=0; j < MAX_ENV_SIZE; j++) 
//			{
//			  if ( ! iso.envelope[j].isEmpty()) 
//			  {
//				idx = (int) roundf(((iso.envelope[j].mass - min_first_mono)*charge));
//				if ( idx >= 0 && idx < MAX_ENV_SIZE ) 
//				{
//				  intensity_sum[idx] += iso.envelope[j].intensity;
//				}
//			  }
//			}
//		  }
//
//		  //copy the first mono into current isoform, if 
//		  //some of the siblings is more intense than the 
//		  //current isoform
//		  Isoform &most_intense_iso = siblings[most_intense_iso_index];
//		  int idx_shift_of_1st_mono = (int) 
//				roundf(((most_intense_iso.first_mono - min_first_mono)*charge));
//
//		  if ( this->total_intensity < 
//			   most_intense_iso.total_intensity) {
//			//update the first mono and each mono mass
//			this->first_mono = most_intense_iso.first_mono;
//			this->last_mono_index = most_intense_iso.last_mono_index;
//			this->top_intensity = most_intense_iso.top_intensity;
//			this->single_peak_top_intensity = most_intense_iso.single_peak_top_intensity;
//
//			for (int j=0; j < MAX_ENV_SIZE; j++) {
//			  if ( (j+idx_shift_of_1st_mono < MAX_ENV_SIZE) ) {
//			if ( most_intense_iso.envelope[j].mass <= 0 ) //some mono peak may be missing
//			  this->envelope[j].mass = this->first_mono + j/((float)this->charge);
//			else this->envelope[j].mass = most_intense_iso.envelope[j].mass;
//			this->envelope[j].intensity = 
//			  intensity_sum[j+idx_shift_of_1st_mono];
//			  }
//			  else {
//			//does not allow missing peak.
//			this->last_mono_index = j;
//			this->envelope[j].setEmpty();
//			break;
//			  }
//			}
//		  }
//
//		  this->first_mono = min_first_mono;
//		  this->start_time = min_start_time;
//		  this->end_time = max_end_time;
//		  this->start_scan = min_start_scan;
//		  this->end_scan = max_end_scan;
//
//		  //combined the screening ion dependency list
//		  if ( glycan_cluster_screening == TURNED_ON ) 
//		  {
//			  for (int i=0; i < (int) siblings.size(); i++) 
//			  {
//				  Isoform &iso = siblings[i]; 
//				  this->addScreeningDependency(
//				  iso.getGlycanClusterScreeningDependency());
//			  }
//		  }
//
//
//		  int test_idx = (int) 
//			roundf(((this->envelope[0].mass - this->first_mono)* this->charge));
//
//		  this->update(0); //set the ratio_bound to 0 to 
//						   //avoid eliminating the first mono
//		  siblings.clear();
//
//		  test_idx = (int) 
//			roundf(((this->envelope[0].mass - this->first_mono)* this->charge));
//
//		  assert ( test_idx >= 0 );
//		}
//
//
//		/*
//		//this commented version of finalize() uses the sophisticated
//		//approach. which chooses the smallest mono as the first mono
//		//of a set of Isoforms from multiple scans. and it takes the avg
//		//mz as the mz of the first mono.
//		void Isoform::finalize() {
//		  float min_first_mono;
//		  int mono_count[MAX_ENV_SIZE];
//		  float intensity_sum[MAX_ENV_SIZE];
//		  float avg_mz[MAX_ENV_SIZE];
//		  Isoform iso;
//		  int idx;
//
//		  if ( siblings.size() <= 0 )
//			return; //there is no sibling
//
//		  min_first_mono = this->first_mono;
//
//		  //find the smallest first mono
//		  for (int i=0; i < (int) siblings.size(); i++) {
//			if ( min_first_mono > siblings[i].first_mono )
//			  min_first_mono = siblings[i].first_mono;
//		  }
//
//		  //count the presence of each mono peak
//		  memset(mono_count, 0, sizeof(int)*MAX_ENV_SIZE);
//		  memset(intensity_sum, 0, sizeof(float)*MAX_ENV_SIZE);
//		  memset(avg_mz, 0, sizeof(float)*MAX_ENV_SIZE);
//		  siblings.insert(siblings.begin(), (*this)); //include self
//		  for (int i=0; i < (int) siblings.size(); i++) {
//			iso = siblings[i];
//			for (int j=0; j < MAX_ENV_SIZE; j++) {
//			  if ( ! iso.envelope[j].isEmpty()) {
//			idx = (int) 
//			  roundf((iso.envelope[j].mass - min_first_mono)*charge);
//			if ( idx >= 0 && idx < MAX_ENV_SIZE ) {
//			  mono_count[idx]++;
//			  if ( iso.envelope[j].intensity > 0 ) {
//				//weighted average of the mz value
//				avg_mz[idx] = (avg_mz[idx]*intensity_sum[idx] + 
//					   iso.envelope[j].mass*iso.envelope[j].intensity)/
//				  (iso.envelope[j].intensity + intensity_sum[idx]);
//			  }
//			  intensity_sum[idx] += iso.envelope[j].intensity;
//			}
//			  }
//			}
//		  }
//		    
//		  //set the first peaks which appears more than 2
//		  //times, as the first mono.
//		  int k = 0;
//		  assert (charge>0 && siblings.size()>0);
//		  for (k=0; k < MAX_ENV_SIZE; k++) {
//			if ( mono_count[k] >= 2 ) {
//			  this->first_mono = avg_mz[k];
//			  break;
//			}
//		  }
//
//		  //update the avg intensity of each mono peak
//		  for (int j=0; j < MAX_ENV_SIZE; j++) {
//			if ( j+k < MAX_ENV_SIZE && 
//			 intensity_sum[j+k] > 0 ) {
//			  this->envelope[j].mass = avg_mz[j+k];
//			  //this->envelope[j].intensity = intensity_sum[j+k]/num_siblings;
//			  this->envelope[j].intensity = intensity_sum[j+k];
//			}
//			else
//			  this->envelope[j].setEmpty();
//		  }
//
//		  this->update(0); //set the ratio_bound to 0 to 
//						   //avoid eliminating the first mono
//		  siblings.clear();
//		}
//		*/
//
//
//		//test if the isoform is covered by the ion acquiztion
//		//window of the MS2. An isoform is covered by the ion
//		//acquizition window if 10% of its intensity is covered
//		//by this window.
//		/*
//		bool Isoform::isCoveredBySubMS(Spectrum* ms2, 
//						   float mz_tolerance) {
//		  float precursor = ms2->precursor_mass;
//		  float total_intensity, covered_intensity;
//		  float min_cover_ratio = 0.1;
//
//		  assert( mz_tolerance > 0 );
//
//		  if ( ms2->scanheader_num < this->start_scan ||
//			   ms2->scanheader_num > this->end_scan + 10 )
//			return false;
//
//		  total_intensity = covered_intensity = 0;
//		  for (int i=0; i< MAX_ENV_SIZE; i++) {
//			if ( envelope[i].isEmpty() ) 
//			  continue;
//			else {
//			  total_intensity += envelope[i].intensity;
//		      
//			  if ( envelope[i].mass >= precursor - mz_tolerance  &&
//			   envelope[i].mass <= precursor + mz_tolerance  ) {
//			covered_intensity += envelope[i].intensity;
//			  }
//			}
//		  }
//
//		  if ( total_intensity == 0 || 
//			   covered_intensity/total_intensity <
//			   min_cover_ratio ) 
//			return false;
//		  else 
//			return true;
//		}
//		*/
//
//
//		//test if the isoform is covered by the ion acquiztion
//		//window of the MS2. An isoform is covered by the ion
//		//acquizition window if its most intensed ion is covered
//		//by this window.
//		bool Isoform::isCoveredBySubMS(Spectrum* ms2, float mz_tolerance) 
//		{
//		  float precursor = ms2->precursor_mass;
//		  assert( mz_tolerance > 0 );
//
//		  if ( ms2->scanheader_num < this->start_scan ||
//			   ms2->scanheader_num > this->end_scan + 10 )
//			return false;
//
//		  float most_intense_ion = this->getMostIntenseMono();
//
//		  if ( most_intense_ion >= precursor - mz_tolerance  &&
//			   most_intense_ion <= precursor + mz_tolerance  )
//			return true;
//		  else
//			return false;
//		}
//
//
//
//		bool Isoform::acceptSubMS(Spectrum* ms2) 
//		{
//		  float allowed_err_ppm = 30;
//		  float allowed_err_da = envelope[0].mass * allowed_err_ppm / 1000000;
//		  float effective_time = 0.1; //minute
//
//		  if ( (ms2->start_time < start_time - effective_time) ||
//			   (ms2->start_time > end_time + effective_time) )
//			return false;
//
//		  return isPartOf(ms2->precursor_mass, allowed_err_da);;
//		}
//
//
//		//test if an peak (m/z) is part of the current
//		//isoform. given (-err, +err) window
//		bool Isoform::isPartOf(float mz, float err) {
//		  bool result = false;
//
//		  for (int i=0; i< MAX_ENV_SIZE; i++) {
//			if ( envelope[i].isEmpty() ) 
//			  continue;
//			else if ( envelope[i].mass - err > mz  )
//			  break;
//			else if ( envelope[i].mass + err >= mz  ) {
//			  result = true;
//			  break;
//			}
//		  }
//
//		  return result;
//		}
//
//
//		//returns the m/z of the last mono peak.
//		float Isoform::getLastMono() {
//		  float result = 0;
//
//		  for (int i=0; i< MAX_ENV_SIZE; i++) {
//			//if ( ! envelope[i].isEmpty() ) 
//			//result = envelope[i].mass;
//			if ( envelope[i].intensity > 0 ) 
//			  result = envelope[i].mass;
//			else
//			  break;
//		  }
//
//		  return result;
//		}
//
//		//returns the index (starting from 1) 
//		//of the last mono peak.
//		int Isoform::getLastMonoIndex() {
//		  int result = -1;
//
//		  for (int i=0; i< MAX_ENV_SIZE; i++) {
//			if ( ! envelope[i].isEmpty() ) 
//			  result = i;
//		  }
//
//		  return result;
//		}
//
//		  
//		//returns the m/z of the most intense peak.
//		float Isoform::getMostIntenseMono() {
//		  float result = 0;
//		  float highest_intensity = -1;
//
//		  //remove this line!!!???
//		  float last_mono = this->getLastMono(); 
//
//		  for (int i=0; i< MAX_ENV_SIZE; i++) {
//			if ( envelope[i].mass > last_mono )
//			  break;
//
//			if ( ! envelope[i].isEmpty()  ) {
//			  if ( envelope[i].intensity > highest_intensity ) {
//			result = envelope[i].mass;
//			highest_intensity = envelope[i].intensity;
//			  }
//			}
//		  }
//
//		  return result;
//		}
//
//		float Isoform::getIntensityOfMostIntenseMono() {
//		  float highest_intensity = -1;
//
//		  //remove this line!!!???
//		  float last_mono = this->getLastMono(); 
//
//		  for (int i=0; i< MAX_ENV_SIZE; i++) {
//			if ( envelope[i].mass > last_mono )
//			  break;
//
//			if ( ! envelope[i].isEmpty()  ) {
//			  if ( envelope[i].intensity > highest_intensity ) {
//			highest_intensity = envelope[i].intensity;
//			  }
//			}
//		  }
//
//		  return highest_intensity;
//		}
//
//		//convert the current isoform to the peak whose
//		// m/z is the first mono.
//		Peak Isoform::toFirstMonoPeak() {
//		  Peak p = Peak(first_mono, total_intensity, charge);
//		  p.parent_isoform_id = this->id;
//
//		  //attach MS2
//		  for (int i=0; 
//			   i< MAX_NUM_SUB_MS && i< (int) attached_sub_ms.size(); i++) {
//			p.sub_ms[i] = attached_sub_ms[i];
//		  }
//		  
//		  return p;
//		}
//
//		//convert the current isoform to the peak whose
//		// m/z is the, possibly, missed first mono.
//		Peak Isoform::toMissedFirstMonoPeak() {
//		  Peak p = Peak(first_mono - 1.0/charge, total_intensity, charge);
//		  p.parent_isoform_id = this->id;
//
//		  //attach MS2
//		  for (int i=0; 
//			   i< MAX_NUM_SUB_MS && i< (int) attached_sub_ms.size(); i++) {
//			p.sub_ms[i] = attached_sub_ms[i];
//		  }
//		  
//		  return p;
//		}
//
//
//		void Isoform::addSubMS(Spectrum *ms) {
//		  assert( ms != NULL );
//		  assert( ms->scanheader_num > 0);
//		  //assert( ms->file_id > 0 );
//
//		  //first check if the ms is already part of
//		  //the attached sub_ms
//		  for (int i=0; i < (int) attached_sub_ms.size(); i++) {
//			if ( attached_sub_ms[i]->scanheader_num == ms->scanheader_num &&
//			 attached_sub_ms[i]->file_id == ms->file_id )
//			  return;
//		  }
//
//		  if ( ms->scanheader_num == 5789 ) 
//		  {
//			  cout << "scan 5789 is attached to iso " << this->id << endl;
//			  cout << "precursor of 5789 is " << ms->precursor_mass << endl;
//			  this->print();
//			  ms->print();
//		  }
//
//		  //insert the ms into the list by the order of its
//		  //ms2 score. (from high score to low score)
//		  if ( attached_sub_ms.size() == 0 )
//			attached_sub_ms.push_back(ms);
//		  else {
//			vector<Spectrum*>::iterator iter;
//			for (iter = attached_sub_ms.begin(); 
//			 iter !=  attached_sub_ms.end(); iter++) {
//			  if ( (*iter)->totalScore() < ms->totalScore() ) {
//			attached_sub_ms.insert(iter, ms);
//			return;
//			  }
//			}
//
//			if ( iter == attached_sub_ms.end() ) {
//			  attached_sub_ms.push_back(ms);
//			  return;
//			}
//			else {
//			  cout << "Unable to insert sub ms in Isoform::addSubMS()." << endl;
//			  assert(0);
//			}
//		  }
//		}
//
//
//		float Isoform::getIncListIon() {
//		  return this->getMostIntenseMono();
//		  //return this->first_mono;
//		}
//
//		//find the best glyco-MS2
//		Spectrum* Isoform::getKeySubMS() {
//		  float ms2_score, best_ms2_score;
//		  Spectrum *key_spec;
//		  best_ms2_score = ms2_score = 0;
//		  key_spec = NULL;
//		  
//		  for (int i=0; i< (int) attached_sub_ms.size();
//			   i++) {
//			assert ( attached_sub_ms[i] != NULL );
//			ms2_score = attached_sub_ms[i]->totalScore();
//			if ( best_ms2_score < ms2_score ) {
//			  best_ms2_score = ms2_score;
//			  key_spec = attached_sub_ms[i];
//			}
//		  }
//		  return key_spec;
//		}
//
//		void Isoform::updateScreeningDependency(vector<float> &screen_ion_id_list) 
//		{
//			if ( screen_ion_id_list.size() <= 0 )
//				return;
//			else {
//				this->glycan_cluster_screening_dependency.insert(
//				this->glycan_cluster_screening_dependency.begin(),
//				screen_ion_id_list.begin(), screen_ion_id_list.end());
//				removeDuplicates(this->glycan_cluster_screening_dependency);
//			}
//		}
//
//
//		void Isoform::addScreeningDependency(vector<float> &screen_ion_id_list) 
//		{
//			if ( screen_ion_id_list.size() <= 0 )
//				return;
//			else {
//				this->glycan_cluster_screening_dependency.insert(
//				this->glycan_cluster_screening_dependency.end(),
//				screen_ion_id_list.begin(), screen_ion_id_list.end());
//				removeDuplicates(this->glycan_cluster_screening_dependency);
//			}
//		}
//
//		vector<float>& Isoform::getGlycanClusterScreeningDependency()
//		{
//			return this->glycan_cluster_screening_dependency;
//		}
//
//		void Isoform::updateFileID(int id) 
//		{
//			assert( id > 0 );
//			this->file_id_in_which_ion_is_found = id;
//
//			if ( this->glycan_cluster_screening_flag ) 
//				this->setGlycanClusteringScreeningHistory(id);
//		}
//
//		bool Isoform::isApprovedInScreening()
//		{
//			if ( this->glycan_cluster_screening_is_approved_by_ion < 0 )
//				return false;
//
//			int ion_id = 
//				getID1FromCombination(this->glycan_cluster_screening_is_approved_by_ion);
//			int data_file_id = 
//				getID2FromCombination(this->glycan_cluster_screening_is_approved_by_ion);
//			return ( ion_id > 0 && data_file_id > 0 );
//		}
//
//		float Isoform::getApprovedInScreeningBy()
//		{
//			return this->glycan_cluster_screening_is_approved_by_ion;
//		}
//
//		void Isoform::setApprovedInScreening(int ion_id,
//											 int data_file_id)
//		{
//			this->glycan_cluster_screening_is_approved_by_ion = 
//				idNumberToCombination(ion_id, data_file_id);
//		}
//
//
//		//check if the current ion is in the screening
//		//ion list for raw of the the data_file_id.
//		//e.g. if ion is put into the screening
//		//ion list, and is subject to MS/MS in the experiment
//		//which generate file i. then
//		// ion.isInScreenIonList(i) should return true.
//		bool Isoform::isInScreenIonList(int file_id)
//		{
//			assert(this->glycan_cluster_screening_history >= 0);
//		    
//			if ( this->glycan_cluster_screening_history > 0 &&
//				this->glycan_cluster_screening_history == file_id )
//				return true;
//			else return false;
//		}
//
//		void Isoform::setInScreenIonList(int data_file_id) 
//		{
//			assert( data_file_id > 0 );
//
//			if ( data_file_id > 0 )
//			{
//				if ( (! this->isIdentifiedAsGlycopeptide() ) &&
//					this->charge > 1 && this->charge < 6 )
//					this->setGlycanClusteringScreeningHistory(data_file_id);
//				else 
//				{
//					cout << "An invalid ion is going to be set as screening ion in setInScreenIonList()." << endl;
//					exit(0);
//				}
//			}
//			else 
//			{
//				cout << "Invalid file id in setInScreenIonList()." << endl;
//				exit(0);
//			}
//		}
//
//		bool Isoform::isIdentifiedAsGlycopeptide() 
//		{
//			Spectrum *spec = this->getKeySubMS();
//
//			float threshold = 6.0;
//			if ( this->ms2_score >= threshold || 
//				( spec != NULL && spec->totalScore() >= threshold ))
//				return true;
//			else return false;
//		}
//
//		void Isoform::clearGlycanClusterScreeningDependency() 
//		{
//			this->glycan_cluster_screening_dependency.clear();
//		}
//
//		int Isoform::getGlycanClusteringScreeningHistory() 
//		{
//			return this->glycan_cluster_screening_history;
//		}
//
//		void Isoform::setGlycanClusteringScreeningHistory(int file_id)
//		{
//			assert(file_id >= 0);
//			this->glycan_cluster_screening_history = file_id;
//		}
//	}
//}

#endif
