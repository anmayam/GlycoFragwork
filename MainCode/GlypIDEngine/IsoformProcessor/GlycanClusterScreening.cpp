// Written by Yin Wu (wuyin@indiana.edu)
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// modified by Anoop for GlypID




#ifndef GLYCAN_CLUSTER_SCREENING_CPP
#define GLYCAN_CLUSTER_SCREENING_CPP

//#include "../Utilities/system.h"

namespace Engine
{
	namespace IsoformProcessor
	{

		////if the screening method is turned on,
		////check identification results from previous
		////screening list. 
		//void checkPendingScreenedIons(Dataset &data, 
		//							  vector<Isoform> &old_ion_list,
		//							  int data_file_id) 
		//{
		//	//get the screening ion list from the old_ion_list
		//	vector<Isoform*> screening_ion_list;
		//	for (int i=0; i< (int) old_ion_list.size(); i++) 
		//	{
		//		if ( old_ion_list[i].isInScreenIonList(data_file_id) )
		//			screening_ion_list.push_back(&(old_ion_list[i]));
		//	}

		//	//for each screening ion, check if it is identified in the
		//	//new file or not. produce an identidifed list.
		//	vector<Isoform> &new_ion_list = data.extracted_isoforms;
		//	vector<Isoform*> identified_screening_ion_list;
		//	for (int i=0; i < (int) screening_ion_list.size(); i++) 
		//	{
		//		Isoform *screening_ion = screening_ion_list[i];
		//		for (int j=0; j< (int) new_ion_list.size(); j++) 
		//		{
		//			if ( new_ion_list[j].isIdentifiedAsGlycopeptide() )
		//				if ( screening_ion->overlap(new_ion_list[j], 
		//					FT_PPM_ERROR, peptide_elute_time_shift) ) 
		//					identified_screening_ion_list.push_back(screening_ion);
		//		}
		//	}

		//	//check each pending-to-be-approved ion, if it is approved
		//	//by any of the identified screening ion, then mark it.
		//	for (int i=0; i< (int) old_ion_list.size(); i++) 
		//	{
		//		Isoform &ion = old_ion_list[i];

		//		if (! ion.isApprovedInScreening() ) 
		//		{
		//			vector<float> &glycan_cluster_screening_dependency = 
		//				ion.getGlycanClusterScreeningDependency();
		//			for (int j=0; 
		//				j < (int) glycan_cluster_screening_dependency.size();
		//				j++) 
		//			{
		//				float ion_id_combination = 
		//					glycan_cluster_screening_dependency[j];
		//				int depending_on_screen_ion_id = 
		//					getID1FromCombination(ion_id_combination);
		//                    
		//				for (int k=0; 
		//					k < (int) identified_screening_ion_list.size();
		//					k++) 
		//				{
		//					if ( identified_screening_ion_list[k]->id ==
		//						depending_on_screen_ion_id ) 
		//					{
		//						ion.setApprovedInScreening(depending_on_screen_ion_id, 1); //need to use a valid data file id here.
		//					}
		//				}
		//			}
		//		}

		//		//clear the screening dependency. since we do not carry it
		//		//over to the next experiment.
		//		ion.clearGlycanClusterScreeningDependency();
		//	}
		//}

		////After merging the new and old ion list, some new ions are
		////merged with old ions, and their ion id change. 
		////this function finds the updated ion id and file id for
		////any given original ion_id and file id.
		////if no updated ion is found, return 0.
		//float getUpdatedIonID(vector<Isoform> &merged_ion_list, 
		//					  int ion_id, int data_file_id) 
		//{
		//	assert( ion_id > 0 && data_file_id > 0 );

		//	float ion_id_combination = 
		//		idNumberToCombination(ion_id, data_file_id);

		//	for (int i=0; i< (int) merged_ion_list.size(); i++) 
		//	{
		//		Isoform &ion = merged_ion_list[i];

		//		//check if 'ion' is the one that is being searched
		//		if ( ion.id == ion_id && ion.file_id_in_which_ion_is_found == data_file_id )
		//		{
		//			return idNumberToCombination(ion_id, data_file_id);
		//		}

		//		//check the merge history of 'ion'
		//		for (int j=0; j< (int) ion.ion_merging_history.size(); 
		//			j++)   
		//		{
		//			if ( ion.ion_merging_history[j] == ion_id_combination )
		//			{
		//				return idNumberToCombination(ion.id, ion.file_id_in_which_ion_is_found);   
		//			}
		//		}
		//	}

		//	return 0;
		//}


		////Another version of getUpdatedIon(), but returns the
		////pointer to the ion. 
		////if no updated ion is found, return NULL.
		//Isoform* getUpdatedIon(vector<Isoform> &merged_ion_list, 
		//					  int ion_id, int data_file_id) 
		//{
		//	assert( ion_id > 0 && data_file_id > 0 );

		//	float ion_id_combination = 
		//		idNumberToCombination(ion_id, data_file_id);

		//	for (int i=0; i< (int) merged_ion_list.size(); i++) 
		//	{
		//		Isoform &ion = merged_ion_list[i];

		//		//check if 'ion' is the one that is being searched
		//		if ( ion.id == ion_id && ion.file_id_in_which_ion_is_found == data_file_id )
		//		{
		//			return &ion;
		//		}

		//		//check the merge history of 'ion'
		//		for (int j=0; j< (int) ion.ion_merging_history.size(); 
		//			j++)   
		//		{
		//			if ( ion.ion_merging_history[j] == ion_id_combination )
		//			{
		//				return &ion;  
		//			}
		//		}
		//	}

		//	return NULL;
		//}

		////After merging the new and old ion list, the dependency list of the
		////new ions should be updated, because some screening ions are merged with
		////old ions, and therefore the new ions id is replaced by the
		////old ion's id. Other ions, which depends on this new ion, still
		////have the original id in there dependency list. therefore, the
		////dependency list of these ions must be updated.
		//void updateScreeningDependency(vector<Isoform> &merged_ion_list) 
		//{
		//	vector<float> new_dependency_list;

		//	for (int i=0; i< (int) merged_ion_list.size(); i++) 
		//	{
		//		Isoform &ion = merged_ion_list[i];

		//		if (!  ion.isApprovedInScreening() ) 
		//		{
		//			vector<float> &glycan_cluster_screening_dependency =
		//				ion.getGlycanClusterScreeningDependency();
		//			new_dependency_list.clear();
		//			for (int j=0; j< (int) glycan_cluster_screening_dependency.size(); 
		//				j++)
		//			{
		//				float ion_id_combination = glycan_cluster_screening_dependency[j];
		//				int ion_id = getID1FromCombination(ion_id_combination);
		//				int data_file_id = getID2FromCombination(ion_id_combination);
		//				ion_id_combination = getUpdatedIonID(merged_ion_list, ion_id, data_file_id);
		//				assert(ion_id_combination > 0 );
		//				new_dependency_list.push_back(ion_id_combination);
		//			}
		//		}

		//		ion.clearGlycanClusterScreeningDependency();
		//		ion.updateScreeningDependency(new_dependency_list);
		//	}
		//}


		//void updateScreeningIon(vector<Isoform> &merged_ion_list, 
		//						vector<Isoform> &ion_list_b4_merge,
		//						int data_file_id) 
		//{
		//	vector<float> new_dependency_list;

		//	int count = 0;
		//	for (int i=0; i< (int) ion_list_b4_merge.size(); i++) 
		//	{
		//		Isoform &ion = ion_list_b4_merge[i];

		//		if ( ion.isInScreenIonList(data_file_id) ) 
		//		{
		//			count++;
		//			int ion_id = ion.id;
		//			Isoform *ion_in_merged_list = getUpdatedIon(merged_ion_list, ion_id, data_file_id);
		//			assert( ion_in_merged_list != NULL );
		//			if ( ion_in_merged_list != NULL && 
		//				(! ion_in_merged_list->isIdentifiedAsGlycopeptide()) )
		//			{   
		//				//if ( ion_in_merged_list->id != ion_id )
		//				//    cout << "screen ion " << ion_id << " => "
		//				//    << ion_in_merged_list->id << endl;
		//		if ( ion_in_merged_list->ms2_score >= 7 )
		//			cout << "ion **" << ion_in_merged_list->id 
		//			<< "," << ion_in_merged_list->ms2_score << endl;

		//				ion_in_merged_list->setInScreenIonList(data_file_id);
		//			}   
		//		}
		//	}

		//	// cout << "roundx " << endl;
		//	//for (int i=0; i< (int) merged_ion_list.size(); i++) 
		//	//{
		//	//    Isoform &iso = merged_ion_list[i];

		//	//    if ( iso.priority < inc_list_priority_threshold && 
		//	//        iso.isInScreenIonList(data.getFileID()) )
		//	//        cout << "ion #" << iso.id << "," << iso.priority << endl;

		//	//    if ( iso.ms2_score >= 7 && 
		//	//        iso.isInScreenIonList(data.getFileID()) )
		//	//        cout << "ion ##" << iso.id << "," << iso.ms2_score << endl;
		//	//}

		//	//cout << "there are " << count << " screening ions" << endl;
		//}
	}
}

#endif