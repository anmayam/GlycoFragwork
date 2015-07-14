/*****************************************************************************\
 * Spectrum.cpp
 *
 * defination of Spectrum class
 *
 * author: Yin Wu (wuyin@indiana.edu)
 *
\*****************************************************************************/

#include "./Spectrum.h"
#include "../Utilities/system.h"



namespace Engine
{
	namespace MS2CIDScoring
	{	

		//the accurate masses of the 4 monosaccharides GlcNac, mannose, fucose, NeuNac
		float Spectrum::monosaccharide_masses[Spectrum::num_monosaccharides]=
						{203.0794f, 162.0528f, 146.0579f, 291.0954f};

		float Spectrum::oxonium_ions[Spectrum::num_oxonium_ions] = {366.33, 657.59, 204, 274, 292} ; // Keep this ordering since the first two are the most likely candidates

		Spectrum::Spectrum(void)
		{			
				mdbl_min_mz = 0;
				mdbl_max_mz = 2000 ; 				
				m_sort_mode = Engine::PeakProcessing::SORT_MODE::SORT_BY_INTENSITY ; 
				m_sort_order = Engine::PeakProcessing::SORT_ORDER::DESCENDING ; 
				this->l_table = NULL;
				

		}

		Spectrum::~Spectrum(void)
		{
			if ( l_table != NULL )
				delete(l_table);
			if ( this != NULL ) 
			{
				this->peaks.clear();				
			}
		}

		void Spectrum::sampleMS2Distribution() 
		{	
			//this is a pesudo one
			mobj_ms2_distribution.confidence_interval_list.clear();

			mobj_ms2_distribution.confidence_interval_list.push_back(Engine::PeakProcessing::GPixel(0.833, 1.0, 0, 0));

			mobj_ms2_distribution.confidence_interval_list.push_back(Engine::PeakProcessing::GPixel(0.15, 2.0, 0, 0));

			mobj_ms2_distribution.confidence_interval_list.push_back(Engine::PeakProcessing::GPixel(0.015, 3.0, 0, 0));

			mobj_ms2_distribution.confidence_interval_list.push_back(Engine::PeakProcessing::GPixel(0.0017, 4.0, 0, 0));

		}
		
		//The recurrence function that does the actual
		//enumeration for allMassDif().
		//Note that: The init value of "mass" should be 0.
		//need improvment!!!
		void Spectrum::enumerateMassDif(vector<float> &list, float mass, int mass_idx, int allowed_missing) 
		{
				int i;
				float new_mass;

				for (i=mass_idx; i<num_monosaccharides; i++) 
				{
					new_mass = mass + monosaccharide_masses[i];
					list.push_back(new_mass);
					if ( allowed_missing > 0 ) 
						enumerateMassDif(list, new_mass, i, allowed_missing-1);
				}
		}

		float Spectrum::MS2ScoreToPValue(int score)
		{
			vector<Engine::PeakProcessing::GPixel>::iterator iter1;
			this->sampleMS2Distribution();
			int size = mobj_ms2_distribution.confidence_interval_list.size();

			//improve the performance of this programm!!!
			for (iter1 = mobj_ms2_distribution.confidence_interval_list.begin();
				iter1 != mobj_ms2_distribution.confidence_interval_list.end(); iter1++) 
			{
				if ( score < (*iter1).x )
					return (*iter1).value;
			}

			return 1;
		}

		//gives a list of all possible featured mass differences
		//betwee the glycan chains, depending on how many missing
		//peaks allowed. this list is derived from the
		//pre-defined (or pre-loeaded) possible masses
		void Spectrum::allMassDif(vector<float> &m_list, int allowed_missing)
		{
				vector<float>::iterator iter;

				assert( num_monosaccharides > 0 );
				enumerateMassDif(m_list, 0.0, 0, allowed_missing);
				sort(m_list.begin(), m_list.end());
				iter = unique(m_list.begin(), m_list.end());

				while ( m_list.end() != iter )
					m_list.pop_back();
		}


		bool Spectrum::addPeak(Engine::PeakProcessing::CID_Peak &p) 
		{
			if ( p.mz <= 0 || p.intensity <= 0 )
				return false;
			else 
			{
				this->peaks.push_back(p);
				return true;
			}
		}
	

		void Spectrum::ConvertPeaksToCIDPeaks(Engine::PeakProcessing::PeakData &peak_data)
		{
			Engine::PeakProcessing::CID_Peak this_cid_peak ; 
			Engine::PeakProcessing::Peak this_peak ; 
			/*for (int i =0; i < peak_data.GetNumPeaks(); i++)
			{
				peak_data.GetNextPeak(mdbl_min_mz, mdbl_max_mz, this_peak) ; 
				this_cid_peak.intensity = (float) this_peak.mdbl_intensity ;
				this_cid_peak.mz = (float) this_peak.mdbl_mz ; 
				addPeak(this_cid_peak) ; 				
			}*/
			peaks.clear() ; 
			for (int i = 0 ; i < peak_data.GetNumPeaks() ; i++)
			{
				peak_data.GetPeak(i, this_peak) ; 
				this_cid_peak.intensity = (float) this_peak.mdbl_intensity ;
				this_cid_peak.mz = (float) this_peak.mdbl_mz ; 
				addPeak(this_cid_peak) ; 	
			}

		}

		bool peakGreaterByMz(Engine::PeakProcessing::CID_Peak p1, Engine::PeakProcessing::CID_Peak p2)
		{
			return p1.mz > p2.mz;
		}

		bool peakGreaterByIntensity(Engine::PeakProcessing::CID_Peak p1, Engine::PeakProcessing::CID_Peak p2)
		{
			return p1.intensity > p2.intensity;
		}

		bool peakSmallerByMz(Engine::PeakProcessing::CID_Peak p1, Engine::PeakProcessing::CID_Peak p2) 
		{
			return p1.mz < p2.mz;			
		}

		bool peakSmallerByIntensity(Engine::PeakProcessing::CID_Peak p1, Engine::PeakProcessing::CID_Peak p2) 
		{
			return p1.intensity < p2.intensity;			
		}
		
		void Spectrum::setSortOrder(Engine::PeakProcessing::SORT_ORDER order) 
		{
			this->m_sort_order = order ; 
		}

		void Spectrum::setSortMode(Engine::PeakProcessing::SORT_MODE mode) 
		{
			this->m_sort_mode = mode  ; 
		}

		void Spectrum::sortPeaks(Engine::PeakProcessing::SORT_MODE mode, Engine::PeakProcessing::SORT_ORDER order) 
		{			
			m_sort_mode = mode ; 
			m_sort_order = order ; 
			switch ( order ) 
			{
				case Engine::PeakProcessing::SORT_ORDER::ASCENDING:
					if (mode == Engine::PeakProcessing::SORT_MODE::SORT_BY_INTENSITY)
						sort( peaks.begin(), peaks.end(), peakSmallerByIntensity);
					else
						sort( peaks.begin(), peaks.end(), peakSmallerByMz);
					break;
				case Engine::PeakProcessing::SORT_ORDER::DESCENDING: default:
					if (mode == Engine::PeakProcessing::SORT_MODE::SORT_BY_INTENSITY)
						sort( peaks.begin(), peaks.end(), peakGreaterByIntensity);
					else
						sort( peaks.begin(), peaks.end(), peakGreaterByMz);
					break ; 
						
			}
		}

		bool Spectrum::containsOxoniumion(Engine::PeakProcessing::PeakData &peak_data)
		{
			//---------- Function to look through the peak list and identify oxonium ions ----------// 
			//	Oxonium ions come from the terminal saccharide product (cut-off from ends of the glycan chains).		
			Engine::PeakProcessing::Peak currentPeak ; 

			bool found = false ; 

			for(int i = 0; i < peak_data.GetNumPeaks() ; i++)
			{
				peak_data.GetPeak(i, currentPeak) ; 				
				for (int j= 0; j < num_oxonium_ions; j++)
				{
					double this_mz = oxonium_ions[j] ; 
					if(abs(currentPeak.mdbl_mz-this_mz) < MZ_ERROR)
					{
						found= true ; 
						return found; 						
					}
				}
			}
			/*for (int i = 0 ; i < num_oxonium_ions; i++)
			{
				//peak_data.FindPeak(oxonium_ions[i]- MZ_ERROR, oxonium_ions[i] + MZ_ERROR, ox_peak) ; 		
				if (ox_peak.mdbl_mz >0)
				{
					found = true; 
					break;
				}
			}*/
			
			return found ; 
		}

		//pick the greatest (intensity) N peaks 
		//in the spectrum and filter the others.
		void Spectrum::topNPeaks(int n) 
		{
			int size = peaks.size();

			if ( size > n ) 
			{
				vector<Engine::PeakProcessing::CID_Peak> sorted_peaks;
				vector<Engine::PeakProcessing::CID_Peak>::iterator iter1;
				this->rankPeaksByIntensity();
				for (iter1 = peaks.begin(); iter1 != peaks.end(); iter1++) 
				{
					if ( iter1->intensity_rank <= n )
						sorted_peaks.push_back( (*iter1) );
				}
				peaks.clear(); 
				peaks.assign(sorted_peaks.begin(), sorted_peaks.end());
				sorted_peaks.clear();
			}
		}

		void Spectrum::rankPeaksByIntensity()
		{
			this->sortPeaks(Engine::PeakProcessing::SORT_MODE::SORT_BY_INTENSITY, Engine::PeakProcessing::SORT_ORDER::DESCENDING);
			for (int i=0; i< (int) peaks.size(); i++) 
				peaks[i].intensity_rank = i+1;
			this->sortPeaks(Engine::PeakProcessing::SORT_MODE::SORT_BY_MZ, Engine::PeakProcessing::SORT_ORDER::ASCENDING);		
		}	

		
		//this methods finds all connected-monosaccharide path in the spectrum,
		//and record them. the score of the spectrum should normally be the size
		//of the longest connected-monoscacchride path. but this method does not 
		//return a score immediately. because, the way for computing the score
		//can be slightly different depending on different considerations. e.g.
		//sometimes you want to count the oxonium ion (mass 366) as an additional
		//point of score etc. Therefore , call getScore(scoreVersion) to get the 
		//score you want.
		// getScore deleted until future needs
		float Spectrum::computeMs2Score(vector<float> &mass_dif, Engine::PeakProcessing::PeakData &pk_data, int charge) 
		{
			
			// Convert to Yin's structure
			ConvertPeaksToCIDPeaks(pk_data) ; 
			// take top 20 peaks.
			topNPeaks(20) ; 

			float dif, max_d ; 
			int step, i, j, k, path_length, end_pos;  //max_path_length, 
			int num_peaks = peaks.size();
			int s = (num_peaks+1)*num_peaks/2;
			int *link_tbl, *gap_too_wide;

			link_tbl = new int[s];
			initTable(link_tbl, num_peaks);		

			gap_too_wide = new int[num_peaks];
			memset(gap_too_wide, 0, num_peaks*sizeof(int));
			max_d = (float)((*mass_dif.rbegin()) + MZ_ERROR) ; 	

			for (i=0; i < num_peaks-1; i++) 
			{
				//AM - Haixu suggested using all cs <= parent charge, but commenting this out for Yin
				//for (int cs = charge ; cs > 0 ; cs--) 
				//{
					for (j = i+1; j < num_peaks; j++)	
					{
						
						dif = ((peaks.at(j)).mz - (peaks.at(i)).mz) * charge;
						//assert(dif >= 0);
						
						//the gap is too wide. no need to continue
						if ( dif > max_d ) 
						{
							break;
						}
						else 
						{
							float accurate_dif = approxFindWithCharge(mass_dif, dif, charge, INTRA_SPECTRUM_ION_TRAP_MZ_ERROR);
							if ( accurate_dif > 0 ) 
							{   
								updateLinkTable(link_tbl, num_peaks, i, j);
								double mz1 ; 
								double mz2 ; 
								mz1 = peaks.at(j).mz ; 
								mz2 = peaks.at(i).mz ; 
							}
						}
					}
				//}
			}

		  //use dynamic programming to find the maximum connected tree.
			int *tree_size_table = new int[num_peaks];
			float *branch = new float[mass_dif.size()]; //maximum branches that a node can have
			int max_tree_size = 0;
			int tree_size;
			int start_pos = end_pos = 0;
			memset(tree_size_table, 0, sizeof(int)*num_peaks);
			
			for (i= num_peaks-2; i >= 0; i--) 
			{
				tree_size = 0; //size of the tree starting from node (peak) i
				memset(branch, 0, sizeof(float)*mass_dif.size());
				for (j=i; j<num_peaks; j++) 
				{
					path_length = getLink(link_tbl, num_peaks, j, i);
					if ( path_length == 1 ) { //only count directly linked nodes 
						float mass_difference = peaks.at(j).mz - peaks.at(i).mz;
						
						//avoid the branches caused by isotopec ions
					for (k=0; k< (int) mass_dif.size(); k++)
					{
						if ( branch[k] == 0 ) 
						{
							// a new branch is found
							branch[k] = mass_difference;
							tree_size += 1 + tree_size_table[j]; //current link plus the size of the connected tree
							break;
						}
						else
						{	    
							if ( branch[k] - mass_difference <= 3 * MZ_ERROR &&  branch[k] - mass_difference >= -3 * MZ_ERROR )
							{	      
								break;
							}
						}
					}
				}
			}

			tree_size_table[i] = tree_size;
			
			if ( tree_size > max_tree_size)
			{
				max_tree_size = tree_size;
				start_pos = i;
			}
		}

		//re-construct the score path and store it
		if ( max_tree_size > 0 || num_peaks > 0 ) 
		{
			mobj_score_path.clear();
			mobj_score_path.push_back(peaks.at(start_pos));
			j = start_pos;
			tree_size = max_tree_size - 1;
			for (i= start_pos+1; i < num_peaks; i++)
			{
				path_length = getLink(link_tbl, num_peaks, start_pos, i);
				if ( path_length > 0 )
				{
					mobj_score_path.push_back(peaks.at(i));
					j = i;
					tree_size--;
				}
			}
		}

		// Anoop March 25 Check
		if (abs(max_tree_size - ((int)mobj_score_path.size()-1)) > 5)
		{
			bool debug = true ; 
		}


		#ifdef DEBUG
		assert ( charge <= MAX_CHARGE && charge > 0 );	
		#endif

		if (this->l_table != NULL)
		{
			delete(this->l_table);
		}
		this->l_table = link_tbl;

		delete(gap_too_wide);
		delete(tree_size_table);
		delete(branch);
		//return (float) max_tree_size;
		// ANOOP MARCH 26, 2012, testing alternate measure of score
		int score = 0 ; 
		if (max_tree_size > 0)
			 score = (int)mobj_score_path.size() -1 ; 
		
		return (float) score ; 

	}

				
		void Spectrum::initTable(int *tbl, int size)
		{
				int i, j;
				for (i=0; i<size; i++)
				{
					for (j=0; j<=i; j++) 
					{
						if ( i == j )
						//a link from one node to itself is 0
							tbl[(i*(i+1))/2+j] = 0;
						else
							tbl[(i*(i+1))/2+j] = -1;
					}
				}
		}
		
		
		float Spectrum::approxFindWithCharge(vector<float> &list, float value, int charge, float error_tolerance) 
		{
			//Find the element in the list matching the value with in the m/z error tolerance.
			//tolerance is manified by the charge of the ion. If there are two
			//matching elements, returns the one
			//with least difference. If no result is found, return
			//a negative value.
			float result = -1;
			for (int i=0; i< (int) list.size(); i++) 
			{
				float dif = fabsf(value - list[i]);
				if (  dif <= error_tolerance * charge ) 
				{
					if ( result < 0 || dif < fabsf(value - result) )
						result = list[i];
				}
			}
			return result;
		}

		void Spectrum::updateLinkTable(int *tbl, int size,int pk1, int pk2)
		{
			int p1, p2, i, j;
			int i_to_p1, p2_to_j;
			
			if ( pk1 < pk2 ) 
			{
				p1 = pk1;
				p2 = pk2;
			}
			else 
			{
				p1 = pk2;
				p2 = pk1;    
			}
			
			if ( getLink(tbl, size, p1, p2) < 1 )
			{
				setLink(tbl, size, p1, p2, 1);
				//update the rest of the table
				//for all the nodes which may have a link to p1 
				//(including p1)
				for (i=0; i<=p1; i++) 
				{
					//if there is a link from i to p1
					i_to_p1 = getLink(tbl, size, i, p1);
					if ( i_to_p1 >= 0 ) 
					{
						//for all the nodes which may be linked by p2 
						//(including p2)
						for (j=p2; j<size; j++) 
						{
							//if there is a link from p2 to j
							p2_to_j = getLink(tbl, size, p2, j);
							if ( p2_to_j >=0 )
							{
								if ( getLink(tbl, size, i, j) < i_to_p1 + 1 + p2_to_j )
									setLink(tbl, size, i, j, (i_to_p1 + 1 + p2_to_j));
							}
						}
					}
				}
			}
		}


		int Spectrum::getLink(int *tbl, int size, int pk1, int pk2)
		{
			assert ( pk1 < size && pk2 < size );
			if ( pk1 <= pk2 )
				return tbl[(pk2*(pk2+1))/2+pk1];
			else 
				return tbl[(pk1*(pk1+1))/2+pk2];    
		}
		
		void Spectrum::setLink(int *tbl, int size, int pk1, int pk2, int val)
		{
			assert ( pk1 < size && pk2 < size );
			if ( pk1 <= pk2 )
				tbl[(pk2*(pk2+1))/2+pk1] = val;
			else 
				tbl[(pk1*(pk1+1))/2+pk2] = val;    
		}


	}
}