// Written by Navdeep Jaitly and Anoop Mayampurath for the Department of Energy (PNNL, Richland, WA)
// Copyright 2006, Battelle Memorial Institute
// E-mail: navdeep.jaitly@pnl.gov
// Website: http://ncrr.pnl.gov/software
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

#include "rawdata.h"
#include <math.h>
#include <windows.h>
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <io.h>
#include <fstream>
#include <errno.h>
#include <fcntl.h>
#include "..\Utilities\Interpolation.h"
#include "..\PeakProcessor\PeakIndex.h" 

namespace Engine 
{
	namespace Readers
	{

		bool IsDir(char *path)
		{
			int fp = _open(path, _O_RDONLY) ;  
			if (fp == -1)
			{
				int *error_num = _errno() ; 
				if (*error_num == EACCES)
				{
					// since we are opening file for read only access, if this happens, this is probably a folder.
					return true ; 
				}
				return false ; 
			}
			_close(fp) ; 
			return false ; 
		}

		RawData::RawData(void)
		{			
		}

		RawData::~RawData(void)
		{
		}
		
	
		void RawData::GetRawData(std::vector <double> &vectMZs, std::vector<double> &vectIntensities, int scan, double min_mz,
			double max_mz)
		{
			Engine::PeakProcessing::PeakIndex<double> peakIndex ; 
			if (max_mz <= min_mz)
				throw "max_mz should be greater than min_mz" ;

			vectMZs.clear() ; 
			vectIntensities.clear() ; 
			std::vector<double> allMZs ; 
			std::vector<double> allIntensities ; 
			GetRawData(&allMZs, &allIntensities, scan) ; 
			int numPts = allMZs.size() ; 
			if (numPts < 0)
			{
				std::cerr<<scan<<std::endl ; 
				//throw "mz_vector empty in GetRawData in RawData.cpp" ;  //Commented out May 16, 2012.  Don't know what this effect might have . Need to test.
				return ; 
			}
			int startIndex = 0 ; 
			int stopIndex = numPts ; 
			if (numPts  > 3)
			{
				startIndex = peakIndex.GetNearestBinary(allMZs, min_mz,0, numPts-1) ; 
				stopIndex = peakIndex.GetNearestBinary(allMZs, max_mz,0, numPts-1) ; 
				if ((stopIndex - startIndex) <= 1) //nothing in this m/z range
					return  ; 
				vectMZs.insert(vectMZs.begin(), &allMZs[startIndex], &allMZs[stopIndex]) ; 
				vectIntensities.insert(vectIntensities.begin(), &allIntensities[startIndex], &allIntensities[stopIndex]) ; 
				allMZs.clear() ; 
				allIntensities.clear() ; 
			}
			else
			{
				vectMZs.insert(vectMZs.begin(), allMZs.begin(), allMZs.end()) ;
				vectIntensities.insert(vectIntensities.begin(), allIntensities.begin(), allIntensities.end()); 
				allMZs.clear() ; 
				allIntensities.clear() ; 
				
			}
		}

		void RawData::GetSummedSpectra(std::vector <double> *mzs, std::vector <double> *intensities, int scan_num, int scan_range) 
		{			
			Engine::Utilities::Interpolation interpolator ; 
			// Get scan info
			GetRawData(mzs, intensities, scan_num) ; 
			int num_mzs = mzs->size() ; 
			const int numZeroFills = 3 ;
			if (num_mzs > numZeroFills)
			{
				try
				{
					interpolator.ZeroFillMissing(*mzs, *intensities, numZeroFills) ;             			
					num_mzs = mzs->size() ; 
					//recheck
                    if (num_mzs > numZeroFills)
					{
						int demarcationPoint = (int)((mzs->size()*1.0)/1.414) ; 
						double mzBin = (*mzs)[demarcationPoint+1] - (*mzs)[demarcationPoint] ; 
						double minMZBin = (*mzs)[1]-(*mzs)[0] ; 
						// need to make sure this bins size never goes below 0.1.
						// we need a better way than this stupid idea of mine. 
						if (mzBin > 5*minMZBin)
							mzBin = minMZBin ; 
						else if (mzBin > 0.1)
							mzBin = 0.1 ; 
						
						//checking range
						double minMZ = mzs->at(0) ; 
						double maxMZ = mzs->at(num_mzs-1) ; 
						mzs->clear() ; 
						intensities->clear() ; 					
						GetSummedSpectra(mzs, intensities, scan_num, scan_range, minMZ, maxMZ,mzBin) ; 
					}
				}
				catch (System::NullReferenceException *err)
				{
					mzs->clear() ; 
					intensities->clear() ; 	
					return ; 
				}
			}
		}

		void RawData::GetSummedSpectra(std::vector <double> *mzs, std::vector <double> *intensities, 
			std::vector<int> *scans_to_sum_over, double min_mz, double max_mz, double mz_bin, bool isMSMS)
		{

			//----- Function to sum over a set of scans indicated in the vector -------//
			std::vector<double> scan_mzs ; 
			std::vector<double> scan_intensities ; 
			std::vector<double> interpolatedIntensities ; 

			if (max_mz <= min_mz)
			{
				return ; 
			}

			int num_bins = (int)((max_mz - min_mz)/mz_bin) ; 
			double mz = min_mz ; 
			// initialize the vector
			try
			{
				for (int bin_num = 0 ; bin_num < num_bins ; bin_num++)
				{
					mzs->push_back(mz) ; 
					intensities->push_back(0) ;
					mz+=mz_bin ; 
				}
			}
			catch (System::NullReferenceException *err)
			{	
				mzs->clear() ; 
				intensities->clear() ; 
				return ; 
			}

			try
			{
				// For every scan
				for (int i = 0 ; i < scans_to_sum_over->size(); i++)
				{
					int currentScan = (*scans_to_sum_over)[i]	;
					scan_mzs.clear() ;
					scan_intensities.clear() ; 	
					interpolatedIntensities.clear() ; 

					// Read in scans
					GetRawData(scan_mzs, scan_intensities, currentScan, min_mz, max_mz) ; 
					if ((isMSMS) && (scan_intensities.size() > 0))
					{
						// MSMS does not require interpolation
						// So simply add the scan intensities to the vector						
						for (int peak_index = 0 ; peak_index < (int)scan_intensities.size(); peak_index++)
						{
							int bin_num = (int) ((scan_mzs[peak_index] - min_mz)/ mz_bin);	
							if ((scan_mzs[peak_index] >= min_mz) && (scan_mzs[peak_index] <= max_mz))
								(*intensities)[bin_num] = (*intensities)[bin_num] + scan_intensities[peak_index] ;									
						}
					}
					else
					{
						// TO be done
					}
				}
			}
			catch (System::NullReferenceException *err)
			{					
				mzs->clear() ; 
				intensities->clear() ; 
				return ; 
			}
		}

		void RawData::GetSummedSpectra(std::vector <double> *mzs, std::vector <double> *intensities, 
			int scan, int scan_range, double min_mz, double max_mz)
		{
			Engine::Utilities::Interpolation interpolator ; 

			// Get scan info
			GetRawData(*mzs, *intensities, scan, min_mz, max_mz) ;
			int num_mzs = mzs->size() ; 
			const int numZeroFills = 3 ; 
			if (num_mzs > numZeroFills)
			{			
				try
				{
					interpolator.ZeroFillMissing(*mzs, *intensities, numZeroFills) ;
					//recheck
					num_mzs = mzs->size() ; 
					if (num_mzs > numZeroFills)
					{
						int demarcationPoint = (int) ((mzs->size()*1.0)/1.414) ; 
						double mzBin = (*mzs)[demarcationPoint+1] - (*mzs)[demarcationPoint] ; 
						double minMZBin = (*mzs)[1]-(*mzs)[0] ; 

						// need to make sure this bins size never goes below 0.1.
						// we need a better way than this stupid idea of mine. 
						if (mzBin > 5*minMZBin)
							mzBin = minMZBin ; 
						else if (mzBin > 0.1)
							mzBin = 0.1 ; 

					
						mzs->clear() ; 
						intensities->clear() ; 				
						GetSummedSpectra(mzs, intensities, scan, scan_range, min_mz, max_mz, mzBin) ; 
					}
				}
				catch (System::NullReferenceException *err)
				{
					mzs->clear() ; 
					intensities->clear() ; 	
					return ; 
				}
			}
		}
		
		void RawData::GetSummedSpectra(std::vector <double> *mzs, std::vector <double> *intensities, int scan, int scan_range,
		double min_mz, double max_mz, double mz_bin)
		{
			int minDatasetScan = GetFirstScanNum() ; 
			int maxDatasetScan = GetLastScanNum() ; 

			std::vector<double> scan_mzs ; 
			std::vector<double> scan_intensities ; 
			std::vector<double> interpolatedIntensities ; 

			if (max_mz <= min_mz)
			{
				return ; 
			}

			int num_bins = (int)((max_mz - min_mz)/mz_bin) ; 
			double mz = min_mz ; 
			try
			{
				for (int bin_num = 0 ; bin_num < num_bins ; bin_num++)
				{
					mzs->push_back(mz) ; 
					intensities->push_back(0) ;
					mz+=mz_bin ; 
				}
			}
			catch (System::NullReferenceException *err)
			{	
				mzs->clear() ; 
				intensities->clear() ; 
				return ; 
			}
	
			// first read current scan and move to the left
			int currentScan = scan ; 
			int numScansSummed = 0 ; 

			try
			{
				// numScansSummed needs to be 1 + the scan range because we are summing the first 
				// scan here. 
				while(currentScan >= minDatasetScan && numScansSummed < scan_range +1 )
				{
					if (!IsMSScan(currentScan))
					{
						currentScan-- ; 
						continue ; 
					}

					scan_mzs.clear() ;
					scan_intensities.clear() ; 	
					interpolatedIntensities.clear() ; 

					GetRawData(scan_mzs, scan_intensities, currentScan, min_mz, max_mz) ; 
					
					if (scan_intensities.size() <= 3) //Keeping the min number of data points same as number of Zero Fill
					{
						currentScan-- ; 
						numScansSummed++ ;
						continue;
					}
					
					Engine::Utilities::Interpolation scan_interpolator; 		
					scan_interpolator.ZeroFillMissing(scan_mzs, scan_intensities, 3) ;             			
					if (scan_intensities.size() <= 3) //re-check
					{
						currentScan-- ; 
						numScansSummed++ ;
						continue;
					}
					scan_interpolator.Spline(scan_mzs, scan_intensities, 0, 0 ) ; 	
					scan_interpolator.Splint(scan_mzs, scan_intensities, *mzs, interpolatedIntensities) ; 					

					double maxScanMz = scan_mzs[(int)scan_mzs.size()-1] ;
					for (int bin_num = 0 ; bin_num < num_bins ; bin_num++)
					{
						if (bin_num * mz_bin + min_mz >= maxScanMz)
							break ; 
						if (interpolatedIntensities[bin_num]>0)
							(*intensities)[bin_num] = (*intensities)[bin_num] + interpolatedIntensities[bin_num] ;									
					}
					currentScan-- ;
					numScansSummed++ ;
				}
				currentScan = scan +1 ; 
				numScansSummed = 0 ; 
				while(currentScan <= maxDatasetScan && numScansSummed < scan_range )
				{
					if (!IsMSScan(currentScan))
					{
						currentScan ++ ; 
						continue ; 
					}

					scan_mzs.clear() ;
					scan_intensities.clear() ; 	
					interpolatedIntensities.clear() ; 

					GetRawData(scan_mzs, scan_intensities, currentScan, min_mz, max_mz) ; 

					if (scan_intensities.size() <= 3) //Keeping the min number of data points same as number of Zero Fill
					{
						currentScan++ ; 
						numScansSummed++ ;
						continue;
					}
					
					Engine::Utilities::Interpolation scan_interpolator ; 
					scan_interpolator.ZeroFillMissing(scan_mzs, scan_intensities, 3) ;             			
					if (scan_intensities.size() <= 3) //re-check
					{
						currentScan++ ; 
						numScansSummed++ ;
						continue;
					}

					scan_interpolator.Spline(scan_mzs, scan_intensities, 0, 0 ) ; 				
					scan_interpolator.Splint(scan_mzs, scan_intensities, *mzs, interpolatedIntensities) ; 
					
					double maxScanMz = scan_mzs[(int)scan_mzs.size()-1]; 
					for (int bin_num = 0 ; bin_num < num_bins ; bin_num++)
					{
						if (bin_num * mz_bin + min_mz >= maxScanMz)
							break ; 
						if (interpolatedIntensities[bin_num]>0)
							(*intensities)[bin_num] = (*intensities)[bin_num] + interpolatedIntensities[bin_num] ;									
					}
					currentScan++ ;
					numScansSummed++ ;
				}

				scan_mzs.clear() ; 
				scan_intensities.clear() ; 
				interpolatedIntensities.clear() ; 	
			}
			catch (System::NullReferenceException *err)
			{					
				mzs->clear() ; 
				intensities->clear() ; 
				return ; 
			}			
		}
	
	}
}