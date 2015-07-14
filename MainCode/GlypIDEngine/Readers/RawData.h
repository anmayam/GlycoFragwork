// Written by Navdeep Jaitly for the Department of Energy (PNNL, Richland, WA)
// Copyright 2006, Battelle Memorial Institute
// E-mail: navdeep.jaitly@pnl.gov
// Website: http://ncrr.pnl.gov/software
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0

#pragma once
//#include "../Calibrations/Calibrator.h"
#include <vector>

namespace Engine {
 namespace Readers
{

	const int MAX_FNAME_LEN = 512 ;
	const int MAX_ERR_LEN = 512 ;

	enum FileType {FINNIGAN = 0, MZXMLRAWDATA, UNDEFINED} ;


	static int BACKGROUND_RATIO_FOR_TIC = 3 ;
	static double MIN_MZ = 400 ;
	static double MAX_MZ = 2000 ;

	bool IsDir(char *path) ;

	class   RawData
	{
		FileType menm_file_type ;
	protected:		
		virtual void GetSummedSpectra(std::vector <double> *bins, std::vector <double> *intensities, int scan, int scan_range,
			double min_mz, double max_mz, double mz_bin) ;
		virtual void GetRawData(std::vector <double> &vectMZs, std::vector<double> &vectIntensities, int scan, double min_mz,
			double max_mz) ;

	public:
		static const int MAX_SCAN_SIZE = 4*1024*1024 ;
		RawData(void)  ;
		virtual ~RawData(void) ;
		virtual const char *GetFileName() = 0 ;
		virtual FileType GetFileType() = 0 ;		
		virtual bool GetRawData(std::vector<double> *mzs, std::vector<double> *intensities, int scan_num)= 0 ;
		virtual bool GetRawData(std::vector<double> *mzs, std::vector<double> *intensities, int scan_num, int num_pts)= 0 ;		
		virtual void Load(char *file_n) = 0 ;
		virtual void Close() {} ;
		virtual int GetNumScans() = 0 ;
		virtual double GetScanTime(int scan_num) = 0 ;
		virtual int GetScanSize() = 0 ;
		virtual int GetNumScansLoaded() { return GetNumScans() ; } ;
		virtual short GetSpectrumType(int scan_num)
		{
			return 0 ;
		}
		virtual void GetScanDescription(int scan, char *description)
		{			
			strcpy(description,  "Scan #") ;
			_itoa(scan, &description[strlen(description)], 10) ;
		}
		virtual double GetSignalRange(int scan_num) = 0 ;
		virtual void GetTicFromFile(std::vector<double> *intensities, std::vector<double> *scan_times, bool base_peak_tic) = 0  ;
		virtual int GetNextScanNum(int current_scan_num) { return current_scan_num + 1 ; } ;
		virtual int GetFirstScanNum() { return 1 ; } ;
		virtual int GetLastScanNum() { return 1 ; } ;
		virtual int GetParentScan(int scan_num) = 0;
		virtual bool IsMSScan(int scan_num) = 0;

		virtual int GetMSLevel(int scan_num) = 0 ;
		virtual bool IsProfileScan(int scan_num) = 0;
		virtual double GetParentMz(int scan_num)  = 0;
		virtual double GetLevelNParentMz(int scan_num, int level_n)
		{
			return  0 ; 
		}
		virtual bool IsHCDScan(int scan_num) = 0 ;

		virtual bool IsFTScan(int scanNum)
		{
			return false ;
		}
		virtual bool IsETDScan(int scanNum)
		{
			return false ;
		}
		virtual bool IsCIDScan(int scanNum)
		{
			return false ;
		}
		virtual double  GetMonoMZFromHeader(int scan_num)
		{
			return 0 ;
		}
		virtual short GetMonoChargeFromHeader(int scan_num)
		{
			return 0 ;
		}		
		virtual void GetSummedSpectra(std::vector <double> *bins, std::vector <double> *intensities,
			int scan, int scan_range) ;
		virtual void GetSummedSpectra(std::vector <double> *mzs, std::vector <double> *intensities,
			int scan, int scan_range, double min_mz, double max_mz)  ;	
		virtual void GetSummedSpectra(std::vector <double> *mzs, std::vector <double> *intensities,
			std::vector<int> *scans_to_sum_over, double min_mz, double max_mz, double mz_bin, bool isMSMS)  ;		
		virtual double GetTICForScan(int scan_num)
		{
			return 0 ;
		}
	};
}
}
