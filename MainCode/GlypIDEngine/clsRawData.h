// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath, Indiana University
// Modified from GlypID code written by Navdeep Jaitly and Anoop Mayampurath 
// for the Department of Energy(PNNL, Richland, WA)

#pragma once
#include "Readers\ReaderFactory.h"
#include "GlypIDEngineUtils.h"


#using <mscorlib.dll>

namespace GlypID
{
 namespace Readers
	{
		public __value enum FileType {FINNIGAN = 0, MZXMLRAWDATA, UNDEFINED} ; 

		public __gc class clsRawData
		{
			Engine::Readers::RawData __nogc *mobj_raw_data ; 
		

		public:
			clsRawData();
			clsRawData(System::String *file_name, FileType file_type);
			int GetFirstScanNum() ; 
			~clsRawData(void);

			__property GlypID::Readers::FileType get_FileType()
			{
				if (mobj_raw_data == NULL)
					return GlypID::Readers::FileType::UNDEFINED ; 
				return (GlypID::Readers::FileType) mobj_raw_data->GetFileType() ; 
			}

			__property System::String* get_FileName()
			{
				System::String *file_name = S"" ; 
				if (mobj_raw_data == NULL)
					return NULL ; 
				GlypID::Utils::GetStr(mobj_raw_data->GetFileName(), &file_name) ; 
				return file_name ; 
			}

			__property int get_PercentDone()
			{
				if (mobj_raw_data == NULL)
					return 0 ; 
				int num_scans = mobj_raw_data->GetNumScans() ; 
				if (num_scans != 0)
				{
					int percent_done = (100 * mobj_raw_data->GetNumScansLoaded())/ num_scans  ; 
					if (percent_done < 0)
						return 0 ; 
					if (percent_done > 100)
						return 100 ; 
					return percent_done ; 
				}
				return 0 ; 
			}

			__property System::String* get_StatusMessage()
			{
				if (mobj_raw_data == NULL)
					return S"" ; 
				int current_scan = mobj_raw_data->GetNumScansLoaded() ; 
				int num_scans = mobj_raw_data->GetNumScans() ; 
				return System::String::Concat(S"Processed :", System::Convert::ToString(current_scan), S" of ", System::Convert::ToString(num_scans), S" scans") ; 
			}

			short GetMonoChargeFromHeader(int scanNum)
			{
				return mobj_raw_data->GetMonoChargeFromHeader(scanNum) ; 
			}

			double GetMonoMzFromHeader(int scanNum)
			{
				return mobj_raw_data->GetMonoMZFromHeader(scanNum) ; 
			}

			bool IsFTScan(int scanNum)
			{
				return mobj_raw_data->IsFTScan(scanNum) ; 
			}

			bool IsMSScan(int scanNum)
			{
				return mobj_raw_data->IsMSScan(scanNum) ; 
			}

			bool IsHCDScan(int scanNum)
			{
				return mobj_raw_data->IsHCDScan(scanNum) ; 
			}
			bool IsCIDScan(int scanNum)
			{
				return mobj_raw_data->IsCIDScan(scanNum) ; 
			}

			bool IsProfileScan(int scanNum) 
			{
				return mobj_raw_data->IsProfileScan(scanNum) ; 
			}

			bool IsETDScan(int scanNum)
			{
				return mobj_raw_data->IsETDScan(scanNum) ; 
			}

			void LoadFile(char *file_name, GlypID::Readers::FileType file_type) ; 
			void LoadFile(char file_name __gc [], GlypID::Readers::FileType file_type) ; 
			void LoadFile(System::String *file_name, GlypID::Readers::FileType file_type) ; 

			void Close() ; 
			int GetNumScans() ;
			double GetScanTime(int scan_num) ; 		
			int GetMSLevel(int scan_num) ; 
			int GetScanSize() ; 
			short GetSpectrumType(int scan_num) ;			
			void GetTicFromFile(float  (&intensities) __gc [], float (&scan_times) __gc [], bool base_peak_tic); 
			void GetSpectrum(int scan_num, float  (&mzs) __gc [], float (&intensities) __gc []) ; 
			int GetParentScan(int scan_num);
			double GetParentMz(int scan_num);
			double GetLevelNParentMZ(int scan_num, int level_n);
			void GetMzsInRange(float (&in_mzs) __gc [], float (&in_intensities) __gc [], float (&out_mzs) __gc [], float (&out_intensities) __gc [], float central_value, float range);
			System::String *GetScanDescription(int scan_num) ; 			
			void GetSummedSpectra(int current_scan , int scan_range, float (&mzs) __gc[], float (&intensities) __gc[]) ; 
			void GetSummedSpectra(int start_scan, int stop_scan, double min_mz, double max_mz, float (&mzs) __gc[], float (&intensities) __gc[]) ; 
			void GetSummedSpectra(int msms_scans_to_sum_over __gc[], double min_mz, double max_mz, double mz_bin, float (&mzs) __gc[], float (&intensities) __gc[]) ; 



		};
	}
}