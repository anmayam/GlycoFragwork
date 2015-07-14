// Written by Navdeep Jaitly for the Department of Energy (PNNL, Richland, WA)
// Copyright 2006, Battelle Memorial Institute
// E-mail: navdeep.jaitly@pnl.gov
// Website: http://ncrr.pnl.gov/software
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

#include "ReaderFactory.h"
#include "FinniganRawData.h" 
#include "MZXmlRawData.h"


namespace Engine 
{
	namespace Readers
	{
		RawData* ReaderFactory::GetRawData(FileType file_type, char *file_name)
		{
			char *header_n = "header" ; 
			switch(file_type)
			{
				case FINNIGAN:
#ifdef XCALIBUR_INSTALLED
					FinniganRawData *finnigan_raw_data ; 
					finnigan_raw_data = new FinniganRawData() ;
					finnigan_raw_data->Load(file_name) ; 
					return finnigan_raw_data ; 
#endif 
					break ; 
				case MZXMLRAWDATA:
					MZXmlRawData *mzxml_raw_data ; 
					mzxml_raw_data = new MZXmlRawData() ; 
					mzxml_raw_data->Load(file_name) ; 
					return mzxml_raw_data ; 
					break ; 
				default:
					break ; 
			}
			return NULL ; 
		}

		void ReaderFactory::GetRawData(RawData **raw_data, FileType file_type)
		{
			char *header_n = "header" ; 
			*raw_data = NULL ; 
			switch(file_type)
			{
				case FINNIGAN:
#ifdef XCALIBUR_INSTALLED
					*raw_data = new FinniganRawData() ; 
					return ; 
#endif 
					break ; 
				case MZXMLRAWDATA:
					*raw_data = new MZXmlRawData() ; 
					return ; 
					break ; 
				default:
					break ; 
			}
			return  ; 
		}

		RawData* ReaderFactory::GetRawData(FileType file_type)
		{
			char *header_n = "acqu" ; 
			switch(file_type)
			{
			
				case FINNIGAN:
#ifdef XCALIBUR_INSTALLED
					FinniganRawData *finnigan_raw_data ; 
					finnigan_raw_data = new FinniganRawData() ;
					return finnigan_raw_data ; 
#endif 
					break ; 
				case MZXMLRAWDATA:
					MZXmlRawData *mzxml_raw_data ; 
					mzxml_raw_data = new MZXmlRawData() ; 
					return mzxml_raw_data ; 
					break ; 
				default:
					break ; 
			}
			return NULL ; 
		}
	}
}