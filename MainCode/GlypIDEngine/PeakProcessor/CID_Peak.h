#ifndef PEAK_H
#define PEAK_H

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <string>
#include <math.h>
using namespace std;

#define WIN_OS

#define HEXNAC_MONO 203.0794
#define HEXNAC_AVG 203.195

#define HEX_MONO 162.0528 
#define HEX_AVG 162.1424 

#define NEUAC_MONO 291.0954
#define NEUAC_AVG 291.2579

//e.g. fucose
#define DEHEX_MONO 146.0579
#define DEHEX_AVG 146.1463

#define MAX_NUM_PEAKS 100
#define FT_MZ_ERROR 0.015
#define FT_PPM_ERROR 5
#define ISOFORM_DIFF 1.00335
#define MAX_CHARGE 4
#define CHARGE_3_LIMIT 1000
#define CHARGE_4_LIMIT 1300
#define PROTON 1.00727638
#define IODOACEDAMIDE_MONO 57.02146 
#define IODOACEDAMIDE_AVG 57.0520 

//when changing the value of MAX_NUM_SUB_MS, 
//clean the .o files and re-make from empty
//to avoid seg fault.
#define MAX_NUM_SUB_MS 20

namespace Engine
{
	namespace PeakProcessing
	{	

		//******** data structure *******/
		enum SORT_MODE {SORT_BY_MZ=100, SORT_BY_INTENSITY};
		enum SORT_ORDER {ASCENDING = 150, DESCENDING};
		
		class CID_Peak 
		{
			public:
				float mz, intensity, original_intensity, original_mass;
				float glycan_mass; //the glycan that is associated with this peak
				int glycan_charge; //charge level of the associated glycan
				float ms2_score, b_ion_in_ms2;
				int ms2_scanheader, charge, intensity_rank;		
				
				
				float flag;
				int parent_isoform_id;				

				void init();
				CID_Peak(float mz, float intensity);
				CID_Peak(float mz, float intensity, int charge);
				CID_Peak(float mz, float intensity, float glycan_mass, int glycan_cs);
				
				 

				CID_Peak(void);
				~CID_Peak(void);
				
				float originalMass();
				


				bool isEmpty();
				void setEmpty();		

		};
	}
}

#endif