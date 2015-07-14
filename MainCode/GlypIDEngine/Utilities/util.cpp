///*****************************************************************************\
// * util.cpp
// *
// * some commonly used utilities
// *
// * author: Yin Wu (wuyin@indiana.edu)
// *
//\*****************************************************************************/
//
//

//
//#include "../MS2CIDScore/Spectrum.h"
#include "system.h"
#include "./util.h"
#include <vector>
#include <algorithm>
//
namespace Engine
{
	namespace Utilities
	{
//
//
//		template <class T>
//		void removeGenericDuplicates(std::vector<T> &t_list) 
//		{
//			if ( t_list.size() == 0 )
//				return;
//			else {
//				sort(t_list.begin(), t_list.end());
//				unique(t_list.begin(), t_list.end());
//			}
//		}
//
//		void removeDuplicates(vector<int> &t_list) 
//		{
//			removeGenericDuplicates(t_list);
//		}
//
//		void removeDuplicates(vector<float> &t_list) 
//		{
//			//removeGenericDuplicates(t_list);
//			if ( t_list.size() == 0 )
//				return;
//			else {
//				sort(t_list.begin(), t_list.end());
//				vector<float> unique_list;
//				float f = 0;
//				for (int i=0; i<(int) t_list.size(); i++) {
//					if ( unique_list.empty() || f != t_list[i] )
//					{
//						f = t_list[i];
//						unique_list.push_back(f);
//					}
//					else continue;
//				}
//				t_list.assign(unique_list.begin(), unique_list.end());
//				unique_list.clear();
//			}
//		}
//
//		bool approxEqual(float v1, float v2) {
//			//Note: not absolute value is needed here.
//			//the binary_search method of STL::algorithm 
//			//only requires a function to determine which
//			//value is less then current value.
//		  //if ( (v1 < v2 - MZ_ERROR * 1.5) )
//			if ( v1 < v2 - INTRA_SPECTRUM_ION_TRAP_MZ_ERROR )
//			return true;
//		  else return false;
//		}
//
//		//using the approxEqual function to find the element
//		//in the list matching the value. If there are two
//		//matching elements, approxEqual returns the one
//		//with least difference. If no result is found, return
//		//a negative value.
//		float approxFind(vector<float> &list, float value) 
//		{
//			float result = -1;
//			for (int i=0; i< (int) list.size(); i++) 
//			{
//				if ( (! approxEqual(list[i], value)) &&
//					(! approxEqual(value, list[i])) ) 
//				{
//					if ( result < 0 ||
//						( fabsf(value - list[i]) < 
//						fabsf(value - result) ) ) 
//						result = list[i];
//				}
//			}
//			return result;
//		}
//
//		//Find the element
//		//in the list matching the value with in the m/z error tolerance.
//		//tolerance is manified by the charge of the ion. If there are two
//		//matching elements, returns the one
//		//with least difference. If no result is found, return
//		//a negative value.
//		float approxFindWithCharge(vector<float> &list, float value, 
//								   int charge, float error_tolerance) 
//		{
//			float result = -1;
//			for (int i=0; i< (int) list.size(); i++) 
//			{
//				float dif = fabsf(value - list[i]);
//				if (  dif <= error_tolerance * charge ) 
//				{
//					if ( result < 0 ||
//						dif < fabsf(value - result) )
//						result = list[i];
//				}
//			}
//			return result;
//		}
//
//		bool peakSmaller(Peak p1, Peak p2) {
//		  switch ( Peak::sort_mode ) {
//		   case SORT_BY_MASS:
//			 return p1.mass < p2.mass;
//			 break;
//
//		   case SORT_BY_INTENSITY: default:
//			 return p1.intensity < p2.intensity;
//		  }
//		}
//		 
//		bool peakGreater(Peak p1, Peak p2) {
//		  switch ( Peak::sort_mode ) {
//		   case SORT_BY_MASS:
//			 return p1.mass > p2.mass;
//			 break;
//		   case SORT_BY_INTENSITY: default:
//			 return p1.intensity > p2.intensity;
//		  }
//		}
//


		float roundf(float val) 
		{
			//double factor = val<0 ?-0.5 : 0.5;
			double factor = 0.5;
			float near_rounded = floorf(val + (float) factor);
			return near_rounded;
		}


		float correlation(float* x, float* y, int delay)
		{
			int i,j ;
			double mx,my,sx,sy,sxy,denom,r;
			double trial1 = 0 ;
			int n= 7 ; 

			// Trying alternative
			for (i=0 ; i <n ; i++)
			{
				 j = i+delay ; 
				 if (j<0 || j >=n)
					 continue ; 
				 else
					 trial1 += x[i] *y[j] ; 
			}
   
			/* Calculate the mean of the two series x[], y[] */
			mx = 0;
			my = 0;   
			for (i=0;i<n;i++)
			{
				mx += x[i];
				my += y[i];
			}
			mx /= n;
			my /= n;

			/* Calculate the denominator */
			sx = 0;
			sy = 0;
			for (i=0;i<n;i++) 
			{
				  sx += (x[i] - mx) * (x[i] - mx);
				  sy += (y[i] - my) * (y[i] - my);
		   }
			denom = sqrt(sx*sy);

			/* Calculate the correlation series */
			sxy = 0;
			for (i=0;i<n;i++) 
			{
				 j = i + delay;
				 if (j < 0 || j >= n)
					continue;
				 else
					sxy += (x[i] - mx) * (y[j] - my);		
			 }
			 r = sxy / denom;
			 //return float(r) ; 

			 // Anoop this def seems better
			 return float(trial1) ; 

		}

		float max(float v1, float v2) {
		  if (v1 > v2)
			return v1;
		  else 
			return v2;
		}

		float min(float v1, float v2) {
		  if (v1 > v2)
			return v2;
		  else 
			return v1;
		}

		long NChooseK(int n, int k) 
		{
		  int factor1 = n;
		  int factor2 = n - k + 1;
		  int factor3 = k;
		  
		  assert( k >=0 && n >= k && n > 0 );

		  if ( k == 0 || n == k )
			return 1;
		  else if ( k > n - k )
			return NChooseK(n, n-k);

		  //the result is PIE(fact1, fact2)/PIE(fact3, 1)
		  double prod1 = 1;
		  double prod2 = 1;
		  for (int i= factor1; i >= factor2; i--)
			prod1 *= i;

		  for (int i= factor3; i > 1; i--)
			prod2 *= i;  

		  return (long) (prod1/prod2);
		}

		/*
		int abs(int v) {
		  if ( v > 0 )
			return v;
		  else return -v;
		} 
		*/

		//sort the source array and copy the
		//sorted result into dest. both src
		//and dest much be pre-allocated.
		//"order" must be either ASSENDING
		//or DESSENDING
		void sort_array(float *src, float *dest, 
				int length, int order) {
		  int i;
		  std::vector<float> v;

		  for (i=0; i< length; i++)
			v.push_back(src[i]);

		  std::sort(v.begin(), v.end());

		  if ( order == ASCENDING ) {
			i=0;
			while (! v.empty()) {
			  dest[i] = v.front();
			  v.erase(v.begin());
			  i++;
			}
		  }
		  else {
			i=length-1;
			while (! v.empty()) {
			  dest[i] = v.front();
			  v.erase(v.begin());
			  i--;
			}
		  }
		}

		long factorial(long n)
		{
			long int k;
			if(n==0)
				return(1);
			else
			{
				k=n*factorial(n-1);
			}
			return(k);
		}

//
//
//		int flip(int i) {
//		  if ( i == 0 )
//			return 1;
//		  else 
//			return 0;
//		}
//
//
//		//convert a spectrum to GPixel List
//		//intensity => value AND mass => x
//		//store the result in plist.
//		void SpectrumToPixelList(Spectrum *s, 
//					vector<GPixel> &plist) {
//		  vector<Peak>::iterator iter1;
//
//		  for ( iter1 = s->peaks.begin(); 
//			iter1 != s->peaks.end(); iter1++ ) {
//			plist.push_back(GPixel((*iter1).intensity,
//					   (*iter1).mass, 0, 0));
//		  }
//		}
//
//
//		//produce a random number between the
//		//value low and high inclusive. with uniform 
//		//distribution.
//		int rangedRandom(int low, int high) {
//		  double fraction = (1.0* rand())/RAND_MAX;
//		  return (int) (ceil((high - low + 1)* fraction) - 1 + low);
//		}
//
//
//		//produce a random number between 0 and 1
//		/*
//		double fractionRandom() {
//		  return (double) (1.0* rand())/RAND_MAX;
//		}
//		*/
//
//		//a function that computs the combination
//		//of N choose k
//		/*
//		long NChooseK(int n, int k) {
//		  int factor1 = n;
//		  int factor2 = n - k + 1;
//		  int factor3 = k;
//		  
//		  assert( k >=0 && n >= k && n > 0 );
//
//		  //the result is PIE(fact1, fact2)/PIE(fact3, 1)
//		  long prod1 = 1;
//		  long prod2 = 1;
//		  for (long i= (long) factor1; i >= (long) factor2; i--)
//			prod1 *= i;
//
//		  for (long i= (long) factor3; i > 1; i--)
//			prod2 *= i;  
//
//		  return prod1/prod2;
//		}
//		*/
//
//		//a function that computs the combination
//		//of N choose k
//		long NChooseK(int n, int k) {
//		  int factor1 = n;
//		  int factor2 = n - k + 1;
//		  int factor3 = k;
//		  
//		  assert( k >=0 && n >= k && n > 0 );
//
//		  if ( k == 0 || n == k )
//			return 1;
//		  else if ( k > n - k )
//			return NChooseK(n, n-k);
//
//		  //the result is PIE(fact1, fact2)/PIE(fact3, 1)
//		  double prod1 = 1;
//		  double prod2 = 1;
//		  for (int i= factor1; i >= factor2; i--)
//			prod1 *= i;
//
//		  for (int i= factor3; i > 1; i--)
//			prod2 *= i;  
//
//		  return (long) (prod1/prod2);
//		}
//
//
//		//convert a number to a format-controled 
//		//string
//		string number2Str(double num, const char *format ) {
//		  char buffer[MAX_LINE];
//		  sprintf(buffer, format, num);
//		  return string(buffer);
//		}
//
//
//		//VC++ does not support roundf() as linux does.
//		//this is a re-implementation of roundf
//		float roundf(float val) {
//			//double factor = val<0 ?-0.5 : 0.5;
//			double factor = 0.5;
//			float near_rounded = floorf(val + (float) factor);
//			return near_rounded;
//		}
//
//
//		//translate the combination of a ion/spectrum
//		//id and an file id into a floating number.
//		float idNumberToCombination(int id1, int id2) 
//		{
//			assert( id1 >= 0 && id2 >= 0 );
//			float f = id2 * 0.01;
//			return id1 + f;
//		}
//
//
//		//reverse-translate the combination of a ion/spectrum
//		//id and an file id into a floating number.
//		int getID1FromCombination(float comb) 
//		{
//			assert( comb >= 0 );
//			return floorf(comb);
//		}
//
//		//reverse-translate the combination of a ion/spectrum
//		//id and an file id into a floating number.
//		int getID2FromCombination(float comb) 
//		{
//			assert( comb >= 0 );
//			return roundf((comb - floorf(comb))*100);
//		}
//	}
//}
}
}
