#include <cmath>
#include <time.h>
using namespace std;

#include "Distribution.h"
#include "../Utilities/system.h"

namespace Engine
{
	namespace MS2CIDScoring
	{
		Distribution::Distribution()
		{
		  mean = 0;
		  variance = 0;
		  confidence_interval_resolution = 0.01;
		  max_x = 0;
		  sample_size = 0;
		  type = DIST_UNKNOWN;
		}


		int Distribution::sampleSize()
		{
		  return (int) this->sample_size;
		}


		void Distribution::updateConfidenceIntervalList(double bound)
		{
		  double conf, temp;
		  double x;

		  if ( confidence_interval_list.size() <= 0 ) {
			conf = 0;
			x = 0;
		  }
		  else {
			conf = confidence_interval_list.back().value;
			x = confidence_interval_list.back().x;
		  }

		  temp = (double) (1.0/sqrt(2*PI));
		  for (double i = x + confidence_interval_resolution; 
			   i <= bound + confidence_interval_resolution + EPSILON; 
			   i += confidence_interval_resolution) {
			conf += temp * expf((float) (-0.5*i*i)) * confidence_interval_resolution;
			confidence_interval_list.push_back(Engine::PeakProcessing::GPixel((float) conf, 
								  (float) i, 0, 0));
		  }
		}

		void Distribution::setConfidenceIntervalResolution(double r)
		{
		  if ( r > EPSILON ) {
			this->confidence_interval_resolution = r;
			confidence_interval_list.clear();
		  }
		}

		//give x, check up its p-value
		/*
		double Distribution::xToPValue(double x, int test_type) {

		  //convert the distribution into a standard 
		  //normal distribution
		  assert(variance >= 0);
		  double standard_deviation = sqrtf(variance);
		  double standard_x = (x - mean) / standard_deviation;

		  if ( x > max_x )
			max_x = x;

		  //flip the negative half into the positive half
		  if ( standard_x < 0 )
			standard_x = - standard_x;

		  //check if updating the current confidence
		  //interval list is needed
		  if ( confidence_interval_list.size() <= 0 )
			this->updateConfidenceIntervalList( standard_x +
							confidence_interval_resolution );
		  else if ( standard_x > confidence_interval_list.back().x ) 
			this->updateConfidenceIntervalList(standard_x);

		  assert(confidence_interval_resolution > 0);  
		  long i = (long) 
			floor(standard_x/confidence_interval_resolution) - 1;

		  if ( i < 0 )
			i = 0;

		  assert( i < (long) confidence_interval_list.size() );

		  for (; i < (long) confidence_interval_list.size(); i++) {
			if ( confidence_interval_list[i].x >= standard_x ) {
			  if ( test_type == ONE_SIDED ) {
			if ( x > mean )
			  //return 0.5 - confidence_interval_list[i].value;
			  return 0.5 - confidence_interval_list[i].value;
			else 
			  return 0.5 + confidence_interval_list[i].value;
			  }
			  else 
			return 1.0 - confidence_interval_list[i].value*2;
			}
		  }

		  return 0;
		}
		*/

		double Distribution::xToPValue(double x, int test_type)
		{

		  //convert the distribution into a standard 
		  //normal distribution in which mean= 0 and var = 1
		  assert(variance >= 0);
		  double standard_deviation = sqrtf(variance);
		  double standard_x = (x - mean) / standard_deviation;

		  if ( x > max_x )
			max_x = x;

		  //cumulative distribution function
		  //to achieve higher accuracy. (of the erf() function).
		  //try to do the integration of the smaller end
		  double cdf = 0;
			
		#ifdef UNIX_OS
			cdf = (1 + erf(standard_x*ONE_OVER_SQRT_TWO))/2;
		#endif

		#ifdef WIN_OS
			cout << "windows OS does not support erf() in math.h" << endl;
			assert(0);
		#endif



		  double p = -1;

		  if ( test_type == ONE_SIDED ) {
			p = 1 - cdf;
		  }
		  else {
			//TWO_SIDED
			if ( cdf > 0.5 )
			  p = (1 - cdf) * 2;
			else
			  p = cdf * 2;
		  }
		    
		  return p;
		}


		double Distribution::xToPValue2(double x, int test_type) 
		{

		  //convert the distribution into a standard 
		  //normal distribution in which mean= 0 and var = 1
		  assert(variance >= 0);
		  double standard_deviation = sqrtf(variance);
		  double standard_x = (x - mean) / standard_deviation;

		  return standard_x*ONE_OVER_SQRT_TWO;
		}


		//sample in one dimension (only the values of a gpixel)
		//update the sample size, mean and variance.
		void Distribution::sample1D(vector<Engine::PeakProcessing::GPixel>& s, 
						bool keep_samples)
		{
		  vector<Engine::PeakProcessing::GPixel>::iterator iter1;
		  double sum = 0;

		  mean = 0;
		  variance = 1;
		  sample_size = s.size();
		  
		  if ( sample_size <= 0 )
			return;

		  //compute sample mean
		  for (iter1=s.begin(); iter1 != s.end(); iter1++) {
			sum += (*iter1).value;
		  }
		  mean = sum/sample_size;

		  //compute sample variance
		  sum = 0;
		  for (iter1=s.begin(); iter1 != s.end(); iter1++) {
			sum += ((*iter1).value - mean)*((*iter1).value - mean);
		  }

		  variance = sum/sample_size;

		  if ( keep_samples ) {
			samples.assign(s.begin(), s.end());
			//GPixel::sort(samples);
		  }
		}

		//initialize the random seed. this
		//must be called before calling produceRand
		void Distribution::initSRand() 
		{
		  srand((unsigned) clock());
		}


		double Distribution::produceRand() 
		{
		  switch ( type ) {
		  case DIST_UNKNOWN:
			return unknownRand();
			break;
		  default: 
			return 0;
		  }
		}


		//the density function of the distribution
		//is unknown, so return a random number
		//according to the current samples. assuming that
		//RAND_MAX 
		double Distribution::unknownRand() 
		{
		  //get a random index in the samples list
		  int rand_idx = (int) ( ((double) rand())/ RAND_MAX ) * ((int) samples.size());

		  //get a random number from -1~1
		  double fraction = 2 * ((double) rand())/ RAND_MAX - 1.0;
		  double delta = 3.0;

		  if ( rand_idx < 0 )
			rand_idx = 0;
		  if ( rand_idx >= (int) samples.size() )
			rand_idx = samples.size() -1;

		  return samples[rand_idx].value + delta * fraction;
		}



		void Distribution::setSampleMean(double mean)
		{
		  this->mean = mean;
		}

		void Distribution::setSampleVariance(double variance) 
		{
		  this->variance = variance;
		}

		void Distribution::setSampleSize(double size)
		{
		  this->sample_size = sample_size;
		}
	}
}