#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <assert.h>
#include "../PeakProcessor/GPixel.h"




#define ONE_OVER_SQRT_TWO 0.707107

enum {DIST_UNKNOWN=1, DIST_NORMAL, DIST_BINOMIAL, DIST_UNIFORM};
enum {ONE_SIDED=10, TWO_SIDED};

namespace Engine
{
	namespace MS2CIDScoring
	{
		class Distribution 
		{
			public:
			  double mean, variance;
			  double sample_size;
			  float max_x;
			  vector<Engine::PeakProcessing::GPixel> samples;
			  vector<Engine::PeakProcessing::GPixel> confidence_interval_list;
			  int type;

			  Distribution();

			  //give x, check up its p-value
			  //test_type can be either ONE_SIDED or TWO_SIDED
			  double xToPValue(double x, int test_type);
			  double xToPValue2(double x, int test_type);
			 
			  void setConfidenceIntervalResolution(double r);

			  //sample in one dimension (only the values of a gpixel)
			  //update the sample size, mean and variance.
			  void sample1D(vector<Engine::PeakProcessing::GPixel>& s, bool keep_samples);

			  //generate a random number from current distribution
			  double produceRand();

			  //initialize the random seed.
			  void initSRand();

			  int sampleSize();

			  void setSampleMean(double mean);
			  void setSampleVariance(double variance);
			  void setSampleSize(double size);

			private:
			  //the list of the confidence interval of
			  //a normal distrubtion of 0 mean and 1 variance.
			  //it is used to calculate the confidence
			  //interval of the distributions.
			  //vector<GPixel> confidence_interval_list;

			  //the step increment of x.
			  double confidence_interval_resolution;

			  void updateConfidenceIntervalList(double bound);

			  double unknownRand();
		};


		//some environment variables for one MS scan
		class Environment
		{
			public:
			  int num_glycoforms, num_ms2, num_ms1_ion;
			  Distribution nglycan_score_dist;

			  Environment(int num_glyco, int num_ms2, int num_ms1_ion) {
				this->num_glycoforms = num_glyco;
				this->num_ms2 = num_ms2;
				this->num_ms1_ion = num_ms1_ion;
			  };

			  Environment() {};
		};
	}
}

#endif
