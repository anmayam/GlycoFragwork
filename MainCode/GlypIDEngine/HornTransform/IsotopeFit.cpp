#include ".\IsotopeFit.h"
#include <iostream>
#include <time.h>
#include <fstream>
namespace Engine
{
	namespace HornTransform
	{
		IsotopeFit::IsotopeFit(void)
		{
			mdbl_cc_mass = 1.00727638 ; 		
			mbln_thrash = false ; 
			mbln_complete_fit = false ; 
			mbln_use_isotope_distribution_caching = false ; 
			mint_find_peak_cached = 0 ; 
			mint_find_peak_calc = 0 ; 
			Init() ; 
			mobj_averagine.SetElementalIsotopeComposition(mobj_isotope_dist.mobj_elemental_isotope_composition) ; 
		}

		void IsotopeFit::Init()
		{
			Reset() ; 
			mint_distribution_processing_time = 0 ; 
			mint_fit_processing_time = 0 ;
		}
		IsotopeFit::~IsotopeFit(void)
		{
			int i = 0 ; 
		}

		IsotopeFit::IsotopeFit(const IsotopeFit &fit)
		{
			// only copies settings not variables.
			this->mbln_complete_fit = fit.mbln_complete_fit ; 
			this->mbln_thrash = fit.mbln_thrash ;
			this->mdbl_cc_mass = fit.mdbl_cc_mass ; 
			this->mobj_averagine = fit.mobj_averagine ;
			this->mobj_isotope_dist = fit.mobj_isotope_dist ; 
			Init() ; 
		}

		IsotopeFit& IsotopeFit::operator=(const IsotopeFit &fit)
		{
			// only copies settings not variables.
			this->mbln_complete_fit = fit.mbln_complete_fit ; 
			this->mbln_thrash = fit.mbln_thrash ;
			this->mdbl_cc_mass = fit.mdbl_cc_mass ; 
			this->mobj_averagine = fit.mobj_averagine ;
			this->mobj_isotope_dist = fit.mobj_isotope_dist ; 
			Init() ; 
			return (*this) ;
		}

		bool IsotopeFit::FindPeak(double min_mz, double max_mz, double &mz_value, double &intensity, bool debug)
		{
			clock_t start_t = clock() ; 
			if (!mbln_last_value_was_cached)
			{
				bool found = mobj_isotope_dist.FindPeak(min_mz, max_mz, mz_value, intensity) ; 
				mint_find_peak_calc += (clock() - start_t) ; 
				return found ; 
			}

			mz_value = 0 ; 
			intensity = 0 ; 

			int num_pts = (int) mvect_distribution_mzs.size() ; 
			if (debug)
			{
				std::cerr<<" Num points = "<<num_pts<<" searching for min_mz = "<<min_mz<<std::endl ; 
			}

			int index = mobj_pk_index.GetNearestBinary(mvect_distribution_mzs, min_mz, 0, num_pts-1) ; 
			if (index >= num_pts)
				return false ; 

			if (mvect_distribution_mzs[index] > min_mz)
			{
				while(index > 0 && mvect_distribution_mzs[index] > min_mz)
					index-- ; 
			}
			else
			{
				while(index < num_pts && mvect_distribution_mzs[index] < min_mz)
					index++ ; 
				index-- ; 
			}

			int max_index = -1 ; 
			for ( ; index < num_pts ; index++)
			{
				double mz = mvect_distribution_mzs[index] ; 
				if (mz > max_mz)
					break ; 
				if (mz > min_mz)
				{
					if (mvect_distribution_intensities[index] > intensity)
					{
						max_index = index ; 
						intensity = mvect_distribution_intensities[index] ; 
					}
				}
			}
			if (max_index == -1)
			{
				mint_find_peak_cached += (clock() - start_t) ; 
				return false ; 
			}

			double X2 = mvect_distribution_mzs[max_index] ; 
			double X1 = X2 - 1.0 / mobj_isotope_dist.mint_points_per_amu ; 
			double X3 = X2 + 1.0 / mobj_isotope_dist.mint_points_per_amu ; 

			if (max_index > 0 && max_index < num_pts-1)
			{
				double Y1 = mvect_distribution_intensities[max_index - 1] ; 
				double Y2 = mvect_distribution_intensities[max_index] ; 
				double Y3 = mvect_distribution_intensities[max_index + 1] ; 

				// if the points are cached, these could be single sticks with surrounding 
				// points below background. To avoid that case, lets just check
				// and see if the differences in the theoretical mz values is as
				// expected or not. 
				if (mvect_distribution_mzs[max_index - 1] > X2 - 2.0/mobj_isotope_dist.mint_points_per_amu 
					&& mvect_distribution_mzs[max_index + 1] < X2 + 2.0/mobj_isotope_dist.mint_points_per_amu)
				{
					double d ;
					d = (Y2 - Y1) * (X3 - X2) ; 
					d = d - (Y3 - Y2) * (X2 - X1) ;
			   
					if (d == 0)
						mz_value = X2 ; 
					else
						mz_value = ((X1 + X2) - ((Y2 - Y1) * (X3 - X2) * (X1 - X3)) / d) / 2.0 ; 
				}
				else
				{
					mz_value = X2 ; 
					intensity = mvect_distribution_intensities[max_index] ; 
				}
				mint_find_peak_cached += (clock() - start_t) ; 
				return true ; 
			}

			mz_value = X2 ; 
			intensity = mvect_distribution_intensities[max_index] ; 
			mint_find_peak_cached += (clock() - start_t) ; 
			return true ; 
		}

		bool IsotopeFit::IsIsotopeLinkedDistribution(double min_threshold)
		{
			int num_isotopes = mvect_isotope_mzs.size()  ; 
			int size = mvect_isotope_intensities.size() ; 
		 
			for (int isotope_num = 3 ; isotope_num < num_isotopes; isotope_num ++) 
			{
				double mz = mvect_isotope_mzs[isotope_num] ; 
				double intensity = mvect_isotope_intensities[isotope_num]  ; 
				if ((intensity-min_threshold) > 50)
					return true  ; 
			}

			return false ; 
		}

		double IsotopeFit::GetFitScore(PeakProcessing::PeakData &pk_data, short cs, PeakProcessing::Peak &pk, IsotopeFitRecord &iso_record, double delete_intensity_threshold, 
				double min_theoretical_intensity_for_score, bool debug)
		{
			if (cs <= 0)
			{
				std::cerr<<"Negative value for charge state. "<<cs<<std::endl ; 
				exit(1) ; 
			}
			//initialize 
			double peak_mass = (pk.mdbl_mz - mdbl_cc_mass)* cs ; 
			// by now the cc_mass, tag formula and media options in Mercury(Isotope generation)
			// should be set. 
			if (debug)
				std::cerr<<"Getting isotope distribution for mass = "<<peak_mass<<" mz = "<<pk.mdbl_mz<<" charge = "<<cs<<std::endl ; 
			double resolution = pk.mdbl_mz / pk.mdbl_FWHM ;
			// DJ Jan 07 2007: Need to get all peaks down to delete interval so that range of deletion is correct. 
			//GetIsotopeDistribution(peak_mass, cs,  resolution, mvect_distribution_mzs, mvect_distribution_intensities, 
			//	min_theoretical_intensity_for_score, debug) ;
			GetIsotopeDistribution(peak_mass, cs,  resolution, mvect_distribution_mzs, mvect_distribution_intensities, 
				delete_intensity_threshold, debug) ;

			// Anoop April 9 2007: For checking if the distribution does not overlap/link with any other distribution
			// Beginnings of deisotoping correction			
			//bool is_linked = false ; 
			//is_linked =  IsIsotopeLinkedDistribution(delete_intensity_threshold) ; 

			double delta =  pk.mdbl_mz - mobj_isotope_dist.mdbl_max_peak_mz ; 
			double fit ; 
			if(debug)
				std::cerr<<"Going for first fit"<<std::endl ; 
			fit = FitScore(pk_data, cs, pk, delta, min_theoretical_intensity_for_score, debug) ; 
			if (debug)
				std::cerr<<"\tFirst fit  = "<<fit<<" Intensity = "<<pk.mdbl_intensity
					<<" Charge = "<<cs<<" FWHM = "<<pk.mdbl_FWHM<<" delta = "<<delta<<std::endl ; 

			if (!mbln_thrash )
			{
				iso_record.mdbl_fit = fit ; 
				iso_record.mdbl_mz = pk.mdbl_mz ; 
				iso_record.mdbl_average_mw = mobj_isotope_dist.mdbl_average_mw + delta * cs ; 
				iso_record.mdbl_mono_mw = mobj_isotope_dist.mdbl_mono_mw + delta*cs ; 
				iso_record.mdbl_most_intense_mw = mobj_isotope_dist.mdbl_most_intense_mw + delta*cs ;
				iso_record.mdbl_delta_mz = delta ; 
				return fit ; 
			}

			double p1fit =-1, m1fit = -1  ; 
			double Mpeak = mobj_isotope_dist.mdbl_max_peak_mz ; 
			PeakProcessing::Peak nxt_peak ; 

			double best_fit = fit ; 
			double best_delta = delta ; 
			double MaxY = pk.mdbl_intensity ; 

			for (double dd = 1.003/cs ; dd <= 10.03/cs ; dd+= 1.003/cs)
			{
				double mz, intensity ; 
				bool found_peak = FindPeak(Mpeak+dd - 0.2/cs, Mpeak+dd+ 0.2 /cs, mz, intensity, debug) ; 
				if (found_peak)
				{
					pk_data.FindPeak(pk.mdbl_mz - dd - pk.mdbl_FWHM, pk.mdbl_mz - dd + pk.mdbl_FWHM, nxt_peak) ; 
				}

				if (debug)
					std::cerr<<"\t\t Move by "<<dd ; 

				if (mz > 0 && nxt_peak.mdbl_mz > 0)
				{
					delta = pk.mdbl_mz - mz ;
					PeakProcessing::Peak current_peak_copy = pk ; 
					current_peak_copy.mdbl_intensity = nxt_peak.mdbl_intensity ; 
					fit = FitScore(pk_data, cs, current_peak_copy, delta, min_theoretical_intensity_for_score) ; 
					if (debug)
						std::cerr<<" isotopes. Fit ="<<fit<<" Charge = "<<cs<<" Intensity = "<<nxt_peak.mdbl_intensity<<" delta = "<<delta<<std::endl ; 			}
				else
				{
					if (debug)
						std::cerr<<" No peak found"<<std::endl ; 
					fit = best_fit + 0.001 ; 
				}
				// should allow thrashing backwards even if fit was greater than best_fit, if its still very small.
				// 26th February 2007 Deep Jaitly
				if (fit <= best_fit)
				{
					if (nxt_peak.mdbl_intensity > pk.mdbl_intensity)
						pk.mdbl_intensity = nxt_peak.mdbl_intensity ; 
					MaxY = pk.mdbl_intensity ; 
					best_fit = fit ; 
					best_delta = delta ; 
				}
				else
				{
					if (p1fit == -1)
						p1fit = fit ; 
					if (!mbln_complete_fit)
						break ; 
				}
			}

			for (double dd = 1.003/cs ; dd <= 10.03/cs ; dd+= 1.003/cs)
			{
				double mz, intensity ;
				mz = 0 ; 
				intensity = 0 ; 

				bool found_peak = FindPeak(Mpeak - dd - 0.2/cs, Mpeak - dd + 0.2 /cs, mz, intensity, debug) ; 
				if (found_peak)
				{
					pk_data.FindPeak(pk.mdbl_mz + dd - pk.mdbl_FWHM, pk.mdbl_mz + dd + pk.mdbl_FWHM, nxt_peak) ; 
				}
				if (debug)
					std::cerr<<"\t\t Move back by "<<dd ;
				if (mz > 0 && nxt_peak.mdbl_mz > 0)
				{
					delta = pk.mdbl_mz - mz ; 
					PeakProcessing::Peak current_peak_copy = pk ; 
					current_peak_copy.mdbl_intensity = nxt_peak.mdbl_intensity ; 
					fit = FitScore(pk_data, cs, current_peak_copy, delta, min_theoretical_intensity_for_score) ; 
					//fit = FitScore(pk_data, cs, nxt_peak.mdbl_intensity, delta) ; 
					if (debug)
						std::cerr<<" isotopes. Fit ="<<fit<<" Charge = "<<cs<<" Intensity = "<<nxt_peak.mdbl_intensity<<" delta = "<<delta<<std::endl ; 
				}
				else
				{
					fit = best_fit +0.001 ;
					if (debug)
						std::cerr<<"No peak found"<<std::endl ; 
				}
				if (fit <= best_fit)
				{
					if (nxt_peak.mdbl_intensity > pk.mdbl_intensity)
						pk.mdbl_intensity = nxt_peak.mdbl_intensity ; 
					MaxY = pk.mdbl_intensity ; 
					best_fit = fit ; 
					best_delta = delta ; 
				}
				else
				{
					if (m1fit == -1)
						m1fit = fit ; 
					if (!mbln_complete_fit)
						break ; 
				}
			}

			delta = best_delta ;

			iso_record.mdbl_fit = best_fit ; 
			iso_record.mshort_cs = cs ; 
			iso_record.mdbl_mz = pk.mdbl_mz ; 
			iso_record.mdbl_delta_mz = delta ; 
			iso_record.mdbl_average_mw = mobj_isotope_dist.mdbl_average_mw + delta * cs ; 
			iso_record.mdbl_mono_mw = mobj_isotope_dist.mdbl_mono_mw + delta*cs ; 
			iso_record.mdbl_most_intense_mw = mobj_isotope_dist.mdbl_most_intense_mw + delta*cs ; 
			//iso_record.mbln_flag_isotope_link = is_linked ; 
			return best_fit ; 
		}

		double IsotopeFit::GetFitScore(PeakProcessing::PeakData &pk_data, short cs, PeakProcessing::Peak &pk, 
			HornTransformTheoreticalProfile::MolecularFormula &formula, double delete_intensity_threshold, 
			double min_theoretical_intensity_for_score, bool debug)
		{
			if (cs <= 0)
			{
				std::cerr<<"Negative value for charge state. "<<cs<<std::endl ; 
				exit(1) ; 
			}
			if (debug)
				std::cerr<<"Getting isotope distribution for formula = "<<formula<<" mz = "<<pk.mdbl_mz<<" charge = "<<cs<<std::endl ; 

			double resolution = pk.mdbl_mz / pk.mdbl_FWHM ;

			mobj_isotope_dist.mdbl_cc_mass = mdbl_cc_mass ;
			mobj_isotope_dist.menm_ap_type = HornTransformTheoreticalProfile::GAUSSIAN ; 

			mvect_distribution_intensities.clear() ; 
			mvect_distribution_mzs.clear() ; 

			mobj_isotope_dist.CalculateDistribution(cs, resolution, formula, mvect_distribution_mzs, mvect_distribution_intensities, 
				delete_intensity_threshold, mvect_isotope_mzs, mvect_isotope_intensities, debug) ; 

			double delta =  pk.mdbl_mz - mobj_isotope_dist.mdbl_max_peak_mz ; 
			double fit ; 
			if(debug)
				std::cerr<<"Going for first fit"<<std::endl ; 
			fit = FitScore(pk_data, cs, pk, delta, min_theoretical_intensity_for_score, debug) ; 
			return fit ; 
		}

		void IsotopeFit::GetOptions(std::string &averagine_mf, std::string &tag_mf, bool &thrash_or_not, bool &complete_fit)
		{
			thrash_or_not = mbln_thrash ; 
			complete_fit = mbln_complete_fit ; 
			mobj_averagine.GetAveragineFormula(averagine_mf) ; 
			mobj_averagine.GetTagFormula(tag_mf) ; 
		}

		void IsotopeFit::SetOptions(std::string averagine_mf, std::string tag_mf, double cc_mass, bool thrash_or_not, bool complete_fit)
		{
			mobj_averagine.SetElementalIsotopeComposition(mobj_isotope_dist.mobj_elemental_isotope_composition) ; 
			mobj_averagine.SetAveragineFormula(averagine_mf) ; 
			mobj_averagine.SetTagFormula(tag_mf) ; 
			mbln_thrash = thrash_or_not ; 
			mbln_complete_fit = complete_fit ; 
			mdbl_cc_mass = cc_mass ; 
			mobj_mercury_cache.SetOptions(cc_mass, mobj_isotope_dist.mint_mercury_size) ; 
		}

		void IsotopeFit::GetZeroingMassRange(double &start_mz, double &stop_mz, double delta, double thresh, bool debug = false)
		{

			// this assumes that the last peak was not changed till now.
			int num_theoretical_points = (int) mvect_distribution_intensities.size() ; 
			double max_pk ; 
			double max_int = 0 ; 

			for(int i = 0 ; i < num_theoretical_points ; i++)
			{
				if (mvect_distribution_intensities[i] > max_int)
				{
					max_pk = mvect_distribution_mzs[i] ; 
					max_int = mvect_distribution_intensities[i] ; 
				}
			}


			for(int i = 0 ; i < num_theoretical_points ; i++)
			{
				double intensity = mvect_distribution_intensities[i] ;
				if (intensity > thresh)
				{
					start_mz = mvect_distribution_mzs[i] ; 
					stop_mz = start_mz ; 
					break ; 
				}
			}

			for(int i = num_theoretical_points-1 ; i > 0 ; i--)
			{
				if (mvect_distribution_intensities[i] > thresh)
				{
					stop_mz = mvect_distribution_mzs[i] ; 
					break ; 
				}
			}
			start_mz = start_mz + delta ; 
			stop_mz = stop_mz + delta ; 
			if (debug)
			{
				std::cerr<<"\t Start MZ for deletion ="<<start_mz<<" Stop MZ for deletion = "<<stop_mz<<std::endl ;
			}

			return ; 
		}

	

		inline void IsotopeFit::GetIsotopeDistribution(double most_abundant_mass, short charge, double resolution, 
			std::vector<double>&mzs, std::vector<double>&intensities, double min_theoretical_intensity, 
			bool debug)
		{
			double current_mz = most_abundant_mass/charge + mdbl_cc_mass ;
			double FWHM = current_mz / resolution ; 
			clock_t start_t = clock() ;

			if (!mbln_use_isotope_distribution_caching ||!mobj_mercury_cache.GetIsotopeDistributionCached(most_abundant_mass, charge, 
				FWHM, min_theoretical_intensity, mzs, intensities))
			{
				mbln_last_value_was_cached = false ;
				mobj_averagine.GetAverageFormulaForMass(most_abundant_mass, mobj_empirical_formula) ; 
				long lng_charge = (long) charge ; 

				mobj_isotope_dist.mdbl_cc_mass = mdbl_cc_mass ;
				mobj_isotope_dist.menm_ap_type = HornTransformTheoreticalProfile::GAUSSIAN ; 

				intensities.clear() ; 
				mzs.clear() ; 

				if (debug)
				{
					std::cerr<<"Getting distribution for chemical ="<<mobj_empirical_formula<<std::endl ;
				}
				mobj_isotope_dist.CalculateDistribution(charge, resolution, mobj_empirical_formula, mzs, intensities, min_theoretical_intensity, 
					mvect_isotope_mzs, mvect_isotope_intensities, debug) ; 
				if (mbln_use_isotope_distribution_caching)
				{
					mobj_mercury_cache.CacheIsotopeDistribution(most_abundant_mass, mobj_isotope_dist.mdbl_most_intense_mw, 
						mobj_isotope_dist.mdbl_mono_mw, mobj_isotope_dist.mdbl_average_mw, mobj_isotope_dist.mdbl_max_peak_mz, 
						charge, FWHM, mobj_isotope_dist.mdbl_mass_variance, mobj_isotope_dist.mint_points_per_amu, min_theoretical_intensity,
						mvect_isotope_mzs, mvect_isotope_intensities) ;
				}
			}
			else
			{
				mbln_last_value_was_cached = true ;
				//mzs, intensities are fetched, get the average and mono mw stats.
				mobj_isotope_dist.mdbl_average_mw = mobj_mercury_cache.mdbl_average_mw ; 
				mobj_isotope_dist.mdbl_mono_mw = mobj_mercury_cache.mdbl_mono_mw ; 
				mobj_isotope_dist.mdbl_most_intense_mw = mobj_mercury_cache.mdbl_most_intense_mw ; 
				mobj_isotope_dist.mdbl_max_peak_mz =  mobj_mercury_cache.mdbl_most_intense_mw/charge + mdbl_cc_mass - ELECTRON_MASS ; 
			}
			clock_t stop_t = clock() ; 
			mint_distribution_processing_time  += (stop_t - start_t) ; 
			
			return ; 
		}


		void IsotopeFit::SetCCMass(double mass)
		{
			mdbl_cc_mass = mass ; 
			mobj_mercury_cache.SetOptions(mass, mobj_isotope_dist.mint_mercury_size) ; 
		}

		void IsotopeFit::SetElementalIsotopeComposition(const HornTransformTheoreticalProfile::AtomicInformation &iso_comp)
		{
			mobj_isotope_dist.SetElementalIsotopeComposition(iso_comp) ; 
			mobj_averagine.SetElementalIsotopeComposition(iso_comp) ; 
		}
		void IsotopeFit::GetElementalIsotopeComposition(HornTransformTheoreticalProfile::AtomicInformation &iso_comp)
		{
			iso_comp = mobj_isotope_dist.mobj_elemental_isotope_composition ;
		}
	}
}