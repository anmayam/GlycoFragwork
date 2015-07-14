// Modified from Tiger's proteomics library
//---------------------------------------------------------------------------

#pragma hdrstop

#include "Aminoacid.h"
#include "../Utilities/String_Utils.h"

#include <iomanip>

//---------------------------------------------------------------------------

namespace Engine
{
	namespace SequenceManager
	{
		TAminoacids aminoacids;
		
		std::ostream& operator << (std::ostream& stream,TAminoacid& aminoacid) 
		{
			aminoacid.print(stream);
			return stream;
		}

		std::ostream& operator << (std::ostream& stream,TAminoacids& aminoacids) 
		{
			for (size_t i = 0;i < aminoacids.size();i++)
				if (aminoacids[i].visible())
					stream << aminoacids[i] << std::endl;
			return stream;
		}

//---------------------------------------------------------------------------

		void TAminoacid::print(std::ostream& stream) const 
		{
			stream << mOneName << ' '<< mThreeName<< Engine::Utilities::StringUtils::float_to_str("%10.4lf", mMonoMass)
				<< Engine::Utilities::StringUtils::float_to_str("%10.4lf", mAverageMass)
				<< Engine::Utilities::StringUtils::int_to_str(mNominalMass)
				<< ' '<< mDescription;
		}

		void TAminoacids::initialize() 
		{
			  mAminoacids.clear();
			  mAminoacids.resize(128);
			  mAminoacids['G'].initialize('G', "Gly", 57.02147, 57.05, "Glycine C2H3NO");
			  mAminoacids['A'].initialize('A', "Ala", 71.03712, 71.08, "Alanine C3H5NO");
			  mAminoacids['S'].initialize('S', "Ser", 87.03203, 87.08, "Serine C3H5NO2");
			  mAminoacids['P'].initialize('P', "Pro", 97.05277, 97.12, "Proline C5H7NO");
			  mAminoacids['V'].initialize('V', "Val", 99.06842, 99.13, "Valine C5H9NO");
			  mAminoacids['T'].initialize('T', "Thr", 101.04768, 101.10, "Threonine C4H7NO2");
			  mAminoacids['C'].initialize('C', "Cys", 103.00919, 103.14, "Cysteine C3H5NOS");
			  mAminoacids['I'].initialize('I', "Ile", 113.08407, 113.16, "Isoleucine C6H11NO");
			  mAminoacids['L'].initialize('L', "Leu", 113.08407, 113.16, "Leucine C6H11NO");
			  mAminoacids['N'].initialize('N', "Asn", 114.04293, 114.10, "Asparagine C4H6N2O2");
			  mAminoacids['D'].initialize('D', "Asp", 115.02695, 115.09, "Aspartic acid C4H5NO3");
			  mAminoacids['Q'].initialize('Q', "Gln", 128.05858, 128.13, "Glutamine C5H8N2O2");
			  mAminoacids['K'].initialize('K', "Lys", 128.09497, 128.17, "Lysine C6H12N2O");
			  mAminoacids['E'].initialize('E', "Glu", 129.04260, 129.12, "Glutamic acid C5H7NO3");
			  mAminoacids['M'].initialize('M', "Met", 131.04049, 131.19, "Methionine C5H9NOS");
			  mAminoacids['H'].initialize('H', "His", 137.05891, 137.14, "Histidine C6H7N3O");
			  mAminoacids['F'].initialize('F', "Phe", 147.06842, 147.18, "Phenylalanine C9H9NO");
			  mAminoacids['R'].initialize('R', "Arg", 156.10112, 156.19, "Arginine C6H12N4O");
			  mAminoacids['Y'].initialize('Y', "Tyr", 163.06333, 163.18, "Tyrosine C9H9NO2");
			  mAminoacids['W'].initialize('W', "Trp", 186.07932, 186.21, "Tryptophan C11H10N2O");
		}

		double TAminoacids::averageResiduesMass(const string& sequence)
		{
			  double dRes(0.0);
			  for (size_t i = 0;i < sequence.length();i++)
				dRes += mAminoacids[sequence[i]].averageMass();
			  return dRes;
		}

		double TAminoacids::monoResiduesMass(const string& sequence) 
		{
			double dRes(0.0);
			for (size_t i = 0;i < sequence.length();i++)
			    dRes += mAminoacids[sequence[i]].monoMass();
			return dRes;
		}

		const double MONO_H = 1.007825035;
		const double MONO_O = 15.99491463;
		const double AVERAGE_H = 1.00794;
		const double AVERAGE_O = 15.9994;
		
		double	TAminoacids::averagePeptideMass(const string& sequence) 
		{
			double dRes(AVERAGE_H * 2 + AVERAGE_O);
			for (size_t i = 0;i < sequence.length();i++)
				dRes += mAminoacids[sequence[i]].averageMass();
			return dRes;
		}

		double TAminoacids::monoPeptideMass(const string& sequence) 
		{
			  double dRes(MONO_H * 2 + MONO_O);
			  for (size_t i = 0;i < sequence.length();i++)
				dRes += mAminoacids[sequence[i]].monoMass();
			  return dRes;
		}

	}
}

#pragma package(smart_init)

