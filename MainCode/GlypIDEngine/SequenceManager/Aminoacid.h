//---------------------------------------------------------------------------
// Modified from Tiger's proteomics library

#ifndef AminoacidH
#define AminoacidH

#include<string>
#include<vector>
#include<iostream>

using std::string;
using std::ostream;
using std::vector;


namespace Engine
{
	namespace SequenceManager
	{
		class TAminoacid 
		{
			public:
			  TAminoacid() : mOneName(0), mThreeName(""), mMonoMass(0.0), mNominalMass(0), mAverageMass(0.0), mDescription(""), 
				  mVisible(false)
			  { }

			  TAminoacid(const char aOneName, const string& aThreeName, const double aMonoMass,
				const double aAverageMass, const string& aDescription)
			  : mOneName(aOneName),
				mAverageMass(aAverageMass),
				mDescription(aDescription),
				mVisible(true)
			  {
				setThreeName(aThreeName);
				setMonoMass(aMonoMass);
			  }

			  TAminoacid(const TAminoacid& source)
			  : mOneName(source.mOneName),
				mThreeName(source.mThreeName),
				mMonoMass(source.mMonoMass),
				mAverageMass(source.mAverageMass),
				mNominalMass(source.mNominalMass),
				mDescription(source.mDescription),
				mVisible(source.mVisible)
			  { }

			  TAminoacid& operator=(const TAminoacid& source)
			  {
				if (this == &source)
				  return *this;

				mOneName = source.mOneName;
				mThreeName = source.mThreeName;
				mMonoMass = source.mMonoMass;
				mAverageMass = source.mAverageMass;
				mNominalMass = source.mNominalMass;
				mDescription = source.mDescription;
				mVisible = source.mVisible;

				return *this;
			  }

			  virtual ~TAminoacid()
			  { }

			  const string& threeName() const 
			  {
				return mThreeName;
			  }

			  char oneName() const 
			  {
				return mOneName;
			  }

			  double monoMass() const
			  {
				return mMonoMass;
			  }

			  double averageMass() const
			  {
				return mAverageMass;
			  }

			  int nominalMass() const
			  {
				return mNominalMass;
			  }

			  const string& description() const 
			  {
				return mDescription;
			  }

			  bool visible() const 
			  {
				return mVisible;
			  }

			  void resetMass(const double monoMass,	const double averageMass)
			  {
				mAverageMass = averageMass;
				setMonoMass(monoMass);
			  }

			  void initialize(const char oneName, const string& threeName, const double monoMass, 
				  const double averageMass, const string& description, const bool aVisible = true) 
			  {
				mOneName = oneName;
				mAverageMass = averageMass;
				setThreeName(threeName);
				setMonoMass(monoMass);
				mDescription = description;
				mVisible = aVisible;
			  }

			 void print(ostream& stream) const;

			 private:
				  void setMonoMass(const double aMonoMass) 
				  {
					mMonoMass = aMonoMass;
					mNominalMass = int(aMonoMass + 0.5);
				  }

                  void setThreeName(const string& aThreeName) 
				  { 
					  mThreeName = aThreeName.substr(0,3);
				  }

				  string  mThreeName;
				  char    mOneName;
				  double  mMonoMass;
				  double  mAverageMass;
				  int     mNominalMass;
				  string  mDescription;
				  bool    mVisible;
		};

		ostream& operator << (ostream& stream,TAminoacid& aminoacid);

		class TAminoacids
		{
			public:
				TAminoacids()
				{
				    initialize();
				}

				TAminoacids(const TAminoacids& source): mAminoacids(source.mAminoacids)
				{ }

				TAminoacids& operator=(const TAminoacids& source) 
				{
					if (&source == this)
						return *this;
        
					mAminoacids = source.mAminoacids;
    
					return *this;
				}

				virtual ~TAminoacids() 
				{
					mAminoacids.clear();
				}

				TAminoacid& operator[](const size_t index) 
				{
					return mAminoacids[index];
				}

				TAminoacid& operator[](const char amino) 
				{
					return mAminoacids[size_t(amino)];
				}

				size_t size() const 
				{
					return mAminoacids.size();
				}
				
				string information();

				void initialize();

				double averageResiduesMass(const string& sequence);
				double monoResiduesMass(const string& sequence);
				double averagePeptideMass(const string& sequence);
				double monoPeptideMass(const string& sequence);
			
			private:
				vector<TAminoacid> mAminoacids;
		};

		ostream& operator << (ostream& stream,TAminoacids& aminoacids);

		extern TAminoacids aminoacids;

	}
}

//---------------------------------------------------------------------------
#endif

