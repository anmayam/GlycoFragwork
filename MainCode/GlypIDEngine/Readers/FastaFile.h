//---------------------------------------------------------------------------

#ifndef FastaFileH
#define FastaFileH

#include <string>
#include <fstream>
#include <vector>
#include "../SequenceManager/Sequence.h"

using std::string;
using std::ifstream;
using std::vector;

namespace Engine
{
	namespace Readers
	{
		class TFastaFile
		{
			public:
				TFastaFile() : mValid(false)  { }

				TFastaFile(const string& aFastaFilename): mValid(false) 
				{
					openFastaFile(aFastaFilename);
				}

				virtual ~TFastaFile() 
				{
					closeFastaFile();
				}

				bool valid() const
				{
					return mValid;
				}

				streampos tellg()
				{
					return mFile.tellg();
				}

				bool openFastaFile(const string& aFastaFilename);

				void closeFastaFile();

				bool readNextRecord(string &information, string &sequence,bool removeTabCharacterInDescription);

				static void readRecords(const string& aFastaFilename, vector<Engine::SequenceManager::Sequence>& aFastaSequences, bool removeTabCharacterInDescription);

				static bool ReadFastaFile(char *fasta_filename,  vector<Engine::SequenceManager::Sequence> &seq_list)  ; 

				
			private:
				string     mFastaFilename;
				ifstream   mFile;
				bool  mValid;
		};
	}
}

//---------------------------------------------------------------------------
#endif

