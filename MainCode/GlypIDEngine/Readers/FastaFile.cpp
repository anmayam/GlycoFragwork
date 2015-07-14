//---------------------------------------------------------------------------


//#pragma hdrstop

#include "FastaFile.h"
#include "../Utilities/String_Utils.h"

namespace Engine
{
	namespace Readers
	{
		bool TFastaFile::openFastaFile(const string& aFastaFilename)
		{ 
			closeFastaFile();
			mFastaFilename = aFastaFilename;
			mFile.open(mFastaFilename.c_str());
			mValid = mFile ; 
			return mValid;
		}

		void TFastaFile::closeFastaFile() 
		{
			if (mValid)
				mFile.close();
			mValid = false;
		}

		bool TFastaFile::readNextRecord(string &information, //the information after '>'
			string &sequence,  bool removeTabCharacterInDescription)
		{
			if (!mValid)
				return false;
			string tmp;

			while (Engine::Utilities::StringUtils::getLine(mFile, tmp))
			{
			    if ((!tmp.empty()) && ('>' == tmp[0]))
					break;
			}

			if (mFile.eof())
				return false;
			
			information = tmp.substr(1,tmp.length() - 1) ;
			if (removeTabCharacterInDescription)
			{
				replace(information.begin(), information.end(), '\t', ' '); 
			}

			sequence.clear();
			while ((!mFile.eof()) && ('>' != mFile.peek()) && Engine::Utilities::StringUtils::getLine(mFile,tmp)) 
			{
				sequence += tmp;
			}
			return true;
		}

		void TFastaFile::readRecords(const string& aFastaFilename, vector<Engine::SequenceManager::Sequence>& aFastaSequences, bool removeTabCharacterInDescription)
		{
			string info;
			string seq;
			aFastaSequences.clear();
			TFastaFile fastaFile(aFastaFilename);
			while(fastaFile.readNextRecord(info, seq, removeTabCharacterInDescription))
			{
				aFastaSequences.push_back(Engine::SequenceManager::Sequence(info, seq));
			}
			fastaFile.closeFastaFile();
		}

		bool TFastaFile::ReadFastaFile(char *fasta_filename,  vector<Engine::SequenceManager::Sequence> &seq_list) 
		{
			if( fasta_filename == NULL )
			{
				seq_list.clear() ; 
				std::cerr<<"Invalid file name" ; 
				return false; 
			}
			string fn = string(fasta_filename);
			readRecords(fn, seq_list, true);
			return true ; 
		}

	}
}

//---------------------------------------------------------------------------
// Ignoring below for now :Anoop
//#pragma package(smart_init)

