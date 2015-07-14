//---------------------------------------------------------------------------
// Modified from Tiger's code
#ifndef BioExceptionH
#define BioExceptionH

#include <string>
#include <stdexcept>
using namespace std;

namespace Engine
{
	namespace SequenceManager
	{
		class BioException:public runtime_error
		{
			public:
			  BioException(const string& msg): runtime_error(msg)
			  { };
		};
	}
}
//---------------------------------------------------------------------------
#endif
