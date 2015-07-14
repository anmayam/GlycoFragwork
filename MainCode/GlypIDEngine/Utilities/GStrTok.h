/*****************************************************************************\
 * GStrTok.h
 *
 * A simple string tokenizer tool. compatible on all c++ versions.
 *
 * author: Yin Wu (wuyin@indiana.edu)
 *
\*****************************************************************************/

#ifndef GSTRTOK_H
#define GSTRTOK_H
#include <vector>

namespace Engine
{
	namespace Utilities
	{
		class GStrTok 
		{
			public:
			  vector<char*> tokens;

			  //GStrTok();
			  ~GStrTok();
			  GStrTok(char* str, const char* delimitors);
			  int size();  
			  void clear();
			  void tokenize(char* str, const char* delimitors);
			  char* at(int i);
		};
	}
}

#endif
