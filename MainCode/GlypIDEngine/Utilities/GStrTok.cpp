/*****************************************************************************\
 * GStrTok.cpp
 *
 * A simple string tokenizer tool. compatible on all c++ versions.
 *
 * author: Yin Wu (wuyin@indiana.edu)
 *
\*****************************************************************************/

#define GSTRTOK_CPP

#include ".\system.h"
#include "GStrTok.h"

namespace Engine
{
	namespace Utilities
	{

		GStrTok::~GStrTok() {
		  clear();
		}

		GStrTok::GStrTok(char* str, const char* delimitors) {
		  char *tok, *str_cpy, *tmp;
		  int len;

		  str_cpy = tmp = NULL;

		  if ( str == NULL || delimitors == NULL)
			clear();
		  else {
			//make a copy of the str (so that the original
			//one will not be destroyed
			len = strlen(str) + 1;
			str_cpy = (char*) malloc(sizeof(char)*len);
		#ifdef DEBUG
			assert( str_cpy != NULL );
		#endif
			memcpy(str_cpy, str, sizeof(char)*len);

			for (tok = strtok(str_cpy, delimitors);
			 tok != NULL; tok = strtok(NULL, delimitors)) {
			  len = strlen(tok) + 1;
			  tmp = (char*) malloc(sizeof(char)*len);
		#ifdef DEBUG
			  assert( tmp != NULL );
		#endif
			  memcpy(tmp, tok, sizeof(char)*len);
		      
			  tokens.push_back(tmp);
			}
		  }

		  free(str_cpy);
		}


		int GStrTok::size() {
		  return (int) tokens.size();
		}



		void GStrTok::clear() {
		  vector<char*>::iterator iter1;

		  for ( iter1 = tokens.begin();
			iter1 != tokens.end(); iter1++ ) {
			free((*iter1));
		  }

		  tokens.clear();
		}


		void GStrTok::tokenize(char* str, const char* delimitors) {
		  char *tok, *str_cpy, *tmp;
		  int len;

		  clear();
		  str_cpy = tmp = NULL;

		  if ( str == NULL || delimitors == NULL)
			clear();
		  else {
			//make a copy of the str (so that the original
			//one will not be destroyed
			len = strlen(str) + 1;
			str_cpy = (char*) malloc(sizeof(char)*len);
		#ifdef DEBUG
			assert( str_cpy != NULL );
		#endif
			memcpy(str_cpy, str, sizeof(char)*len);

			for (tok = strtok(str_cpy, delimitors);
			 tok != NULL; tok = strtok(NULL, delimitors)) {
			  len = strlen(tok) + 1;
			  tmp = (char*) malloc(sizeof(char)*len);
		#ifdef DEBUG
			  assert( tmp != NULL );
		#endif
			  memcpy(tmp, tok, sizeof(char)*len);
		      
			  tokens.push_back(tmp);
			}
		  }

		  free(str_cpy);
		}

		//returns the token at position i
		//if i excceeds the range of tokens,
		//return NULL. 
		//NOTE: this function does not allocate
		//new memory space for the returned token
		//so, do not write to the returned token
		//directly!!!
		char* GStrTok::at(int i) {

		  if ( i < 0 || i > size() -1 )
			return NULL;
		  else 
			return tokens.at(i);
		}
	}	
}
