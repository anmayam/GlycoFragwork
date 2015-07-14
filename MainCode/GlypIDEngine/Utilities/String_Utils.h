//---------------------------------------------------------------------------
#ifndef string_utilsH
#define string_utilsH

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

//---------------------------------------------------------------------------
using namespace std;

namespace Engine 
{
	namespace Utilities
	{
		class StringUtils
		{
		public:
	  static string getValue(const string& sLine);

	  static string getName(const string& sLine);
	  static bool getLine(istream& stream,
		string& sLine);

	  static void trim_left(string& sLine);
	  static void trim_right(string& sLine);
	  static void trim(string& sLine);
	  static string trim_left_copy(const string& sLine);
	  static string trim_right_copy(const string& sLine);
	  static string trim_copy(const string& sLine);
	  static void to_upper(string& sLine);
	  static void to_lower(string& sLine);
	  static string to_upper_copy(const string& sLine);
	  static string to_lower_copy(const string& sLine);
	  static void strToLower(char *str) ; 

	  static bool starts_with(const string& line, const string& substr);
	  static bool ends_with(const string& line, const string& substr);
	  static string int_to_str(int value);
	 static void StringUtils::str_to_char(std::string *strValue, const char *charValue) ; 	  
	  static string float_to_str(const char* formatStr, double value);
	  
	  static string float_to_str_trim(const char* formatStr, double value);

	  static bool is_trim_char(char c) ; 


	  static int str_to_int(const string& value);
	  static double str_to_float(const string& value);
	  static void printString(ostream& stream,
		const string& astr,
		const int lineLength);

	  static void remove_blank_line(vector<string>& lines);
	}; 
	}
}
//---------------------------------------------------------------------------
#endif


