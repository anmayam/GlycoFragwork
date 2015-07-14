//---------------------------------------------------------------------------
#pragma hdrstop
#include "String_Utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
//---------------------------------------------------------------------------
//namespace proteomics {

namespace Engine
{
	namespace Utilities
	{
		const int GapAa = 'A' - 'a';

		string 	StringUtils::getValue(const string& line)
		{
		  size_t ipos = line.find("=");
		  if (ipos == string::npos)
			return "";
		  else
			return line.substr(ipos+1,line.length());
		}

		void  StringUtils::strToLower(char *str) 
	{
		int i;
		if ( str == NULL )
			return;
		else 
		{
			for (i=0; str[i] != '\0'; i++)
				str[i] = tolower(str[i]);
		}
	}


		string StringUtils::getName(const string& line) {
	  size_t ipos = line.find("=");
	  if (ipos == string::npos)
		return "";
	  else
		return line.substr(0, ipos-1);
	}

		bool StringUtils::getLine(istream& stream,
	  string& sLine) {
	  sLine.clear();
	  if (std::getline(stream,sLine))  {
		if ((sLine.length() > 0) && (sLine[sLine.length()-1] == '\r'))
		  sLine.erase(sLine.length()-1);
		return true;
	  } else {
		return false;
	  }
	}

		bool StringUtils::is_trim_char(char c) {
	  return ' ' == c || '\t' == c || '\n' == c || '\r' == c;
	}

	void
	StringUtils::
	trim_left(string& sLine) {
	  while((sLine.size() > 0) && is_trim_char(sLine[0])){
		sLine.erase(0,1);
	  }
	}

	void
	StringUtils::
	trim_right(string& sLine) {
	  while((sLine.size() > 0) && is_trim_char(sLine[sLine.size()-1])){
		sLine.erase(sLine.size()-1,1);
	  }
	}

	void
	StringUtils::
	trim(string& sLine) {
	  trim_left(sLine);
	  trim_right(sLine);
	}

	string
	StringUtils::
	trim_left_copy(const string& sLine) {
	  string result = sLine;
	  trim_left(result);
	  return result;
	}

	string
	StringUtils::
	trim_right_copy(const string& sLine) {
	  string result = sLine;
	  trim_right(result);
	  return result;
	}

	string
	StringUtils::
	trim_copy(const string& sLine) {
	  string result = sLine;
	  trim(result);
	  return result;
	}

	void
	StringUtils::
	to_upper(string& sLine) {
	  for (size_t i = 0;i < sLine.length();i++) {
		if ('a' <= sLine[i] && sLine[i] <= 'z')
		  sLine[i] = sLine[i] + GapAa;
	  }
	}

	void
	StringUtils::
	to_lower(string& sLine) {
	  for (size_t i = 0;i < sLine.length();i++) {
		if ('A' <= sLine[i] && sLine[i] <= 'Z')
		  sLine[i] = sLine[i] - GapAa;
	  }
	}

	string
	StringUtils::
	to_upper_copy(const string& sLine) {
	  string result = sLine;
	  to_upper(result);
	  return result;
	}

	string
	StringUtils::
	to_lower_copy(const string& sLine) {
	  string result = sLine;
	  to_lower(result);
	  return result;
	}

	bool
	StringUtils::
	starts_with(const string& line, const string& substr) {
	  if (line.empty() || line.length() < substr.length()){
		return false;
	  }

	  return 0 == strncmp(line.c_str(), substr.c_str(), substr.size());
	}

	bool
	StringUtils::
	ends_with(const string& line, const string& substr) {
	  if (line.empty() || line.length() < substr.length()){
		return false;
	  }

	  const char* rightLine = &(line.c_str()[line.length() - substr.length()]);
	  return 0 == strncmp(rightLine, substr.c_str(), substr.size());
	}

	string
	StringUtils::
	int_to_str(int value) {
	  ostringstream format_message;
	  format_message << value;
	  return format_message.str();
	}

	string
	StringUtils::
	float_to_str(const char* formatStr, double value) {
	  char buf[256];
	  sprintf(buf, formatStr, value);
	  return string(buf);
	}

	string
	StringUtils::
	float_to_str_trim(const char* formatStr, double value) {
	  char buf[256];
	  sprintf(buf, formatStr, value);
	  return trim_copy(string(buf));
	}

	int
	StringUtils::
	str_to_int(const string& value) {
	  try {
		istringstream input_string(value.c_str());
		int result;
		input_string >> result;
		return result;
	  }
	  catch(...){
		throw runtime_error(value + " is not a valid integer.");
	  }
	}

	void StringUtils::str_to_char(std::string *strValue, const char *charValue)
	{
		charValue =strValue->c_str() ; 		
	}

	double
	StringUtils::
	str_to_float(const string& value) {
	  try {
		istringstream input_string(value.c_str());
		double result;
		input_string >> result;
		return result;
	  }
	  catch(...){
		throw runtime_error(value + " is not a valid integer.");
	  }
	}

	void
	StringUtils::
	printString(ostream& stream,
	  const string& astr,
	  const int lineLength) {
	  size_t ibegin = 0;
	  while (ibegin < astr.length()) {
		size_t ilength = lineLength;
		if (ibegin + ilength < astr.length()) {
		  stream << astr.substr(ibegin,ilength) << std::endl;
		  ibegin += ilength;
		}
		else {
		  ilength = astr.length() - ibegin;
		  stream << astr.substr(ibegin,ilength);
		  break;
		}
	  }
	}

	void
	StringUtils::
	remove_blank_line(vector<string>& lines) {
	  vector<string>::iterator iter;
	  while(lines.end() != (iter = find(lines.begin(), lines.end(), ""))){
		lines.erase(iter);
	  }
	}

	}
}
//---------------------------------------------------------------------------
#pragma package(smart_init)
