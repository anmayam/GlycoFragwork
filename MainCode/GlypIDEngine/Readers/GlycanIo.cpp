///*****************************************************************************\
// * io.cpp
// *
// * File reader/writer
// *
// * author: Yin Wu (wuyin@indiana.edu)
// *
//\*****************************************************************************/
//

// Modified by Anoop for GlypID 2.0

#include "GlycanIo.h"

#include "../Utilities/String_Utils.h"
#include "../Utilities/system.h"
#include "../Utilities/GStrTok.h"
#include <vector>

namespace Engine
{
	namespace Readers
	{
		GlycanIo::~GlycanIo(void)
		{		
			glycan_compositions_in_search.clear() ; 			
		};

		GlycanIo::GlycanIo(void)
		{
			min_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_TOTAL] = 0;
			min_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_DEHEX] = 0;
			min_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_HEX] = 0;
			min_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_HEXNAC] = 0;
			min_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_NEUAC] = 0;
			min_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_NEUGC] = 0;
			min_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_MOD_SULFATE] = 0 ; 
			
			max_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_TOTAL] = 20;
			max_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_DEHEX] = 20;
			max_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_HEX] = 20;
			max_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_HEXNAC] = 20;
			max_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_NEUAC] = 20;
			max_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_NEUGC] = 20;
			max_glycan_composition.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_MOD_SULFATE] = 0 ; 

		};		

		bool GlycanIo::LoadGlycanListFromFile(char *glycanFileName, std::vector<Engine::GlycanTheoretical::GlycanComposition> &glycan_composition_list)
		{	
			Engine::GlycanTheoretical::GlycanComposition gc;
			
						
			if(glycanFileName == NULL )
			{				
				std::cerr<<"Invalid file name" ; 
				return false; 
				glycan_composition_list.clear() ; 
			}

			//read from glycan library
			std::vector<Engine::GlycanTheoretical::GlycanComposition> gc_list ; 
			ifstream *in_f = new ifstream(glycanFileName);
			readGlycanCompositionLib(in_f, gc_list);		


			// aug 2012 change putting in code to remove duplicate
			std::vector<Engine::GlycanTheoretical::GlycanComposition> distinct_gc_list ; 
			for (int i = 0 ; i < (int) gc_list.size(); i++)
			{
				string &str1 = gc_list[i].getCompositionString() ; 
				bool is_duplicated = false ; 
				for (int j= 0 ; j <(int) distinct_gc_list.size(); j++)
				{
					string &str2 = distinct_gc_list[j].getCompositionString() ; 					
					if (str1.compare(str2) ==0)
					{
						is_duplicated = true ; 
						break ; 
					}
				}

				if (! is_duplicated)
					distinct_gc_list.push_back(gc_list[i]) ; 
			}
			
			//validate the glycan compositions by user defined filter //Anoop disabled						
			glycan_compositions_in_search.clear();			
			int numGlycans = 0 ; 
			numGlycans = (int) distinct_gc_list.size(); 

			for (int i=0; i< numGlycans ; i++)
			{  
				gc = distinct_gc_list[i];

				gc.glycan_accurate_mass = gc.accurateMass() ; 
				gc.glycan_avg_mass = gc.averageMass() ; 

				//Anoop commenting this out for nowif ( gc.validate(min_glycan_composition, max_glycan_composition,5, 20, glycopeptide_matching_glycan_type) )
				//{
					glycan_compositions_in_search.push_back(gc);
					glycan_composition_list.push_back(gc) ; 

					
				//}
			}

			//convert the map into a sorted vector.
			/*vector<Engine::GlycanTheoretical::GlycanComposition>::iterator iter;
			vector<Engine::GlycanTheoretical::GlycanComposition>::iterator iter2;
			glycan_composition_list.clear();
			for (iter = glycan_compositions_in_search.begin(); iter != glycan_compositions_in_search.end(); iter++ )
			{
				for (iter2 = glycan_composition_list.begin(); iter2 != glycan_composition_list.end(); iter2++) 
				{
				  if ( (*iter).accurateMass() <= (*iter2).accurateMass() )
					  break;
				}
				glycan_composition_list.insert(iter2, *iter);
			}	*/
			return true ; 
		}  

  

		void GlycanIo::readGlycanCompositionLib(ifstream *in_f, vector<Engine::GlycanTheoretical::GlycanComposition> &gc_list) 
		{
			  char buffer[MAX_LINE];
			  int line_count, last_data_file_id;
			  int max_line_count = 100000;
			  int max_line_len = MAX_LINE;

			  gc_list.clear();

			  last_data_file_id = 0;
			  if ( in_f != NULL && in_f->is_open() ) 
			  {
				//read the headers
				vector<string> header_list;
				if ( ! in_f->eof() )
				{
				  memset(buffer, 0, max_line_len*sizeof(char));
				  in_f->getline (buffer,max_line_len);
				  Engine::Utilities::GStrTok *gtok1 = new Engine::Utilities::GStrTok(buffer, ","); 
				  for (int i=0; i< gtok1->size(); i++)
					  header_list.push_back(string(gtok1->at(i)));
				  delete(gtok1);
				}
				assert(header_list.size() > 0);

				line_count = 1;
				while (! in_f->eof() )
				{
				  memset(buffer, 0, MAX_LINE*sizeof(char));
				  in_f->getline (buffer,max_line_len);
				  line_count++;
				  int num = 0 ; 

				  if ( line_count >= max_line_count ) 
				  {
					cout << "total file size exceeds " << max_line_count << " lines." 
					 << "consider dividing the file"
					 << endl;
					assert( 0 );
				  }

				  if ( buffer[0] != '#' )
				  { //lines starts with # are comment lines
					  Engine::Utilities::GStrTok *gtok1 = new Engine::Utilities::GStrTok(buffer, ","); 
					  Engine::GlycanTheoretical::GlycanComposition gc = Engine::GlycanTheoretical::GlycanComposition();
					  gc.initNGlycanConstants() ; 

					  
					  
					  for (int i=0; i< gtok1->size(); i++) 
					  {
						  if ( i < (int) header_list.size() ) 
						  {
							  string &header = header_list[i];
							  if ( header.find("type", 0) != string::npos )
							  {
								  Engine::Utilities::StringUtils::strToLower(gtok1->at(i));
								  string type_value = string(gtok1->at(i)); 
								  if ( type_value.find("hybrid") != string::npos ) 
									  gc.glycan_type = Engine::GlycanTheoretical::HYBRID_GLYCAN;
								  else if ( type_value.find("complex") != string::npos ) 
									  gc.glycan_type = Engine::GlycanTheoretical::COMPLEX_GLYCAN;
								  else if ( type_value.find("highman") != string::npos ) 
									  gc.glycan_type = Engine::GlycanTheoretical::HIGH_MANNOSE_GLYCAN;
								  else 
								  {
									  cout << "unknown glycan type at line: " << line_count << endl;
									  assert(0);
								  }
							  }
							  else if ( header.find("mass", 0) != string::npos ) 
							  {
								  //do nothing
							  }
							  else if ( (header.find("glcnac_number", 0) != string::npos ) || (header.find("NumHexNAC",0) != string::npos))
							  {
								  gc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_HEXNAC] = atoi(gtok1->at(i));
								  num= num+ atoi(gtok1->at(i));
							  }
							  else if (header.find("NumHex", 0)!= string::npos)
							  {
								 gc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_HEX] += atoi(gtok1->at(i));
								 num= num+ atoi(gtok1->at(i));
							  }
							  else if ( header.find("man_number", 0) != string::npos) 
							  {
								  gc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_HEX] += atoi(gtok1->at(i));
								  num= num+ atoi(gtok1->at(i));

							  }
							  else if ( header.find("gal_number", 0) != string::npos ) 
							  {
								  gc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_HEX] += atoi(gtok1->at(i));
								  num= num+ atoi(gtok1->at(i));

							  }
							  else if ( (header.find("fuc_number", 0) != string::npos ) || (header.find("NumDeHex",0) != string::npos))
							  {
								  gc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_DEHEX] = atoi(gtok1->at(i));
								  num= num+ atoi(gtok1->at(i));

							  }
							  else if (( header.find("neunac_number", 0) != string::npos ) || (header.find("NumSia", 0) != string::npos) ) 
							  {
								  gc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_NEUAC] = atoi(gtok1->at(i));
								  num= num+ atoi(gtok1->at(i));
							  }
							  else if ( header.find("neugc_number", 0) != string::npos ) 
							  {
								  gc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_NEUGC] = atoi(gtok1->at(i));
								  num= num+ atoi(gtok1->at(i));

							  }
							  else if (header.find("so3_number", 0) != string::npos)
							  {
								  gc.consider_sulfated = true; 
								  gc.monosaccharide_count[Engine::GlycanTheoretical::MONOSACCHARIDE_MOD_SULFATE] = atoi(gtok1->at(i)) ; 
							  }
							  else 
							  {
								  cout << "unknown field: " << header << endl;
								  assert(0);
							  }
						  }
					  }
					delete(gtok1);					
					if (num > 0)
					{
						gc_list.push_back(gc);						
					}
				  }
				}
			  }
			  else 
			  {
				cout << "unable to read glycan library file" << endl;
				assert(0);
			  }
		}			
	}
}