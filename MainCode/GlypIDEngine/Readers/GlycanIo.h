#include <string>
#include <fstream>
#include <vector>
#include <map>
#include "../GlycanTheoretical/Nglycan.h"

using std::string;
using std::ifstream;
using std::vector;

namespace Engine
{
	namespace Readers
	{
		class GlycanIo
		{
			
			std::vector<Engine::GlycanTheoretical::GlycanComposition> glycan_compositions_in_search;

			Engine::GlycanTheoretical::GlycanComposition max_glycan_composition, min_glycan_composition;


			public:
				GlycanIo(void); 
				~GlycanIo(void) ; 			

				bool LoadGlycanListFromFile(char *glycanFileName, vector<Engine::GlycanTheoretical::GlycanComposition> &glycan_composition_list);															
				void readGlycanCompositionLib(ifstream *in_f, vector<Engine::GlycanTheoretical::GlycanComposition> &gc_list)  ; 
				
		};
	}
}



