// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// Written by Anoop Mayampurath and Yin Wu, Indiana University

namespace Engine
{
	namespace GlycanCompositionManager
	{
		class GlycanComposition
		{
		public:
			//! default constructor
			GlycanComposition() ; 
			//! destructor.
			~GlycanComposition() ;
			//! search ions for peptides
			void SearchIonList() ; 
		};
	}
}