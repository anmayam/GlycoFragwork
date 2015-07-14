//the pixel data structure used for 
//glycosylation applications


#include <algorithm>
#include <vector>


namespace Engine
{
	namespace PeakProcessing
	{
		class GPixel 
		{
		
		public:
			  float value;
			  float x, y, z; //the coordinates of the pixel

			  GPixel() 
			  {
				  value = 0;
				  x = y = z = 0;
			  };
			  
			  GPixel(float v, float x, float y, float z)
			  {
				  value = v;
				  this->x = x;
				  this->y = y;
				  this->z = z;
			  };
		};
	}
}