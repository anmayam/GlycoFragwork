////---------------------------------------------------------------------------
//
//#ifndef RangeLocationFilterH
//#define RangeLocationFilterH
//
//#include "Sequence.h"
//#include "RangeLocation.h"
//#include <vector>
//
//using namespace std;
//
//namespace Engine
//{
//	namespace RangeLocation
//	{
//		class IRangeLocationFilter
//		{
//			public:
//			  IRangeLocationFilter(){};
//
//			  virtual ~IRangeLocationFilter(){};
//
//			  virtual void setSequence(Sequence* sequence) = 0;
//
//			  virtual bool accept(RangeLocation& rl) = 0;
//		};
//
//		class NGlycanFilter : public IRangeLocationFilter 
//		{
//			public:
//			  NGlycanFilter() 
//			  {
//				this->isNglycan = NULL;
//			  }
//
//			  virtual ~NGlycanFilter()
//			  {
//				delete isNglycan;
//			  }
//
//			  virtual void setSequence(Sequence* sequence);
//
//			  virtual bool accept(RangeLocation& rl);
//			private:
//			  bool* isNglycan;
//		};
//
//		class PeptideWeightFilter : public IRangeLocationFilter
//		{
//			public:
//			  PeptideWeightFilter(double minWeight, double maxWeight)
//			  {
//				this->minWeight = minWeight;
//				this->maxWeight = maxWeight;
//			  }
//
//			  virtual ~PeptideWeightFilter()
//			  { }
//
//			  virtual void setSequence(Sequence* sequence) {
//				this->sequence = sequence->getSeqString();
//			  }
//
//			  virtual bool accept(RangeLocation& rl);
//			private:
//			  string sequence;
//			  double minWeight;
//			  double maxWeight;
//		};
//
//		class AndRangeLocationFilter : public IRangeLocationFilter
//		{
//			public:
//			  AndRangeLocationFilter()
//			  { }
//			  
//			  virtual ~AndRangeLocationFilter(){ }
//
//			  virtual void setSequence(Sequence* sequence) {
//				for(size_t i = 0;i < filters.size();i++){
//				  filters[i]->setSequence(sequence);
//				}
//			  }
//
//			  virtual bool accept(RangeLocation& rl) {
//				for(size_t i = 0;i < filters.size();i++){
//				  if(!filters[i]->accept(rl)) {
//					return false;
//				  }
//				}
//
//				return true;
//			  }
//
//			  void addFilter(IRangeLocationFilter* filter){
//				filters.push_back(filter);
//			  }
//			private:
//			  vector<IRangeLocationFilter*> filters;
//		};
//
//	}
//}
////---------------------------------------------------------------------------
//#endif
