/// MixedRadix.h
/// MIT LICENSE 2016 Shaun Harker
/// 2016-12-15-2305

typedef std::vector<uint64_t> MixedRadixNumber;

class MixedRadixSystem{
public:
  MixedRadixNumber add ( MixedRadixNumber const& lhs, MixedRadixNumber const& rhs ) const;
  uint64_t digit ( MixedRadixNumber const& n, uint64_t i ) const;
  
  MixedRadixNumber 
  make ( uint64_t n ) const;
  
  uint64_t 
  unmake ( MixedRadixNumber const& mrn ) const;
private:
  std::shared_pointer<MixedRadixNumber> base_;
};
