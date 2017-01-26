/// CubicalComplex.h
/// Shaun Harker 2016-12-15-2159
/// MIT LICENSE

#ifndef CUBICAL_COMPLEX_H
#define CUBICAL_COMPLEX_H

#include "CubicalCell.h"
#include "SmartConstIterator.h"

/// CubicalComplex
///   Implements a full cubical complex with Z_2 coefficients
///   Methods:
///     ImageComplex : Initialize the complex with width N and height M
///     boundary : given a cell, returns an array of boundary cells
///     coboundary : given a cell, returns an array of coboundary cells
///     cells : return the set of cells in complex

/// CubicalComplex_
///   Underlying data for CubicalComplex class
struct CubicalComplex_ {
  uint64_t size_;
  std::vector<uint64_t> sizes_;
};

/// CubicalComplex
class CubicalComplex {
public:
  /// CubicalComplex
  ///   Default constructor
  CubicalComplex ( void ) {}

  /// data
  ///   mutator accessor
  CubicalComplex_ &
  data ( void ) {
    if ( data_ . use_count() > 1 ) {
      // Must make a copy before in-place modification
      // or will invalidate other references
      auto newdata_ = std::make_shared<CubicalComplex_>();
      *newdata_ = *data_;
      data_ = newdata_;
    } 
    // assert(data_.use_count() == 1);
    return *data_;
  }

  /// clone
  ///   Make a copy of the complex
  CubicalComplex
  clone ( void ) const {
    return CubicalComplex ( sizes() );
  }

  /// CubicalComplex
  ///   Initialize the complex that is sizes[i] boxes across 
  ///   in dimensions d = 0, 1, ..., sizes.size() - 1
  CubicalComplex ( std::vector<uint64_t> const& sizes ) {
    assign ( sizes );
  }

  /// assign
  ///   Initialize the complex that is sizes[i] boxes across 
  ///   in dimensions d = 0, 1, ..., sizes.size() - 1
  void
  assign ( std::vector<uint64_t> const& sizes ) {
    data_ = std::make_shared<CubicalComplex_>();
    data_ -> size_ = 1;
    data_ -> sizes_ = sizes;
    for ( uint64_t d = 0; d < dimension(); ++ d ) {
      data_ -> size_ *= 2*sizes[d] + 1;
    }
  }
  
  /// sizes
  std::vector<uint64_t> const&
  sizes ( void ) const {
    return data_ -> sizes_;
  }

  /// Iterator
  typedef SmartConstIterator<CubicalCell> Iterator;

  /// begin
  Iterator
  begin ( void ) const {
    auto increment = [=](CubicalCell & cell){next(cell);};
    return Iterator(CubicalCell(std::vector<uint64_t>(dimension(), 0), 0), increment);
  }

  /// end
  Iterator
  end ( void ) const {
    auto increment = [=](CubicalCell & cell){next(cell);};
    return Iterator(CubicalCell(std::vector<uint64_t>(dimension(), 0), 1 << dimension()), increment);
  }

  /// size
  ///   Return number of cells in complex
  uint64_t
  size ( void ) const {
    return data_ -> size_;
  }

  /// boundary
  ///   Return array of boundary cells
  std::vector<CubicalCell>
  boundary ( CubicalCell const& cell ) const {
    std::vector<CubicalCell> bd;
    uint64_t type = cell.type();
    for ( uint64_t d = 0; d < dimension(); ++ d ) {
      if ( cell.has_extent(d) ) {
        CubicalCell left = cell.left(d);
        CubicalCell right = cell.right(d);
        if ( valid(left) ) bd.push_back(left);
        if ( valid(right) ) bd.push_back(right);     
      }
    }
    return bd;
  }

  /// coboundary
  ///   Return array of coboundary cells
  std::vector<CubicalCell>
  coboundary ( CubicalCell const& cell ) const {
    std::vector<CubicalCell> cbd;
    uint64_t type = cell.type();
    for ( uint64_t d = 0; d < dimension(); ++ d ) {
      if ( not cell.has_extent(d) ) {
        CubicalCell left = cell.left(d);
        CubicalCell right = cell.right(d);
        if ( valid(left) ) cbd.push_back(left);
        if ( valid(right) ) cbd.push_back(right);     
      }
    }
    return cbd;
  }    

  /// dimension
  ///   Return dimension of complex
  uint64_t
  dimension ( void ) const {
    return sizes().size();
  }

  /// next
  ///   Given a reference to a cell, make it the "next" cell 
  ///   in an iteration pattern
  void
  next ( CubicalCell & cell ) const {
    CubicalCell_ & cell_data = cell.data();
    std::vector<uint64_t> & coordinates = cell_data.coordinates_;
    uint64_t & type = cell_data.type_;
    if ( type == 1 << dimension() ) return; // next is idempotent on end.
    do {
      for ( uint64_t d = 0; d <= dimension(); ++ d ) {
        if ( d == dimension() ) {
          type += 1;
          if ( type == 1 << dimension() ) return;
        } else {
          coordinates[d] += 1;
          if ( coordinates[d] == sizes()[d] ) {
            coordinates[d] = 0;
          } else { 
            break;
          }
        }
      }
    } while ( not valid ( cell ) ); // TODO (can this O(d) check be removed?)
  }

  /// valid
  ///   Return true if cell is in complex with given sizes
  ///   (This is used by "next" to get the correct iteration
  ///    pattern)
  bool
  valid ( CubicalCell const& cell ) const {
    for ( uint64_t d = 0; d <= dimension(); ++ d ) {
      if ( cell.coordinates()[d] > sizes()[d] ) return false;
      if ( cell.coordinates()[d] == sizes()[d] && cell.has_extent(d) ) return false; 
    }
    return true;
  }

  /// operator ==
  bool
  operator == ( CubicalComplex const& rhs ) const {
    if ( sizes() != rhs.sizes() ) return false;
    return true;
  }

  /// operator <
  bool
  operator < ( CubicalComplex const& rhs ) const {
    return std::lexicographical_compare(sizes().begin(), sizes().end(), rhs.sizes().begin(), rhs.sizes().end());
  }

  /// operator <<
  friend std::ostream & operator << ( std::ostream & stream, CubicalComplex const& stream_me ) {
    stream << "CubicalComplex([";
    for ( auto x : stream_me.sizes() ) stream << x << ",";
    return stream <<  "])";
  }

private:
  std::shared_ptr<CubicalComplex_> data_;
};

/// std::hash<CubicalComplex>
namespace std {
  template<> struct hash<CubicalComplex> {
    typedef CubicalComplex argument_type;
    typedef std::size_t result_type;
    result_type operator()(argument_type const& complex) const {
      using boost::hash_value;
      using boost::hash_combine;
      std::size_t seed = 0;
      for ( auto x : complex.sizes() ) {
        hash_combine(seed,hash_value(x));
      }      
      return seed;
    }
  };
}

#endif
