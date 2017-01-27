/// CubicalComplex.h
/// Shaun Harker 2016-12-15-2159
/// MIT LICENSE

#ifndef CUBICAL_COMPLEX_H
#define CUBICAL_COMPLEX_H

#include "common.h"

/// CubicalComplex
///   Implements a trivial cubical complex with Z_2 coefficients
///   Methods:
///     ImageComplex : Initialize the complex with width N and height M
///     boundary : given a cell, returns an array of boundary cells
///     coboundary : given a cell, returns an array of coboundary cells
///     cells : return the set of cells in complex

/// CubicalComplex_
///   Underlying data for CubicalComplex class
struct CubicalComplex_ {
  // Defining
  std::vector<uint64_t> sizes_; // All other data follows from this.
  // Derived
  std::vector<uint64_t> place_values_;
  std::vector<uint64_t> shape_from_type_;
  std::vector<uint64_t> type_from_shape_;
  std::vector<uint64_t> begin_;
  std::vector<uint64_t> end_;
  uint64_t dimension_;
  uint64_t num_types_;
  uint64_t type_size_;
  uint64_t size_;
};

/// CubicalComplex
class CubicalComplex {
public:
  /// Types
  typedef uint64_t CellIndex;
  typedef boost::counting_iterator<CellIndex> Iterator;

  /// CubicalComplex
  ///   Default constructor
  CubicalComplex ( void ) {}

  /// CubicalComplex
  ///   Initialize the complex that is boxes[i] boxes across 
  ///   in dimensions d = 0, 1, ..., boxes.size() - 1
  ///   Note: The cubical complex does not have cells on the 
  ///         far right, so to have a "full" cubical 
  ///         complex as a subcomplex, pad with an extra box.
  CubicalComplex ( std::vector<uint64_t> const& sizes ) {
    assign ( sizes );
  }

  /// assign
  ///   Initialize the complex that is boxes[i] boxes across 
  ///   in dimensions d = 0, 1, ..., boxes.size() - 1
  void
  assign ( std::vector<uint64_t> const& boxes ) {
    data_ = std::make_shared<CubicalComplex_>();

    // Get dimension
    uint64_t D = boxes.size();

    // Compute PV := [1, boxes[0], boxes[0]*boxes[1], ..., boxes[0]*boxes[1]*...*boxes[D-1]]
    auto & PV = data_ -> place_values_;
    PV.resize ( D+1 );
    PV[0] = 1;
    std::partial_sum (boxes.begin(), boxes.end(), PV.begin() + 1, std::multiplies<uint64_t>()); 

    uint64_t M = 1L << D; // number of types
    uint64_t L = PV[D]; // Number of cells per shape/type
    uint64_t N = L * M; // cells per type * number of types

    // Set data
    data_ -> sizes_ = boxes;
    data_ -> dimension_ = D;
    data_ -> num_types_ = M;
    data_ -> type_size_ = L;
    data_ -> size_ = N;

    // Generate shapes and then sort them by dimension (which is bit popcount) to create types.
    // Implement a bijection between shapes and types via the following arrays:
    auto & ST = data_ -> shape_from_type_;
    auto & TS = data_ -> type_from_shape_;

    ST.resize ( M );
    TS.resize ( M );
    std::iota(ST.begin(), ST.end(), 0 );
    auto compare = [&](uint64_t x, uint64_t y) { return popcount_(x) < popcount_(y); };
    std::stable_sort(ST.begin(), ST.end(), compare);
    for ( uint64_t type = 0; type < M; ++ type) TS[ST[type]] = type; 

    // Set up iterator bounds for every dimension
    auto & begin = data_ -> begin_;
    auto & end = data_ -> end_;
    begin . resize ( dimension(), N );
    end . resize ( dimension(), 0 );
    for ( uint64_t type = 0, idx = 0; type < M; ++ type, idx += L ) {
      uint64_t shape = ST[type];
      uint64_t dim = popcount_(shape);
      begin[dim] = std::min(begin[dim], idx);
      end[dim] = std::max(end[dim], idx + L);
    }
  }
  
  /// sizes
  std::vector<uint64_t> const&
  sizes ( void ) const {
    return data_ -> sizes_;
  }

  /// begin
  Iterator
  begin ( void ) const {
    return Iterator(0);
  }

  /// end
  Iterator
  end ( void ) const {
    return Iterator(size());
  }

  /// size
  ///   Return number of cells in complex
  uint64_t
  size ( void ) const {
    return data_ -> size_;
  }

  /// begin
  ///   begin iterator for dimension dim
  Iterator
  begin ( uint64_t dim ) const {
    return Iterator(data_ -> begin_[dim]);
  }

  /// end
  ///   end iterator for dimension dim
  Iterator
  end ( uint64_t dim ) const {
    return Iterator(data_ -> end_[dim]);
  }

  /// size
  ///   Return number of dim dimensional cells in complex
  uint64_t
  size ( uint64_t dim ) const {
    return (uint64_t) (end(dim) - begin(dim));
  }

  /// shape_begin
  uint64_t
  shape_begin ( uint64_t shape ) const {
    return TS_()[shape] * type_size_();
  }

  /// shape_end
  uint64_t
  shape_end ( uint64_t shape ) const {
    return TS_()[shape] * ( type_size_() + 1);
  }

  /// operator ()
  ///   Return an iterator pair for traversal 
  ///   of dim dimensional cells
  std::pair<Iterator, Iterator>
  operator () ( uint64_t dim ) const {
    return { begin(dim), end(dim) };
  }

  /// coordinates
  ///   Given a cell index, 
  ///   returns ( x_0, x_1, ..., x_{dim-1} )
  std::vector<uint64_t>
  coordinates ( CellIndex cell ) const {
    std::vector<uint64_t> result ( dimension() );
    for ( uint64_t d = 0; d < dimension(); ++ d ) {
      result[d] = cell % sizes()[d];
      cell /= sizes()[d];
    }
    return result;
  }

  /// cell_shape
  ///   Give shape code
  ///   Interpretation: if ( shape & ( 1 << i ) ) { then the cell has extent in dimension i }
  uint64_t
  cell_shape ( CellIndex cell ) const {
    uint64_t type = cell / type_size_();
    uint64_t shape = ST_() [ type ]; 
    return shape;
  }

  /// boundary
  ///   Return array of boundary cells
  std::vector<CellIndex>
  boundary ( CellIndex cell ) const {
    std::vector<CellIndex> bd;
    uint64_t shape = cell_shape(cell);
    uint64_t type = TS_() [ shape ];
    lldiv_t coordinate = {(int64_t)cell, 0}; // (quotient, remainder), see std::div
    for ( uint64_t d = 0, bit = 1; d < dimension(); ++ d, bit <<= 1L ) {
      // If cell has no extent in this dimension, no boundaries.
      if ( not (shape & bit) ) continue;
      uint64_t offset_cell = cell + type_size_() * ( TS_() [ shape ^ bit ] - type );
      // Otherwise, the cell does have extent in this dimension.
      // It is always the case that such a cell has a boundary to the left.
      bd.push_back( offset_cell );
      // Check if there is a boundary to the right:
      coordinate = std::div(coordinate.quot, sizes()[d] ); 
      if ( coordinate.rem + 1 < sizes()[d]) { 
        bd.push_back( offset_cell + PV_()[d]);
      }
    }
    return bd;
  }

  /// coboundary
  ///   Return array of coboundary cells
  std::vector<CellIndex>
  coboundary ( CellIndex cell ) const {
    std::vector<CellIndex> cbd;
    uint64_t shape = cell_shape(cell);
    uint64_t type = TS_() [ shape ];
    lldiv_t coordinate = {(int64_t)cell, 0};
    for ( uint64_t d = 0, bit = 1; d < dimension(); ++ d, bit <<= 1L ) {
      // If cell has extent in this dimension, no coboundaries.
      if ( shape & bit ) continue;
      uint64_t offset_cell = cell + type_size_() * ( TS_() [ shape ^ bit ] - type );
      // Otherwise, the cell does not have extent in this dimension.
      // It is always the case that such a cell has a coboundary to the right.
      cbd.push_back( offset_cell );
      // Check if there is a coboundary to the left:
      coordinate = std::div(coordinate.quot, sizes()[d] ); 
      if ( coordinate.rem > 0 ) { 
        cbd.push_back( offset_cell - PV_()[d]);
      }
    }
    return cbd;
  }    

  /// dimension
  ///   Return dimension of complex
  uint64_t
  dimension ( void ) const {
    return data_ -> dimension_;
  }

  /// dimension
  ///   Return dimension of cell
  uint64_t
  dimension ( CellIndex cell ) const {
    return popcount_(cell_shape(cell));
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

  std::vector<uint64_t> const&
  TS_ ( void ) const {
    return data_ -> type_from_shape_;
  }

  std::vector<uint64_t> const&
  ST_ ( void ) const {
    return data_ -> shape_from_type_;
  }

  std::vector<uint64_t> const&
  PV_ ( void ) const {
    return data_ -> place_values_;
  }

  uint64_t
  type_size_ ( void ) const {
    return data_ -> type_size_;
  }

  uint64_t
  popcount_ ( uint64_t x ) const {
    // http://lemire.me/blog/2016/05/23/the-surprising-cleverness-of-modern-compilers/
    uint64_t pcnt = 0; 
    while(x != 0) { x &= x - 1; ++pcnt; } 
    return pcnt;
  }
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
