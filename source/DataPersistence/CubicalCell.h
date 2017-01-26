/// CubicalCell.h
/// Shaun Harker 2017-01-1944
/// MIT LICENSE

#ifndef CUBICALCELL_H
#define CUBICALCELL_H

/// CubicalCell
/// TODO: update to reflect higher dimensions
///   Represents a cell in a cubical complex
///   A cell is represented as a triple (x,y,type),
///      x : x-coordinate 0 <= x <= N
///      y : y-coordinate 0 <= x <= M
///   type : either 0,1,2,3 which when written as a two-bit binary 
///          number gives the shape of the cell:
///          As in the following diagram:
///           |
///          (10)    (11)
///           |               (types of 0, 1, 2, 3 written in binary as 00, 01, 10, 11)
///           |
///          (00) ---(01)----
///         Type 0 is vertices, Type 1 is Horizontal Edges, Type 2 is Vertical Edges, Type 3 is 2-cells
///   The positioning of a vertex at position (x,y) is understood to be at the lower-left corner of the (x,y) pixel
///   The positioning of a horizontal edge at position (x,y) is understood to be the bottom side of the (x,y) pixel
///   The positioning of a vertical edge at position (x,y) is understood to be the left side of the (x,y) pixel
///   Due to these conventions, for an (N,M) image we note the following:
///     2-cells will have (x,y) values satisfying            0 <= x < N,  0 <= y < M
///     Horizontal 1-cells will have (x,y) values satisfying 0 <= x < N,  0 <= y <= M
///     Vertical 1-cells will have (x,y) values satisfying   0 <= x <= N, 0 <= y < M
///     Vertices (0-cells) will have (x,y) values satisfying 0 <= x <= N, 0 <= y <= M

/// CubicalCell_
///   Underlying data for CubicalCell class
struct CubicalCell_ {
  std::vector<uint64_t> coordinates_;
  uint64_t type_;
};

/// CubicalCell 
class CubicalCell {
public:
  /// CubicalCell
  ///   Default constructor
  CubicalCell ( void ) {}

  /// data
  ///   mutator accessor
  CubicalCell_ &
  data ( void ) {
    if ( data_ . use_count() > 1 ) {
      // Must make a copy before in-place modification
      // or will invalidate other references
      auto newdata_ = std::make_shared<CubicalCell_>();
      *newdata_ = *data_;
      data_ = newdata_;
    } 
    // assert(data_.use_count() == 1);
    return *data_;
  }

  /// clone
  ///   Make a copy of the cell
  CubicalCell
  clone ( void ) const {
    return CubicalCell ( coordinates(), type() );
  }

  /// CubicalCell
  ///   Construct a cell with given position and type
  CubicalCell ( std::vector<uint64_t> const& coordinates, uint64_t type ) {
    assign(coordinates,type);
  }

  /// assign
  ///   Construct a cell with given position and type
  void
  assign ( std::vector<uint64_t> const& coordinates, uint64_t type ) {
    data_ = std::make_shared<CubicalCell_>();
    data_ -> coordinates_ = coordinates;
    data_ -> type_ = type;
  }

  /// coordinates
  ///  Return coordinates
  std::vector<uint64_t>
  coordinates ( void ) const {
    return data_ -> coordinates_;
  }

  /// type
  ///   Return type  ( type & (1 << d) iff cell has extent in variable d )
  uint64_t
  type ( void ) const {
    return data_ -> type_;
  }

  /// dimension
  ///   Return dimension of cell
  uint64_t
  dimension ( void ) const {
    uint64_t v = 0; 
    uint64_t x = type();
    while(x != 0) { x &= x - 1; v++; } // http://lemire.me/blog/2016/05/23/the-surprising-cleverness-of-modern-compilers/
    return v;
  }

  /// has_extent
  ///   Return true if cell is a 1-cell with width in dth-dimension
  ///   or if it is a two-cell
  bool 
  has_extent ( uint64_t d ) const {
    return type() & (1LL << d);
  }

  /// left
  ///   Return cell at left (face or coface)
  CubicalCell
  left ( uint64_t d ) const {
    CubicalCell result = clone();
    if ( not has_extent(d) ) result.data_->coordinates_[d] -= 1;
    result.data_->type_ ^= 1 << d;
    return result;
  }

  /// right
  ///   Return cell at right (face or coface)
  CubicalCell
  right ( uint64_t d ) const {
    CubicalCell result = clone();
    if ( has_extent(d) ) result.data_->coordinates_[d] += 1;
    result.data_->type_ ^= 1 << d;
    return result;
  }

  /// operator ==
  ///   Returns true if (coordinates, type) of two cells are all the same
  bool
  operator == ( CubicalCell const& rhs ) const {
    if ( type() != rhs.type() ) return false;
    if ( coordinates() != rhs.coordinates() ) return false;
    return true;
  }

  /// operator <
  ///   Lexicographical comparison on (type, coordinates).
  bool
  operator < ( CubicalCell const& rhs ) const {
    if ( type() == rhs.type() ) {
      if ( coordinates() == rhs.coordinates() ) {
        return false;
      } else {
        return coordinates() < rhs.coordinates();
      }
    } else {
      return type() < rhs.type ();
    }
  }

  /// operator <<
  friend std::ostream & operator << ( std::ostream & stream, CubicalCell const& stream_me ) {
    stream << "CubicalCell([";
    for ( auto coordinate : stream_me.coordinates() ) stream << coordinate << ",";
    return stream <<  "]," << (uint64_t) stream_me.type() << ")";
  }
private:
  std::shared_ptr<CubicalCell_> data_;
};

/// std::hash<CubicalCell>
namespace std {
  template<> struct hash<CubicalCell> {
    typedef CubicalCell argument_type;
    typedef std::size_t result_type;
    result_type operator()(argument_type const& cell) const {
      using boost::hash_value;
      using boost::hash_combine;
      std::size_t seed = 0;
      for ( auto coordinate : cell.coordinates() ) {
        hash_combine(seed,hash_value(coordinate));
      }      
      hash_combine(seed,hash_value(cell.type()));
      return seed;
    }
  };
}

#endif
