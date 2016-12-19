/// CubicalComplex.h
/// Shaun Harker 2016-12-15-2159
/// MIT LICENSE

#ifndef CUBICAL_COMPLEX_H
#define CUBICAL_COMPLEX_H

/// Cell
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
class Cell {
public:
  /// Cell
  ///   Default constructor
  Cell ( void ) {}
  /// Cell
  ///   Construct a cell with given position and type
  Cell ( std::vector<uint64_t> const& coordinates, uint64_t type ) : coordinates_(coordinates), type_(type) {}
  /// dimension
  ///   Return dimension of cell
  uint64_t
  dimension ( void ) const {
    uint64_t v = 0; 
    uint64_t x = type_;
    while(x != 0) { x &= x - 1; v++; } // http://lemire.me/blog/2016/05/23/the-surprising-cleverness-of-modern-compilers/
    return v;
  }
  /// coordinates
  ///  Return coordinates
  std::vector<uint64_t>
  coordinates ( void ) const {
    return coordinates_;
  }
  /// type
  ///   Return type  ( type & (1 << d) iff cell has extent in variable d )
  uint64_t
  type ( void ) const {
    return type_;
  }
  /// has_extent
  ///   Return true if cell is a 1-cell with width in dth-dimension
  ///   or if it is a two-cell
  bool 
  has_extent ( uint64_t d ) const {
    return type_ & (1LL << d);
  }
  /// left
  ///   Move to cell at left (face or coface)
  Cell
  left ( uint64_t d ) const {
    Cell result = *this;
    if ( not has_extent(d) ) result.coordinates_[d] -= 1;
    result.type_ ^= 1 << d;
    return result;
  }
  /// right
  ///   Move to cell at right (face or coface)
  right ( uint64_t d ) const {
    Cell result = *this;
    if ( has_extent(d) ) result.coordinates_[d] += 1;
    result.type_ ^= 1 << d;
    return result;
  }
  /// operator ==
  ///   Returns true if (coordinates, type) of two cells are all the same
  bool
  operator == ( Cell const& rhs ) const {
    if ( type() != rhs.type() ) return false;
    if ( coordinates() != rhs.coordinates() ) return false;
    return true;
  }
  /// operator <
  ///   Lexicographical comparison on (type, coordinates).
  bool
  operator < ( Cell const& rhs ) const {
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
  friend std::ostream & operator << ( std::ostream & stream, Cell const& stream_me ) {
    stream << "([";
    for ( auto coordinate : coordinates() ) stream << stream_me.coordinates() << ",";
    return stream <<  "]," << (uint64_t) stream_me.type() << ")";
  }
private:
  std::vector<uint64_t> coordinates_;
  uint64_t type_;
};

/// CellHasher
///   Provides a hash function for Cell objects
///   This is required in order to use hash tables (e.g. std::unordered_map)
struct CellHasher {
  std::size_t operator()(const Cell& cell) const {
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

class CellIterator {
  CellIterator ( CubicalComplex complex, Cell cell ) : complex_(complex) {
    coordinates_ = cell.coordinates();
    type_ = cell.type();
  }
  /// operator *
  ///   Dereference iterator
  Cell
  operator * ( void ) {
    return Cell(coordinates_, type_);
  }
  /// operator ++
  ///   Increment to next cell in complex with widths "sizes[i]"
  ///   Iterates through cartesian product of coordinate possibilities
  ///   and types
  void
  operator ++ ( void ) {
    do {
      for ( uint64_t d = 0; d <= dimension(); ++ d ) {
        if ( d == dimension() ) {
          type_ += 1;
        } else {
          coordinates_[d] += 1;
          if ( coordinates_[d] == sizes[d] ) {
            coordinates_[d] = 0;
          } else { 
            break;
          }
        }
      }
    } while ( not valid ( sizes ) );
  }
  /// operator ==
  bool
  operator == ( CellIterator const& rhs ) const {
    return ( coordinates_ == rhs.coordinates_ && type_ == rhs.type_ );
  }
  /// operator !=
  bool
  operator != ( CellIterator const& rhs ) const {
    return not (*this == rhs);
  }
  /// valid
  ///   Return true if cell is in complex with given sizes
  bool
  valid ( void ) const {
    for ( uint64_t d = 0; d <= dimension(); ++ d ) {
      if ( coordinates_[d] > sizes[d] ) return false;
      if ( coordinates_[d] == sizes[d] && cell.has_extent(d) ) return false; 
    }
    return true;
  }
private:
  CubicalComplex complex_;
  std::vector<uint64_t> coordinates_;
  uint64_t type_;
};

/// CubicalComplex
///   Implements a full cubical complex with Z_2 coefficients
///   Methods:
///     ImageComplex : Initialize the complex with width N and height M
///     boundary : given a cell, returns an array of boundary cells
///     coboundary : given a cell, returns an array of coboundary cells
///     cells : return the set of cells in complex
class CubicalComplex {
public:
  /// CubicalComplex
  ///   Default constructor
  CubicalComplex ( void ) {}
  /// ImageComplex
  ///   Initialize the complex that is sizes[i] boxes across 
  ///   in dimensions d = 0, 1, ..., sizes.size() - 1
  CubicalComplex ( std::vector<uint64_t> const& sizes ) : sizes_(sizes) {
    size_ = 1;
    for ( uint64_t d = 0; d < dimension(); ++ d ) {
      size_ *= 2*sizes_[d] + 1;
    }
  }
  /// boundary
  ///   Return array of boundary cells
  std::vector<Cell>
  boundary ( Cell const& cell ) const {
    std::vector<Cell> bd;
    uint64_t type = cell.type();
    for ( uint64_t d = 0; d < dimension(); ++ d ) {
      if ( cell.has_extent(d) ) {
        Cell left = cell.left();
        Cell right = cell.right();
        if ( left.valid(sizes()) ) bd.push_back(left);
        if ( right.valid(sizes()) ) bd.push_back(right);     
      }
    }
    return bd;
  }
  /// coboundary
  ///   Return array of coboundary cells
  std::vector<Cell>
  coboundary ( Cell const& cell ) const {
    std::vector<Cell> cbd;
    uint64_t type = cell.type();
    for ( uint64_t d = 0; d < dimension(); ++ d ) {
      if ( not cell.has_extent(d) ) {
        Cell left = cell.left();
        Cell right = cell.right();
        if ( left.valid(sizes()) ) cbd.push_back(left);
        if ( right.valid(sizes()) ) cbd.push_back(right);     
      }
    }
    return cbd;
  }    

  /// begin
  CellIterator
  begin ( void ) const {
    return CellIterator(*this, Cell(std::vector<uint64_t>(dimension(), 0), 0));
  }

  /// end
  CellIterator
  end ( void ) const {
    return CellIterator(*this, Cell(std::vector<uint64_t>(dimension(), 0), 1 << dimension()));
  }

  /// size
  ///   Return number of cells in complex
  uint64_t
  size ( void ) const {
    return size_;
  }
  
  /// sizes
  std::vector<uint64_t> const&
  sizes ( void ) const {
    return sizes_;
  }

private:
  std::vector<uint64_t> sizes_;
};

#endif
