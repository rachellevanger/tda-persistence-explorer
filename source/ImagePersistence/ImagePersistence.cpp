// ImagePersistence.cpp
//  Compute persistence homology of 2D images using PHAT.
//  MIT LICENSE 2016 Jacek Cyranka
// 
// Revision History:
//   2016-11-21 Shaun Harker 
//      * Refactored for readability
//      * fixed non-square image bug 
//      * Improved make system

#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <queue>
#include <boost/functional/hash.hpp>
#include "phat/representations/vector_vector.h"
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/twist_reduction.h>
#include <phat/compute_persistence_pairs.h>
#include "CImg.h"

using namespace cimg_library;

typedef std::function<uint8_t(uint64_t, uint64_t)> Image;

/// Cell
///   Represents a 0, 1, or 2D cell in a 2D cubical complex
///   Implementation:
///      x,y give pixel coordinates
///      type of 0,1,2,3 gives dimension and shape of cell
///      vertices are understood to be at the bottom left of (x,y) pixel
///      vertical edges are understood to be at bottom of (x,y) pixel
///      horizontal edges are understood to be at left of (x,y) pixel
///    As in the following diagram:
///       |
///      (10)    (11)
///       |               (types of 0, 1, 2, 3 written in binary as 00, 01, 10, 11)
///       |
///      (00) ---(01)----
class Cell {
public:
  /// Cell
  Cell ( uint64_t x, uint64_t y, uint8_t type ) : x_(x), y_(y), type_(type) {}
  /// dimension
  ///   Return dimension of cell
  uint64_t
  dimension ( void ) const {
    uint64_t result;
    switch ( type_ ) {
      case 0: result = 0; break; // vertex
      case 1: result = 1; break; // horizontal edge
      case 2: result = 1; break; // vertical edge
      case 3: result = 2; break; // 2-cell
      default: throw std::runtime_error("Invalid cell type");
    }
    return result;
  }
  /// x
  ///  Return x coordinate
  uint64_t
  x ( void ) const {
    return x_;
  }
  /// y
  ///  Return y coordinate
  uint64_t
  y ( void ) const {
    return y_;
  }
  /// type
  ///   Return type  ( 0 = vertex, 1 = horizontal 1-cell, 2 = vertical 1-cell, 3 = 2-cell )
  uint8_t
  type ( void ) const {
    return type_;
  }
  /// has_horizontal_extent
  ///   Return true if cell is a 1-cell with width in x-dimension
  ///   or if it is a two-cell
  bool 
  has_horizontal_extent (void) const {
    return type_ & 1;
  }
  /// has vertical_extent
  ///   Return true if cell is a 1-cell with height in y-dimension
  ///   or if it is a two-cell
  bool 
  has_vertical_extent (void) const {
    return type_ & 2;
  }
  /// operator ==
  bool
  operator == ( Cell const& rhs ) const {
    if ( x() != rhs.x() ) return false;
    if ( y() != rhs.y() ) return false;
    if ( type() != rhs.type() ) return false;
    return true;
  }
  /// operator <
  ///   Sorts by dimension, sorting horizontal before vertical
  ///   Then sorts by x, then by y
  bool
  operator < ( Cell const& rhs ) const {
    if ( type() == rhs.type() ) {
      if ( x() == rhs.x() ) {
        return y() < rhs.y();
      } else {
        return x() < rhs.x();
      }
    } else {
      return type() < rhs.type ();
    }
  }
  /// operator <<
  friend std::ostream & operator << ( std::ostream & stream, Cell const& stream_me ) {
    return stream << "(" << stream_me.x() << ", " << stream_me.y() << ", " << (int) stream_me.type() << ")";
  }
private:
  uint64_t x_;
  uint64_t y_;
  uint8_t type_; // 0 for lower-left vertex, 1 for lower horizontal edge
};


struct CellHasher {
  std::size_t operator()(const Cell& cell) const {
      using boost::hash_value;
      using boost::hash_combine;
      std::size_t seed = 0;
      hash_combine(seed,hash_value(cell.x()));
      hash_combine(seed,hash_value(cell.y()));
      hash_combine(seed,hash_value(cell.type()));
      return seed;
  }
};

/// ImageComplex
///   Implements a 2D full cubical complex that is N wide and M tall.
class ImageComplex {
public:
  /// ImageComplex
  ImageComplex ( uint64_t N, uint64_t M ) : N(N), M(M) {}
  /// boundary
  ///   Return array of boundary cells
  std::vector<Cell>
  boundary ( Cell const& cell ) const {
    std::vector<Cell> bd;
    uint64_t x = cell.x();
    uint64_t y = cell.y();
    uint64_t type = cell.type();
    if ( cell.has_horizontal_extent() ) {
      Cell left_bd = Cell(x,y,type & 2);
      Cell right_bd = Cell(x+1,y, type & 2);
      if ( _valid(left_bd) ) bd.push_back(left_bd);
      if ( _valid(right_bd) ) bd.push_back(right_bd);
    }
    if ( cell.has_vertical_extent() ) {
      Cell below_bd = Cell(x,y,type & 1);
      Cell above_bd = Cell(x,y+1, type & 1);
      if ( _valid(below_bd) ) bd.push_back(below_bd);
      if ( _valid(above_bd) ) bd.push_back(above_bd);
    }
    return bd;
  }
  /// coboundary
  ///   Return array of coboundary cells
  std::vector<Cell>
  coboundary ( Cell const& cell ) const {
    std::vector<Cell> cbd;
    uint64_t x = cell.x();
    uint64_t y = cell.y();
    uint64_t type = cell.type();
    if ( not cell.has_horizontal_extent() ) {
      Cell left_cbd = Cell(x-1,y,type | 1);
      Cell right_cbd = Cell(x,y, type | 1);
      if ( _valid(left_cbd) ) cbd.push_back(left_cbd);
      if ( _valid(right_cbd) ) cbd.push_back(right_cbd);
    }
    if ( not cell.has_vertical_extent() ) {
      Cell below_cbd = Cell(x,y-1,type | 2);
      Cell above_cbd = Cell(x,y, type | 2);
      if ( _valid(below_cbd) ) cbd.push_back(below_cbd);
      if ( _valid(above_cbd) ) cbd.push_back(above_cbd);
    }
    return cbd;
  }     
private:
  uint64_t N;
  uint64_t M;
  /// valid
  ///   Return true if cell is in 2D N x M complex
  bool
  _valid ( Cell const& cell ) const {
    //if ( cell.x() < 0 || cell.y() < 0 ) return false;
    if ( cell.x() > N || cell.y() > M ) return false;
    if ( cell.x() == N && cell.has_horizontal_extent() ) return false;
    if ( cell.y() == M && cell.has_vertical_extent() ) return false;
    return true;
  }
};

void
ComputePersistencePairsViaPHAT ( ImageComplex const& complex,
                                 std::vector<std::pair<Cell, int64_t>> const& cell_filtration,
                                 std::string const& outfile_name ) {
  // Index the cells
  uint64_t num_cells = cell_filtration.size();
  std::cout << "Number of cells = " << num_cells << "\n";
  std::unordered_map<Cell, uint64_t, CellHasher> cell_indexing;
  for ( uint64_t i = 0; i < num_cells; ++ i ) {
    cell_indexing[cell_filtration[i].first] = i;
  }

  // Create boundary matrix object
  phat::boundary_matrix< phat::vector_vector > boundary_matrix;

  // Set the number of columns (i.e. number of cells)
  boundary_matrix . set_num_cols(num_cells);

  // Set column for each cell
  for ( phat::index i = 0; i < num_cells; ++ i ) {
    Cell const& cell = cell_filtration[i].first;
    std::vector<phat::index> boundary;
    for ( Cell const& bd_cell : complex.boundary(cell) ) {
      boundary.push_back(cell_indexing[bd_cell]);
    }
    std::sort(boundary.begin(), boundary.end()); // is this required?
    boundary_matrix . set_dim( i, cell.dimension());    
    boundary_matrix . set_col( i, boundary );
  }
  // Call PHAT
  phat::persistence_pairs pairs;  
  phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );

  // Output result
  std::ofstream outfile ( outfile_name );
  outfile << "dim,birth,b_x,b_y,b_z,death,d_x,d_y,d_z\n" ;
  for( int i = 0; i < pairs.get_num_pairs(); ++ i ){
    auto birth_cell_index = pairs.get_pair(i).first;
    auto death_cell_index = pairs.get_pair(i).second;
    auto birth_cell = cell_filtration[birth_cell_index].first;
    auto birth_value = cell_filtration[birth_cell_index].second;
    auto death_cell = cell_filtration[death_cell_index].first;
    auto death_value = cell_filtration[death_cell_index].second;
    if ( birth_value == death_value ) continue;
    outfile << birth_cell.dimension() << ", "
            << birth_value << ", "
            << birth_cell.x() << ", "
            << birth_cell.y() << ", "
            << "0" << ", "
            << death_value << ", "
            << death_cell.x() << ", "
            << death_cell.y() << ", "
            << "0" << "\n";
  }
}

void
sublevel_persistence ( uint64_t N, uint64_t M, Image const& pixel_data, std::string const& outfile_name ) {
  // Create N x M full 2D cubical complex
  ImageComplex complex(N,M);
  std::unordered_map<Cell, uint64_t, CellHasher> cell_values;
  std::queue<Cell> work_queue;

  // Create 2-cells
  for ( int x = 0; x < N; ++ x ) {
    for ( int y = 0; y < M; ++ y ) {
      Cell cell = Cell(x,y,3);
      cell_values[cell] = pixel_data(x,y);
      work_queue.push(cell);
    }
  }

  // Recursively determine boundary cells and assign values
  while ( not work_queue . empty () ) {
    Cell work_cell = work_queue . front (); 
    work_queue . pop ();
    for ( Cell const& bd_cell : complex.boundary(work_cell) ) {
      if ( cell_values.count(bd_cell) == 0 ) { 
        work_queue.push(bd_cell);
        cell_values[bd_cell] = cell_values[work_cell];
      } else {
        cell_values[bd_cell] = std::min(cell_values[bd_cell], cell_values[work_cell]);
      }
    }
  }

  // Assign an indexing to cells in complex
  std::vector<std::pair<Cell, int64_t>> cell_filtration ( cell_values.begin(), cell_values.end() );
  std::sort(cell_filtration.begin(), cell_filtration.end(), 
    [](std::pair<Cell, int64_t> const& a, std::pair<Cell, int64_t> const& b){
      if ( a.second == b.second ) return a.first < b.first;
      return a.second < b.second;
    });

  ComputePersistencePairsViaPHAT(complex, cell_filtration, outfile_name);
}

void
superlevel_persistence ( uint64_t N, uint64_t M, Image const& pixel_data, std::string const& outfile_name ) {
  // Create a (N-1) x (M-1) cubical complex
  ImageComplex complex(N-1,M-1);
  std::unordered_map<Cell, uint64_t, CellHasher> cell_values;
  std::queue<Cell> work_queue;

  // Create vertices
  for ( int x = 0; x <= N; ++ x ) {
    for ( int y = 0; y <= M; ++ y ) {
      Cell cell = Cell(x,y,0);
      cell_values[cell] = pixel_data(x,y);
      work_queue.push(cell);
    }
  }

  // Recursively determine coboundary cells and assign values
  while ( not work_queue . empty () ) {
    Cell work_cell = work_queue . front (); 
    work_queue . pop ();
    for ( Cell const& cbd_cell : complex.coboundary(work_cell) ) {
      if ( cell_values.count(cbd_cell) == 0 ) { 
        work_queue.push(cbd_cell);
        cell_values[cbd_cell] = cell_values[work_cell];
      } else {
        cell_values[cbd_cell] = std::min(cell_values[cbd_cell], cell_values[work_cell]);
      }
      // Note: if we had negated image, we would want max of boundary values.
      // As an implementation detail, we do not negate the image, so we use min here.
      // Below, we sort the values in reverse order to accommodate this.
    }
  }

  // Assign an indexing to cells in complex
  std::vector<std::pair<Cell, int64_t>> cell_filtration ( cell_values.begin(), cell_values.end() );
  std::sort(cell_filtration.begin(), cell_filtration.end(), 
    [](std::pair<Cell, int64_t> const& a, std::pair<Cell, int64_t> const& b){
      if ( a.second == b.second ) return a.first < b.first;
      return a.second > b.second; // reverse order
    });

  ComputePersistencePairsViaPHAT(complex, cell_filtration, outfile_name);
}

int main(int argc, char *argv[]) {
  std::string infile_name(argv[1]);
  std::string outfile_name(argv[2]);
  std::string mode(argv[3]);
  
  // Check arguments
  if ( mode != "sub" && mode != "super" ) {
    std::cerr << "Usage: \n  ImagePersistence <image_filename> <output_filename> <mode> \n mode is either 'sub' or 'super'\n";
    return 1;
  }

  // Load image file with CImg
  CImg<unsigned char> image( argv[1] );

  // Persistence calculation
  uint64_t N = image . width ();
  uint64_t M = image . height ();
  auto pixel_data = [&](int x, int y){return image(x,y,0,1);};
  if ( mode == "sub" ) {
    sublevel_persistence(N,M,pixel_data,outfile_name);
  } 
  if ( mode == "super" ) {
    superlevel_persistence(N,M,pixel_data,outfile_name);
  }
  return 0;
}
