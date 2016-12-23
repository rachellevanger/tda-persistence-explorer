// ImagePersistence.cpp
//  Compute persistence homology of 2D images using PHAT.
//  MIT LICENSE 2016 Jacek Cyranka
// 
// Revision History:
//   2016-11-21 Shaun Harker 
//      * refactored
//      * fixed non-square image bug 
//      * improved make system

#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
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

std::string help_string = 
" Usage: ImagePersistence <image_filename> <output_filename> <mode> \n"
"   <image_filename>  : input image filename \n"
"   <output_filename> : output filename (a CSV file)\n"
"   <mode>            : A mode of execution (either \"sub\" or \"super\") \n"
"                       to indicate whether to compute\n"
"                       sublevel or superlevel persistence.\n";


// Overview:
//  1. The input filename is loaded into an "Image" object
//  2. Based on "mode", either a sublevel or superlevel filtration is constructed.
//  3. Persistence is calculated for the filtration
//  4. Results are saved to the requested output file.

// Implementation overview:
// Classes
//   * `Image`        : class for loading and accessing image data
//   * `Cell`         : data type for cell in 2D complex
//   * `ImageComplex` : implements boundary, coboundary for 2D complex
//   * `Filtration`   : class for storing filtrations of complexes
// Functions
//   * `SublevelFiltration`     : Given an Image object, construct a Filtration object
//   * `SuperlevelFiltration`   : Given an Image object, construct a Filtration object
//   * `PersistenceViaPHAT`     : Given a Filtration object, compute persistence pairs
//   * `SavePersistenceResults` : Given Filtration and PersistencePairs, save result to file                              
//   * `main`                   : Parse command line arguments, create filtration, 
//                                compute persistence, and write output to file

/// Image
///   This class is used to load and access 2D image data
class Image {
public:
  /// Image
  ///   Default constructor
  Image ( void ) {}
  /// Image
  ///   Load image given image_filename
  Image ( std::string const& image_filename ) : image_(image_filename.c_str()) {}
  /// width
  ///   Return width of image in pixels
  uint64_t
  width ( void ) const {
    return image_ . width ();
  }
  /// height
  ///   Return height of image in pixels
  uint64_t
  height ( void ) const {
    return image_ . height ();
  }
  /// data
  ///   Given 0 <= x < width() and 0 <= y < height(), 
  ///   return pixel data at position (x,y)
  uint64_t
  data ( uint64_t x, uint64_t y ) const {
    return image_(x,y,0,1);
  }
private:
  CImg<unsigned char> image_;
};

/// Cell
///   Represents a 0, 1, or 2D cell in a 2D cubical complex
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
  Cell ( uint64_t x, uint64_t y, uint64_t type ) : x_(x), y_(y), type_(type) {}
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
  uint64_t
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
  ///   Returns true if (x, y, type) of two cells are all the same
  bool
  operator == ( Cell const& rhs ) const {
    if ( x() != rhs.x() ) return false;
    if ( y() != rhs.y() ) return false;
    if ( type() != rhs.type() ) return false;
    return true;
  }
  /// operator <
  ///   Lexicographical comparison on (type, x, y).
  ///   This makes vertices < horizontal edges < vertical edges < 2-cells
  ///   and ties are broken by comparing x and y
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
  uint64_t type_; // 0 for lower-left vertex, 1 for lower horizontal edge
};

/// CellHasher
///   Provides a hash function for Cell objects
///   This is required in order to use hash tables (e.g. std::unordered_map)
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
///   Implements a 2D full cubical complex with Z_2 coefficients that is N wide and M tall.
///   Methods:
///     ImageComplex : Initialize the complex with width N and height M
///     boundary : given a cell, returns an array of boundary cells
///     coboundary : given a cell, returns an array of coboundary cells
///     cells : return the set of cells in complex
class ImageComplex {
public:
  /// ImageComplex
  ///   Default constructor
  ImageComplex ( void ) {}
  /// ImageComplex
  ///   Initialize the complex with width N and height M
  ImageComplex ( uint64_t N, uint64_t M ) : N_(N), M_(M) {}
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
  /// cells
  ///   Return set of cells in the complex
  std::unordered_set<Cell, CellHasher>
  cells ( void ) const {
    std::unordered_set<Cell, CellHasher> result;
    for ( uint64_t x = 0; x <= N_; ++ x ) {
      for ( uint64_t y = 0; y <= M_; ++ y ) {
        for ( uint64_t type = 0; type < 4; ++ type ) {
          Cell cell = Cell(x,y,type);
          if ( _valid(cell) ) result . insert ( cell );
        }
      }
    }
    return result;
  } 
  /// size
  ///   Return number of cells in complex
  uint64_t
  size ( void ) const {
    return (N_+1)*(M_+1) + N_*M_ + (N_+1)*M_ + N_*(M_+1); // == 4NM + 2N + 2M + 1
  }
private:
  uint64_t N_;
  uint64_t M_;
  /// valid
  ///   Return true if cell is in 2D N x M complex
  bool
  _valid ( Cell const& cell ) const {
    //if ( cell.x() < 0 || cell.y() < 0 ) return false;
    if ( cell.x() > N_ || cell.y() > M_ ) return false;
    if ( cell.x() == N_ && cell.has_horizontal_extent() ) return false;
    if ( cell.y() == M_ && cell.has_vertical_extent() ) return false;
    return true;
  }
};

/// Filtration
///   A filtration of a complex is a total ordering of its cells such 
///   that for any lower set of the cells furnishes a closed subcomplex.
///   This class implements a special case of filtrations arising from
///   ascending or descending sorts of values associated to cells.
///   The following methods are provided:
///     * complex     : access to underlying complex
///     * cell(i)     : give the ith cell in the total ordering
///     * index(cell) : given a cell in the complex, give its position in the ordering
///     * value(i)    : given a positionin the filtration, give the value of the associated cell
class Filtration {
public:
  /// Filtration
  ///   Default constructor
  Filtration ( void ) {}
  /// Filtration
  ///   Overview:
  ///     Construct a filtration of the given complex according 
  ///     to either an ascending or descending sort of provided values
  ///   Inputs:
  ///     complex   : complex associated with filtration
  ///     valuator  : a function which takes a Cell and returns a value
  ///     direction : either "ascending" or "descending" (gives desired ordering of values)
  Filtration ( ImageComplex complex,
               std::function<int64_t(Cell const&)> valuator,
               std::string const& direction ) : complex_(complex) {
    // Initialize the filtration
    typedef std::pair<Cell, int64_t> CellValue;
    filtration_ . reset ( new std::vector<CellValue> );
    for ( auto cell : complex.cells () ) filtration_ -> push_back ( { cell, valuator(cell) } );
    // Prepare the sort
    bool ascending_or_descending;
    if ( direction == "ascending" ) ascending_or_descending = false;
    if ( direction == "descending" ) ascending_or_descending = true;
    if ( direction != "ascending" && direction != "descending" ) {
      throw std::runtime_error("Filtration:" + direction + " is not a valid mode");
    }
    auto comparator = [&](CellValue const& a, CellValue const& b){
        if ( a.second == b.second ) return a.first < b.first;
        return ascending_or_descending != (a.second < b.second);
      };
    std::sort(filtration_->begin(), filtration_->end(), comparator);
    // Index the cells
    cell_indexing_ . reset ( new std::unordered_map<Cell, uint64_t, CellHasher> );
    for ( uint64_t i = 0; i < complex.size(); ++ i ) (*cell_indexing_)[cell(i)] = i; 
  }
  /// complex
  ///   Return the complex the filtration is associated to
  ImageComplex const&
  complex ( void ) const {
    return complex_;
  }
  /// cell
  ///   Return the ith cell in the filtration
  Cell const&
  cell ( uint64_t i ) const {
    return (*filtration_)[i].first;
  }
  /// index
  ///   Given a cell, return its position in the filtration
  uint64_t
  index ( Cell const& cell) const {
    return (*cell_indexing_)[cell];
  } 
  /// value
  ///   Return the value associated with a cell in the filtration
  int64_t 
  value ( uint64_t i ) const {
    return (*filtration_)[i].second;
  }
private:
  ImageComplex complex_;
  std::shared_ptr<std::vector<std::pair<Cell, int64_t>>> filtration_;
  std::shared_ptr<std::unordered_map<Cell, uint64_t, CellHasher>> cell_indexing_;
};

/// SublevelFiltration
///   Overview:
///     Creates a complex and filtration of its cells based on sublevel sets
///   Inputs:
///     image: image object providing width(), height() and data(x,y) methods
///            for 0 <= i < width(), 0 <= j < height() 
///   Outputs:
///     a Filtration object
Filtration
SublevelFiltration ( Image const& image ) {
  // Create N x M full 2D cubical complex
  uint64_t N = image.width();
  uint64_t M = image.height();
  ImageComplex complex(N,M);
  std::unordered_map<Cell, uint64_t, CellHasher> cell_values;
  std::queue<Cell> work_queue;

  // Create 2-cells and assign pixel data to them
  for ( int x = 0; x < N; ++ x ) {
    for ( int y = 0; y < M; ++ y ) {
      Cell cell = Cell(x,y,3); // type 3 = 2-cell
      cell_values[cell] = image.data(x,y);
      work_queue.push(cell);
    }
  }

  // Recursively determine boundary cells and assign values
  //   * Each vertex and edge is assigned the minimum of the values on its coboundary
  //   * The queueing prevents lower dimensional cells from being processed until after
  //     all higher dimensional cells are processed
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

  // Return filtration object
  auto valuation = [&](Cell const& cell){return cell_values[cell];};
  return Filtration ( complex, valuation, "ascending" );
}

/// SuperlevelFiltration
///   Overview:
///     Creates a complex and filtration of its cells based on superlevel sets
///   Inputs:
///     image: image object providing width(), height() and data(x,y) methods
///            for 0 <= i < width(), 0 <= j < height() 
///   Outputs:
///     a Filtration object
Filtration
SuperlevelFiltration ( Image const& image ) {
  // Create a (N-1) x (M-1) cubical complex
  uint64_t N = image.width();
  uint64_t M = image.height();
  ImageComplex complex(N-1,M-1);
  std::unordered_map<Cell, uint64_t, CellHasher> cell_values;
  std::queue<Cell> work_queue;

  // Create vertices and assign pixel data to them
  for ( int x = 0; x <= N; ++ x ) {
    for ( int y = 0; y <= M; ++ y ) {
      Cell cell = Cell(x,y,0); // type 0 = vertex
      cell_values[cell] = image.data(x,y);
      work_queue.push(cell);
    }
  }

  // Recursively determine coboundary cells and assign values
  //   * Each edge and 2-cell is assigned the minimum of the values on its boundary
  //   * The queueing prevents higher dimensional cells from being processed until after
  //     all lower dimensional cells are processed
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
    }
  }

  // Return filtration object
  auto valuation = [&](Cell const& cell){return cell_values[cell];};
  return Filtration ( complex, valuation, "descending" );
}

/// PersistenceViaPHAT
///   Overview:
///     Given a filtration, compute persistent homology using PHAT. 
///   Inputs:
///     filtration : A filtration object
phat::persistence_pairs
PersistenceViaPHAT ( Filtration const& filtration ) {
  // Get reference to complex
  ImageComplex const& complex = filtration.complex();

  // Create boundary matrix object
  phat::boundary_matrix< phat::vector_vector > boundary_matrix;

  // Set the number of columns (i.e. number of cells)
  uint64_t num_cells = complex.size();
  boundary_matrix . set_num_cols(num_cells);

  // Set column for each cell
  for ( phat::index i = 0; i < num_cells; ++ i ) {
    Cell const& cell = filtration.cell(i) ;
    std::vector<phat::index> boundary;
    for ( Cell const& bd_cell : complex.boundary(cell) ) {
      boundary.push_back(filtration.index(bd_cell));
    }
    std::sort(boundary.begin(), boundary.end()); // is this required?
    boundary_matrix . set_dim( i, cell.dimension());    
    boundary_matrix . set_col( i, boundary );
  }
  // Call PHAT
  phat::persistence_pairs pairs;  
  phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );

  return pairs;
}

/// SavePersistenceResults
///   Overview:
///     Save computed persistence results of a filtration to a file
///     This is a specialized function which assumes the associated complex
///     is an "ImageComplex" and saves feature data related to x and y coordinates.
///   Inputs:
///     filtration   : a filtration of a complex
///     pairs        : PHAT output of birth-death cell persistence pairs
///     outfile_name : the name of the file in which to save the results of the calculation
void
SavePersistenceResults ( Filtration const& filtration,
                         phat::persistence_pairs const& pairs, 
                         std::string const& outfile_name ) {
  // Output result
  std::ofstream outfile ( outfile_name );
  outfile << "dim,birth,b_x,b_y,b_z,death,d_x,d_y,d_z\n" ;
  for( int i = 0; i < pairs.get_num_pairs(); ++ i ){
    auto birth_cell_index = pairs.get_pair(i).first;
    auto death_cell_index = pairs.get_pair(i).second;
    auto birth_cell = filtration.cell(birth_cell_index);
    auto birth_value = filtration.value(birth_cell_index);
    auto death_cell = filtration.cell(death_cell_index);
    auto death_value = filtration.value(death_cell_index);
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

/// main
///  Usage: ImagePersistence <image_filename> <output_filename> <mode>
///   <image_filename>  : input image filename
///   <output_filename> : output filename (a CSV file)
///   <mode>            : A mode of execution (either "sub" or "super")
///                       to indicate whether to compute
///                       sublevel or superlevel persistence.
int main(int argc, char *argv[]) {
  // Check command line arguments
  if ( argc != 4 ) {
    std::cerr << help_string << "\n";
    return 1;
  }

  // Read arguments
  std::string infile_name(argv[1]);
  std::string outfile_name(argv[2]);
  std::string mode(argv[3]);
  
  // Check mode argument
  if ( mode != "sub" && mode != "super" ) {
    std::cerr << help_string << "\n";
    return 1;
  }

  // Load image file
  Image image ( infile_name );

  // Construct filtration, an object that stores an ordering of cells in a complex
  Filtration filtration;
  if ( mode == "sub" ) filtration = SublevelFiltration(image);
  if ( mode == "super" ) filtration = SuperlevelFiltration(image);

  // Compute persistence given filtration
  phat::persistence_pairs pairs = PersistenceViaPHAT(filtration);

  // Save results to file
  SavePersistenceResults ( filtration, pairs, outfile_name );

  return 0;
}
