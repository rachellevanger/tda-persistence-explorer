// DataPersistence.cpp
//  Compute persistent homology of multidimensional array data using PHAT.
//  MIT LICENSE 2016 Jacek Cyranka, Shaun Harker
// 
// Revision History:
//   2016-11-21 Shaun Harker 
//      * refactored
//      * fixed non-square image bug 
//      * improved make system
//   2016-12-15 Shaun Harker
//      * generalized to higher dimensions

// USE FLOATING POINT DATA
// x y z value \n FORMAT
// DON'T ASSUME SORTED PROPERLY.

#include <memory>
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

#include <stdint.h>

using namespace cimg_library;

std::string help_string = 
" Usage: ImagePersistence <input_filename> <output_filename> <mode> \n"
"   <input_filename>  : input filename \n"
"   <output_filename> : output filename (a CSV file)\n"
"   <mode>            : A mode of execution (either \"sub\" or \"super\") \n"
"                       to indicate whether to compute\n"
"                       sublevel or superlevel persistence.\n";


// Overview:
//  1. The input filename is loaded into a "Data" object
//  2. Based on "mode", either a sublevel or superlevel filtration is constructed.
//  3. Persistence is calculated for the filtration
//  4. Results are saved to the requested output file.

// Implementation overview:
// Classes
//   * `Data`        : class for loading and accessing multidimensional data
//   * `Filtration`   : class for storing filtrations of complexes
// Functions
//   * `SublevelFiltration`     : Given an Data object, construct a Filtration object
//   * `SuperlevelFiltration`   : Given an Data object, construct a Filtration object
//   * `PersistenceViaPHAT`     : Given a Filtration object, compute persistence pairs
//   * `SavePersistenceResults` : Given Filtration and PersistencePairs, save result to file                              
//   * `main`                   : Parse command line arguments, create filtration, 
//                                compute persistence, and write output to file

/// Data
///   This class is used to load and access data
class Data {
public:
  /// Data
  ///   Default constructor
  Data ( void ) {}

  /// Data
  ///   Construct Data object from file.
  ///   Notes: If the filename extension looks like an image, it uses CImg to load.
  ///          Otherwise, it uses the `loadData` method
  Data ( std::string const& input_filename ) {
    const auto list_of_image_extensions = { "png", "bmp", "jpeg", "jpg", "gif" };
    // Get filename extension (part after last ".")
    auto ext = input_filename.substr(fn.find_last_of(".") + 1);
    // We check if `ext` is in `list_of_image_extensions`
    // C++ has atrocious semantics for this:
    auto begin = list_of_image_extensions.begin();
    auto end = list_of_image_extensions.end();
    bool input_file_is_an_image = std::find(begin, end, ext) != end;
    // Has an image extension
    if ( input_file_is_an_image ) {
      loadImage(input_filename);
    } else {
      loadData(input_filename);
    }
  }

  /// loadImage
  ///   Load image given image_filename
  void
  loadImage ( std::string const& image_filename ) { 
    auto image = std::make_shared<CImg<unsigned char>> (image_filename.c_str());
    resolution_ . resize ( 2 );
    resolution_[0] = image -> width();
    resolution_[1] = image -> height();   
    data_ = [=](std::vector<uint64_t> const& coordinates) {return (*image)(coordinates[0],coordinates[1],0,1);} 
  }

  /// loadData
  ///   Load data given data_filename
  void
  loadData ( std::string const& data_filename ) { 
    uint64_t dimension = 0;
    // TODO: get input dimension
    resolution_ . resize ( dimension );
    // TODO: set resolution_[d] for d = 0,1,...,D-1
    // TODO: load data from file and put it in a structure wrapped in a 
    //       shared pointer (use std::make_shared) and capture it in the
    //       following lambda
    data_ = [=](std::vector<uint64_t> const& coordinates) {return 0;}
  }
  /// resolution
  ///   Give extent of data in each dimension
  ///   i.e.  0 <= x_i < resolution()[i] for i = 0 ... D-1
  std::vector<uint64_t>
  resolution ( void ) const {
    return resolution_;
  }
  /// data
  ///   Given 0 <= x < width() and 0 <= y < height(), 
  ///   return pixel data at position (x,y)
  uint64_t
  data ( std::vector<uint64_t> const& coordinates ) const {
    return data_(coordinates);
  }
private:
  std::function<uint64_t(std::vector<uint64_t>const&)> data_;
  std::vector<uint64_t> resolution_;
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
  Filtration ( CubicalComplex complex,
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
  CubicalComplex const&
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
  CubicalComplex complex_;
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
SublevelFiltration ( Data const& data ) {
  // Create cubical complex
  CubicalComplex complex(data.resolution());
  uint64_t D = complex.dimension();

  std::unordered_map<Cell, uint64_t, CellHasher> cell_values;
  std::queue<Cell> work_queue;

  // Create top-cells and assign pixel data to them
  for ( cell : complex.cells(D) ) {
    cell_values[cell] = image.data(x,y);
    work_queue.push(cell);
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
  // Create cubical complex
  CubicalComplex complex(data.sizes());
  uint64_t D = complex.dimension();
  std::unordered_map<Cell, uint64_t, CellHasher> cell_values;
  std::queue<Cell> work_queue;

  // Create vertices and assign pixel data to them
  for ( cell : complex.cells(0) ) {
    cell_values[cell] = image.data(x,y);
    work_queue.push(cell);
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
///  TODO: make non-specific to 2D and 3D
///   Overview:
///     Save computed persistence results of a filtration to a file
///     This is a specialized function which assumes the associated complex
///     is an "ImageComplex" and saves feature data related to x, y, and z coordinates (as available)
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
    auto birth_coordinates = [&](uint64_t d) { return birth_cell.coordinates().size() <= d ? 0 : birth_cell.coordinates()[i]; }
    auto death_coordinates = [&](uint64_t d) { return death_cell.coordinates().size() <= d ? 0 : death_cell.coordinates()[i]; }

    if ( birth_value == death_value ) continue;
    outfile << birth_cell.dimension() << ", "
            << birth_value << ", "
            << birth_coordinates[0] << ", "
            << birth_coordinates[1] << ", "
            << birth_coordinates[2] << ", "
            << death_value << ", "
            << death_coordinates[0] << ", "
            << death_coordinates[1] << ", "
            << death_coordinates[2] << "\n";
  }
}

/// main
///   Entry point of program
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
  Data data ( infile_name );

  // Construct filtration, an object that stores an ordering of cells in a complex
  Filtration filtration;
  if ( mode == "sub" ) filtration = SublevelFiltration(data);
  if ( mode == "super" ) filtration = SuperlevelFiltration(data);

  // Compute persistence given filtration
  phat::persistence_pairs pairs = PersistenceViaPHAT(filtration);

  // Save results to file
  SavePersistenceResults ( filtration, pairs, outfile_name );

  return 0;
}
