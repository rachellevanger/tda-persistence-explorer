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

#include "CubicalComplex.h"

#include <stdint.h>

using namespace cimg_library;

std::string help_string = 
" Usage: DataPersistence <input_filename> <output_filename> <mode> \n"
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
    auto ext = input_filename.substr(input_filename.find_last_of(".") + 1);
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
    data_ = [=](std::vector<uint64_t> const& coordinates) {return (*image)(coordinates[0],coordinates[1],0,1);}; 
  }

  /// loadData
  ///   Load data given data_filename
  void
  loadData ( std::string const& data_filename ) { 
    // Open the file for reading
    std::ifstream infile (data_filename);

    // Check to see if file was loaded:
    if ( ! infile.good() ) throw std::runtime_error("Error parsing " + data_filename + ": could not open file");
    typedef std::vector<double> Point;
    std::vector<Point> data;

    // Read line-by-line
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        // Convert line to a point
        Point p;
        double x; while (iss >> x) p.push_back(x);
        // Add point to data
        data.push_back(p);
    }

    // Check to see that data has been parsed
    if ( data.size() == 0 ) throw std::runtime_error("Error parsing " + data_filename + ": empty or incorrect format");

    // Extract the dimension of the data. The last coordinate is value, not space,
    // so is not counted:
    uint64_t dimension = data[0].size() - 1;

    // Construct the totally ordered sets X_0, ..., X_{dim-1}
    // where X_i is the set of values occurring in the ith coordinate
    std::vector<std::set<double>> X(dimension);
    for ( auto point : data ) {
      for ( uint64_t i = 0; i < dimension; ++ i ) {
        X[i].insert(point[i]);
      }
    }

    // Create indexing functions I_i : X_i -> {0,1,2,..., card X_i - 1} 
    std::vector<std::unordered_map<double, uint64_t>> I(dimension);
    for ( uint64_t i = 0; i < dimension; ++ i ) {
      uint64_t j = 0;
      // Iterate through X_i from least to greatest.
      for ( auto value : X[i] ) {
        I[i][value] = j;
        ++ j;
      }
    }

    // Set resolution sizes "resolution_"
    // and set the place_values
    // and also compute N := card \Prod_{i=0}^{dim-1} X_i
    resolution_ . resize ( dimension );
    std::vector<uint64_t> place_values ( dimension );
    uint64_t N = 1;
    for ( uint64_t i = 0; i < dimension; ++ i ) {
      resolution_[i] = X[i].size();
      place_values[i] = (i == 0) ? 1 : resolution_[i-1] * place_values[i-1];
      N *= resolution_[i];
    }

    // Check if there actually are N points of data in file
    if ( data.size() < N ) {
      throw std::runtime_error("Error parsing " + data_filename + ": not every grid point is given a value.");
    }
    if ( data.size() > N ) {
      throw std::runtime_error("Error parsing " + data_filename + ": redundant data points.");
    }

    // Sort the data
    // A point is (x_0, x_1, ..., x_{dim-1}, value)
    // and we want to sort lexicographically with "value" ignored
    // and  x_{dim-1} being most significant.
    //  (Implementation is to use reverse iterators and increment the 'rbegin' iterators)
    auto point_compare = [](Point const& lhs, Point const& rhs) {
      return std::lexicographical_compare( ++lhs.rbegin(), lhs.rend(), ++ rhs.rbegin(), rhs.rend());
    };
    std::sort ( data.begin(), data.end(), point_compare );

    // Check for redundant data points
    for ( uint64_t i = 1; i < N; ++ i ) {
      Point p = data[i-1];
      Point q = data[i];
      // Compare p[0:dim] to q[0:dim] and throw if equal
      p.pop_back();
      q.pop_back();
      if ( p == q ) {
        throw std::runtime_error("Error parsing " + data_filename + ": redundant data points.");
      }
    }

    // Store the values in an array.
    // Since the data is sorted, the indexing of this array is meaningful.
    auto values = std::make_shared<std::vector<double>>();
    for ( auto point : data ) {
      values -> push_back ( point[dimension] );
    }

    // Store the data
    data_ = [=](std::vector<uint64_t> const& coordinates) {
      // Convert coordinates to index
      uint64_t idx = 0;
      for ( int i = 0; i < dimension; ++ i ) {
        idx += coordinates[i] * place_values[i];
      }
      return (*values)[idx];
    };
  }

  /// resolution
  ///   Give extent of data in each dimension
  ///   i.e.  0 <= x_i < resolution()[i] for i = 0 ... D-1
  std::vector<uint64_t>
  resolution ( void ) const {
    return resolution_;
  }

  /// operator ()
  ///   Given 0 <= x < width() and 0 <= y < height(), 
  ///   return pixel data at position (x,y)
  uint64_t
  operator () ( std::vector<uint64_t> const& coordinates ) const {
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
  ///     valuator  : a function which takes a CubicalCell and returns a value
  ///     direction : either "ascending" or "descending" (gives desired ordering of values)
  Filtration ( CubicalComplex complex,
               std::function<int64_t(CubicalCell const&)> valuator,
               std::string const& direction ) : complex_(complex) {
    // Initialize the filtration
    typedef std::pair<CubicalCell, int64_t> CubicalCellValue;
    filtration_ . reset ( new std::vector<CubicalCellValue> );
    for ( auto cell : complex ) filtration_ -> push_back ( { cell, valuator(cell) } );
    // Prepare the sort
    bool ascending_or_descending;
    if ( direction == "ascending" ) ascending_or_descending = false;
    if ( direction == "descending" ) ascending_or_descending = true;
    if ( direction != "ascending" && direction != "descending" ) {
      throw std::runtime_error("Filtration:" + direction + " is not a valid mode");
    }
    auto comparator = [&](CubicalCellValue const& a, CubicalCellValue const& b){
        if ( a.second == b.second ) return a.first < b.first;
        return ascending_or_descending != (a.second < b.second);
      };
    std::sort(filtration_->begin(), filtration_->end(), comparator);
    // Index the cells
    cell_indexing_ . reset ( new std::unordered_map<CubicalCell, uint64_t> );
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
  CubicalCell const&
  cell ( uint64_t i ) const {
    return (*filtration_)[i].first;
  }

  /// index
  ///   Given a cell, return its position in the filtration
  uint64_t
  index ( CubicalCell const& cell) const {
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
  std::shared_ptr<std::vector<std::pair<CubicalCell, int64_t>>> filtration_;
  std::shared_ptr<std::unordered_map<CubicalCell, uint64_t>> cell_indexing_;
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

  std::unordered_map<CubicalCell, uint64_t> cell_values;
  std::queue<CubicalCell> work_queue;

  // Create top-cells and assign pixel data to them
  for ( auto cell : complex ) { // A more efficient iteration pattern would avoid needing "if dim == D"
    if ( cell.dimension() == D ) { 
      cell_values[cell] = data(cell.coordinates());
      work_queue.push(cell);
    }
  }

  // Recursively determine boundary cells and assign values
  //   * Each vertex and edge is assigned the minimum of the values on its coboundary
  //   * The queueing prevents lower dimensional cells from being processed until after
  //     all higher dimensional cells are processed
  while ( not work_queue . empty () ) {
    CubicalCell work_cell = work_queue . front (); 
    work_queue . pop ();
    for ( CubicalCell const& bd_cell : complex.boundary(work_cell) ) {
      if ( cell_values.count(bd_cell) == 0 ) { 
        work_queue.push(bd_cell);
        cell_values[bd_cell] = cell_values[work_cell];
      } else {
        cell_values[bd_cell] = std::min(cell_values[bd_cell], cell_values[work_cell]);
      }
    }
  }

  // Return filtration object
  auto valuation = [&](CubicalCell const& cell){return cell_values[cell];};
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
SuperlevelFiltration ( Data const& data ) {
  // Create cubical complex
  auto decremented_resolution = data.resolution();
  for ( auto & size : decremented_resolution ) --size;
  CubicalComplex complex(decremented_resolution);
  uint64_t D = complex.dimension();
  std::unordered_map<CubicalCell, uint64_t> cell_values;
  std::queue<CubicalCell> work_queue;

  // Create vertices and assign pixel data to them
  for ( auto cell : complex ) { // A more efficient iteration pattern would avoid needing "if dim == 0"
    if ( cell.dimension() == 0 ) {
      cell_values[cell] = data(cell.coordinates());
      work_queue.push(cell);
    }
  }

  // Recursively determine coboundary cells and assign values
  //   * Each edge and 2-cell is assigned the minimum of the values on its boundary
  //   * The queueing prevents higher dimensional cells from being processed until after
  //     all lower dimensional cells are processed
  while ( not work_queue . empty () ) {
    CubicalCell work_cell = work_queue . front (); 
    work_queue . pop ();
    for ( CubicalCell const& cbd_cell : complex.coboundary(work_cell) ) {
      if ( cell_values.count(cbd_cell) == 0 ) { 
        work_queue.push(cbd_cell);
        cell_values[cbd_cell] = cell_values[work_cell];
      } else {
        cell_values[cbd_cell] = std::min(cell_values[cbd_cell], cell_values[work_cell]);
      }
    }
  }

  // Return filtration object
  auto valuation = [&](CubicalCell const& cell){return cell_values[cell];};
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
  CubicalComplex const& complex = filtration.complex();

  // Create boundary matrix object
  phat::boundary_matrix< phat::vector_vector > boundary_matrix;

  // Set the number of columns (i.e. number of cells)
  uint64_t num_cells = complex.size();
  boundary_matrix . set_num_cols(num_cells);

  // Set column for each cell
  for ( phat::index i = 0; i < num_cells; ++ i ) {
    CubicalCell const& cell = filtration.cell(i) ;
    std::vector<phat::index> boundary;
    for ( CubicalCell const& bd_cell : complex.boundary(cell) ) {
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
///     is an "CubicalComplex" and saves feature data related to x, y, and z coordinates (as available)
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
    auto birth_coordinates = [&](uint64_t d) { return birth_cell.coordinates().size() <= d ? 0 : birth_cell.coordinates()[i]; };
    auto death_coordinates = [&](uint64_t d) { return death_cell.coordinates().size() <= d ? 0 : death_cell.coordinates()[i]; };

    if ( birth_value == death_value ) continue;
    outfile << birth_cell.dimension() << ", "
            << birth_value << ", "
            << birth_coordinates(0) << ", "
            << birth_coordinates(1) << ", "
            << birth_coordinates(2) << ", "
            << death_value << ", "
            << death_coordinates(0) << ", "
            << death_coordinates(1) << ", "
            << death_coordinates(2) << "\n";
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
