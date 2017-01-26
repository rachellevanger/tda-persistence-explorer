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

#include "common.h"

#include "phat/representations/vector_vector.h"
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/twist_reduction.h>
#include <phat/compute_persistence_pairs.h>

#include "CubicalCell.h"
#include "CubicalComplex.h"
#include "Filtration.h"
#include "Data.h"

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
//   * `CubicalCell`   : class for storing cubical cells
//   * `CubicalComplex`   : class for storing cubical complexes

// Functions
//   * `SublevelFiltration`     : Given an Data object, construct a Filtration object
//   * `SuperlevelFiltration`   : Given an Data object, construct a Filtration object
//   * `PersistenceViaPHAT`     : Given a Filtration object, compute persistence pairs
//   * `SavePersistenceResults` : Given Filtration and PersistencePairs, save result to file                              
//   * `main`                   : Parse command line arguments, create filtration, 
//                                compute persistence, and write output to file

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
  std::cout << "Creating cubical complex.\n";
  CubicalComplex complex(data.resolution());
  std::cout << "Cubical complex created.\n";
  uint64_t D = complex.dimension();

  std::unordered_map<CubicalCell, uint64_t> cell_values;
  std::queue<CubicalCell> work_queue;

  // Create top-cells and assign pixel data to them
  for ( auto cell : complex ) { // A more efficient iteration pattern would avoid needing "if dim == D"
    std::cout << "Inspect cell " << cell << "\n";
    if ( cell.dimension() == D ) { 
      cell_values[cell] = data(cell.coordinates());
      work_queue.push(cell);
    }
  }

  std::cout << "Pixel data assigned to top cells.\n";

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

  std::cout << "Pixel data assigned to all cells.\n";
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
  std::cout << "Data loaded.\n";
  std::cout << data << "\n";

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
