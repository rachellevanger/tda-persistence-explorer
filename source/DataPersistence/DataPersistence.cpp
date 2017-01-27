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

#include "CubicalComplex.h"
#include "Filtration.h"
#include "Data.h"
#include "Slice.h"

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
  auto incremented_resolution = data.resolution();
  for ( auto & size : incremented_resolution ) ++ size;
  CubicalComplex complex(incremented_resolution);
  //CubicalComplex complex(data.resolution());
  uint64_t D = complex.dimension();
  std::vector<double> cell_values (complex.size(), std::numeric_limits<double>::infinity());
  std::queue<uint64_t> shape_queue;
  std::vector<uint64_t> parent_shape ( 1L << D );
  shape_queue.push( (1L << D) - 1);
  while ( not shape_queue . empty () ) {
    uint64_t shape = shape_queue . front ();
    shape_queue . pop ();
    if ( parent_shape[shape] == 0 ) { 
      // Initialize data.
      uint64_t shape_begin = complex.shape_begin(shape);
      std::vector<uint64_t> L(D);
      std::vector<uint64_t> U = complex.sizes();
      for ( auto & x : U ) --x;
      uint64_t i = 0;
      for ( auto offset : Slice(L, L, U, complex.sizes())) {
        cell_values[shape_begin + offset] = data[i];
        ++ i;
      }
    } else {
      uint64_t parent = parent_shape[shape];
      uint64_t bd_shape_begin = complex.shape_begin(shape);
      uint64_t cbd_shape_begin = complex.shape_begin(parent);
      std::vector<uint64_t> bottom (D);
      std::vector<uint64_t> bd_L(D);
      std::vector<uint64_t> cbd_U = complex.sizes();
      std::vector<uint64_t> top = complex.sizes();
      // Step 1. Vanilla step
      for ( auto offset : Slice(bottom, bottom, top, top)) {
        cell_values[bd_shape_begin + offset] = std::min ( cell_values[ bd_shape_begin + offset ], cell_values[ cbd_shape_begin + offset ]);
      }
      // Step 2. Offset step
      // So the bd gets a 1-value on its L, and the cbd gets a -1 on its U.
      uint64_t d = 0; { uint64_t t = parent ^ shape; while ( t ) {t >>= 1; ++d;} } // collapse dimension
      ++ bd_L[d]; -- cbd_U[d];
      auto bd_slice = Slice(bottom, bd_L, top, top);
      auto cbd_slice = Slice(bottom, bottom, cbd_U, top);
      for ( auto it1 = bd_slice.begin(), it2 = cbd_slice.begin(); 
            it1 != bd_slice.end() && it2 != cbd_slice.end(); ++it1, ++it2) {
        cell_values[ bd_shape_begin + *it1 ] = std::min ( cell_values[ bd_shape_begin + *it1 ], cell_values[ cbd_shape_begin + *it2 ]);
      }
    }
    for ( uint64_t d = 0, bit = 1 ; d < D; ++ d, bit <<= 1 ) {
      if ( not ( shape & bit ) ) continue;
      uint64_t child_shape = shape ^ bit;
      if ( parent_shape[child_shape] == 0 ) {
        parent_shape[child_shape] = shape;
        shape_queue.push(child_shape);
      }
    }
  }

  // Return filtration object
  return Filtration ( complex, cell_values, "ascending" );
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
  // auto decremented_resolution = data.resolution();
  // for ( auto & size : decremented_resolution ) --size;
  // CubicalComplex complex(decremented_resolution);
  CubicalComplex complex(data.resolution());
  uint64_t D = complex.dimension();
  std::vector<double> cell_values(complex.size());

  // // Create vertices and assign pixel data to them
  // for ( auto cell : complex(0) ) {
  //   cell_values[cell] = data(cell.coordinates()); // don't want to pass through coordinates here
  //   work_queue.push(cell);
  // }

  // // Recursively determine coboundary cells and assign values
  // //   * Each edge and 2-cell is assigned the minimum of the values on its boundary
  // //   * The queueing prevents higher dimensional cells from being processed until after
  // //     all lower dimensional cells are processed
  // while ( not work_queue . empty () ) {
  //   uint64_t work_cell = work_queue . front (); 
  //   work_queue . pop ();
  //   for ( uint64_t cbd_cell : complex.coboundary(work_cell) ) {
  //     if ( cell_values.count(cbd_cell) == 0 ) { 
  //       work_queue.push(cbd_cell);
  //       cell_values[cbd_cell] = cell_values[work_cell];
  //     } else {
  //       cell_values[cbd_cell] = std::min(cell_values[cbd_cell], cell_values[work_cell]);
  //     }
  //   }
  // }

  // Return filtration object
  return Filtration ( complex, cell_values, "descending" );
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
    uint64_t original_index = filtration.original(i) ;
    std::vector<phat::index> boundary;
    for ( uint64_t bd_cell : complex.boundary(original_index) ) {
      boundary.push_back(filtration.filtered(bd_cell));
    }
    std::sort(boundary.begin(), boundary.end()); // is this required?
    boundary_matrix . set_dim( i, complex.dimension(original_index));    
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
    auto birth_cell = filtration.original(birth_cell_index);
    auto birth_value = filtration.value(birth_cell_index);
    auto death_cell = filtration.original(death_cell_index);
    auto death_value = filtration.value(death_cell_index);
    auto coordinates = [&](uint64_t cell) { return filtration.complex().coordinates(cell); };
    auto birth_coordinates = [&](uint64_t d) { auto c = coordinates(birth_cell); return c.size() <= d ? 0 : c[d]; };
    auto death_coordinates = [&](uint64_t d) { auto c = coordinates(death_cell); return c.size() <= d ? 0 : c[d]; };
    if ( birth_value == death_value ) continue;
    outfile << filtration.complex().dimension(birth_cell) << ", ";
    outfile << birth_value << ", ";
    outfile << birth_coordinates(0) << ", ";
    outfile << birth_coordinates(1) << ", ";
    outfile << birth_coordinates(2) << ", ";
    outfile << death_value << ", ";
    outfile << death_coordinates(0) << ", ";
    outfile << death_coordinates(1) << ", ";
    outfile << death_coordinates(2) << "\n";
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
