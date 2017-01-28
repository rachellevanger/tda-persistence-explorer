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
//      * generalized to higher dimensions (WIP)
//   2017-01-27 Shaun Harker
//      * optimized using arrays

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

  // Algorithm Description. 
  //
  //   CubicalComplex
  //     Our CubicalComplex structure is constructed given a list of sizes, indicating
  //     "boxes across", for each dimension. There is a "gotcha" here:
  //     the CubicalComplex does not encode the cells on the right-boundary of
  //     the right-most boxes. The reason for this is because if we did, then
  //     cells of different "shapes" (i.e. signature of which dimensions have extent)
  //     would come in different numbers. We can see this in the 1D example.
  //     Rather than encode   * ----- * ----- *,  which has 3 0-cells and 2 1-cells,
  //     we instead have only * ----- * ------ , and there is a 1-to-1 correspondence
  //     between 0-cells and 1-cells. In higher dimensions, we no longer have
  //     a one-to-one correspondence between cells of different dimensions, but we
  //     do keep a one-to-one correspondence between cells of different types. It is
  //     just that some dimensions have more than one corresponding type. For instance
  //     in a 2D complex we have both horizontal edges and vertical edges. So if we make
  //     a CubicalComplex with sizes (2,2) we end up with the following:
  //
  //     | ### | ###
  //     | ### | ###
  //     * --- * ---
  //     | ### | ### 
  //     | ### | ###
  //     * --- * ---
  //
  //    Note the homology of such complexes is trivial. 
  //    In order to use these complexes to handle the case where we actually do want the
  //    right edges present, we use what we call the "wrapping trick": we make the
  //    CubicalComplex one box wider in each dimension. 
  //    For the current application we want to find filtrations; we will choose the 
  //    valuation of the extra cells on the wrapper edge that we don't want to be infinity,
  //    and then explicitly ignore such cells later on.
  //
  //   Cells
  //     A cell is described by is coordinates, a D-tuple, and its type.
  //     The "type" is a numeric code which indicates which dimensions it has extent in.
  //     However, is does not do this in the most convenient way possible. 
  //     An alternative notation we call "shape" is more convenient. 
  //     The "shape" of a cell is an integer in [0, 1 << D) with the ith bit on 
  //       iff there is extent in the ith dimension.
  //     Since the numeric ordering of shapes does not respect the ordering of
  //     dimensions we do not use "shape" directly for indexing. Instead we
  //     use "type", which is the reindexing of shape by performing a stable sort
  //     of shapes when we compare by dimension.
  //
  //    Cell Index
  //      For speed we refer to cells via an integer. This integer is the rank
  //      of its tuple (coord[0], coord[1], ..., coord[D}, type) in the lexicographically
  //      sorted list of all such tuples (where the latter coordinates in the tuple are
  //      sorted first, i.e. are more significant)
  //
  //    Cell Values
  //      To associate a value to each cell, as we need to do for the purposes of
  //      computing a sub- or super- level filtration, we need to create a map
  //      from cells to doubles. The fastest thing to do is create an array of doubles
  //      which is indexed in precisely the same manner that cells are indexed.
  //      Accordingly, may regard the values V as multi-dimensional array
  //      of size [S[0], S[1], ..., S_D] :=  
  //        [complex.sizes()[0], ..., complex.sizes()[D-1], 1 << D].
  //
  //    Slices
  //      In order to perform operations efficiently we exploit "slice operations"
  //      on multidimensional arrays.  A "slice" of a multidimensional array is 
  //      the set of array indices when we restrict each component x[i] to be 
  //      between a[i] and b[i], i.e 0 <= a[i] <= x[i] < b[i] <= V[i].
  //      This is often written as V[a[0]:b[0], a[1]:b[1], ..., a[D]:b[D]]
  //      A neat notation trick is to write V[:,:,...,:] for the "total" slice,
  //      i.e. the whole array, and use negative numbers on the right bound to
  //      count backwards, e.g. for a vector V[:-1] is all but the last entry.

  // We propagate the value information from "data" to
  // "V" in the following manner.

  // 0. We initialize the entire array with +inf.
  // 1. First, we copy the data directly onto the D-dimensional cells
  //    in the slice V[:-1,:-1,...,:-1]  (using python-esque slice notation)
  //    The reason for the -1's is because of the wrap; we do not have data to 
  //    put into the D-dimensional extra cells in our wrapping.
  // 2. Second, we do a sequence of operations which are all of the form
  //     V[Slice1] = std::min ( V[Slice1], V[Slice2] )     (*)
  //    For this to make sense the slices Slice1 and Slice2 must be commensurate
  //    in the obvious way.
  // 
  //    By a judicious pattern of choosing such pairs of slices we can ultimately
  //    propagate the value data on the D-cells to all the lower dimensional cells
  //    with the end goal of achieving that
  //     V[cell] = {  data[cell]  if x is a D-cell not in the wrapped border
  //               {  min { V[y] : x is a face of y } otherwise
  //
  //    It is hopefully straightforward to see that one way to accomplish this
  //    would be to perform, for each type of cell, a (*) operation correspond to
  //    each type of collapse. For example a 3-cell can be collapsed onto 6 faces.
  //    We could then do 6 operations (*) in order to ensure each of those boundaries
  //    underwent a boundary_cell_value <- min(boundary_cell_value, coboundary_value).
  //    operation.
  //
  //    In fact, this is overkill. We leave it as an exercise, but it can be shown that
  //    it suffices that for each class of cells of a given shape, there exists a 
  //    "parent" shape which collapses onto the "child" shape in dimension d, and we
  //    have performed two (*) operations corresponding to propagating value data
  //    from the parent shape cells to the child shape cells in both the left and right
  //    collapses of dimension d. This ends up giving us
  //
  //     (LEFT-PROPAGATE)
  //     V[child-slice] = std::min(V[child-slice], V[parent-slice])  
  //
  //     (RIGHT-PROPAGATE)
  //     V[child-slice \ left-most-cells-in-dim-d] = 
  //               std::min(V[child-slice \ left-most-cells-in-dim-d], 
  //                        V[parent-slice \ right-most-cells-in-dim-d])
  //
  //    It isn't hard to arrange such a sequence of operations; we can do any kind of
  //    spanning tree algorithm on the hypercube of shapes (each vertex corresponds to 
  //    the binary bit code which encodes shape as described above) starting from the 
  //    top dimensional shape (will all-1's binary code). 
  // 
  //    The following algorithm implements these ideas.
  //    It could stand to be written more cleanly. 
  // Create cubical complex using "wrapping trick"
  auto incremented_resolution = data.resolution();
  for ( auto & size : incremented_resolution ) ++ size; // wrapping trick
  CubicalComplex complex(incremented_resolution);
  uint64_t D = complex.dimension();

  // Initialize array to hold cell values for 
  // sublevel filtration with "+inf" values:
  std::vector<double> V (complex.size(), std::numeric_limits<double>::infinity());

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
        V[shape_begin + offset] = data[i];
        ++ i;
      }
    } else {
      // Propagate step
      uint64_t parent = parent_shape[shape];
      uint64_t bd_shape_begin = complex.shape_begin(shape);
      uint64_t cbd_shape_begin = complex.shape_begin(parent);
      std::vector<uint64_t> bottom (D);
      std::vector<uint64_t> bd_L(D);
      std::vector<uint64_t> cbd_U = complex.sizes();
      std::vector<uint64_t> top = complex.sizes();
      // LEFT-PROPAGATE
      for ( auto offset : Slice(bottom, bottom, top, top)) {
        V[bd_shape_begin + offset] = std::min ( V[ bd_shape_begin + offset ], V[ cbd_shape_begin + offset ]);
      }
      // RIGHT-PROPAGATE
      uint64_t d = 0; { uint64_t t = parent ^ shape; while ( t >>= 1 ) ++d; } // collapse dimension
      ++ bd_L[d]; -- cbd_U[d];
      auto bd_slice = Slice(bottom, bd_L, top, top);
      auto cbd_slice = Slice(bottom, bottom, cbd_U, top);
      for ( auto it1 = bd_slice.begin(), it2 = cbd_slice.begin(); 
            it1 != bd_slice.end() && it2 != cbd_slice.end(); ++it1, ++it2) {

        V[ bd_shape_begin + *it1 ] = std::min ( V[ bd_shape_begin + *it1 ], V[ cbd_shape_begin + *it2 ]);
      }
    }
    // Determine other propagation steps (Hypercube Traversal, breadth-first)
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
  return Filtration ( complex, V, "ascending" );
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
  std::vector<double> V(complex.size());

  // // Create vertices and assign pixel data to them
  // for ( auto cell : complex(0) ) {
  //   V[cell] = data(cell.coordinates()); // don't want to pass through coordinates here
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
  //     if ( V.count(cbd_cell) == 0 ) { 
  //       work_queue.push(cbd_cell);
  //       V[cbd_cell] = V[work_cell];
  //     } else {
  //       V[cbd_cell] = std::min(V[cbd_cell], V[work_cell]);
  //     }
  //   }
  // }

  // Return filtration object
  return Filtration ( complex, V, "descending" );
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
  uint64_t num_cells = filtration.finite_size();
  boundary_matrix . set_num_cols(num_cells);

  // Set column for each cell
  for ( phat::index i = 0; i < num_cells; ++ i ) {
    uint64_t original_index = filtration.original(i) ;
    std::vector<phat::index> boundary;
    for ( uint64_t bd_cell : complex.boundary(original_index) ) {
      boundary.push_back(filtration.filtered(bd_cell));
    }
    std::sort(boundary.begin(), boundary.end()); // is this required?
    boundary_matrix . set_dim (i, complex.dimension(original_index));    
    boundary_matrix . set_col (i, boundary );
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
