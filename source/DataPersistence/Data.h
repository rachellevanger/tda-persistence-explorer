/// Data.h
/// Shaun Harker 2017-01-25-2321
/// MIT LICENSE

#ifndef DATAPERSISTENCE_DATA_H
#define DATAPERSISTENCE_DATA_H

#include "common.h"

#include "CImg.h"
using namespace cimg_library;

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
  double
  operator () ( std::vector<uint64_t> const& coordinates ) const {
    return data_(coordinates);
  }

  /// operator <<
  ///   print summary
  friend std::ostream & operator << ( std::ostream & stream, Data const& stream_me ) {
    stream << "Data([";
    for ( auto x : stream_me.resolution() ) stream << x << ",";
    stream << "])";
    return stream;
  }

private:
  std::function<double(std::vector<uint64_t>const&)> data_;
  std::vector<uint64_t> resolution_;
};

#endif
