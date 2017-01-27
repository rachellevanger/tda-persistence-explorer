/// Filtration.h
/// Shaun Harker 2017-01-25-2318
/// MIT LICENSE

#ifndef DATAPERSISTENCE_FILTRATION_H
#define DATAPERSISTENCE_FILTRATION_H

#include "common.h"

#include "CubicalComplex.h"

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

/// Filtration_
///   Underlying data for Filtration class
struct Filtration_ {
  CubicalComplex complex_;
  std::vector<uint64_t> original_index_from_filtration_index_;
  std::vector<uint64_t> filtration_index_from_original_index_;
  std::vector<double> values_; // indexed by filtration ordering.
};

/// Filtration
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
  ///     direction : either "ascending" or "descending" (gives desired ordering of values)
  Filtration ( CubicalComplex complex,
               std::vector<double> original_values,
               std::string const& direction ) {
    assign(complex, original_values, direction);
  }

  // Filtration
  //   constructor
  void
  assign ( CubicalComplex complex_in,
           std::vector<double> original_values,
           std::string const& direction ) {
    data_ = std::make_shared<Filtration_>();
    complex_() = complex_in;
    // Prepare the sort
    bool ascending_or_descending;
    if ( direction == "ascending" ) ascending_or_descending = false;
    if ( direction == "descending" ) ascending_or_descending = true;
    if ( direction != "ascending" && direction != "descending" ) {
      throw std::runtime_error("Filtration:" + direction + " is not a valid mode");
    }
    auto comparator = [&](uint64_t a, uint64_t b) {
      auto value_a = original_values[a];
      auto value_b = original_values[b];
      if ( value_a == value_b ) return a < b;
      return ascending_or_descending != (value_a < value_b);
    };
    // Setup "original_index_from_filtration_index_"
    auto & X = original_index_from_filtration_index_();
    X.resize(complex_().size());     // X == [0,0,0,...,0]
    std::iota(X.begin(), X.end(), 0 ); // X == [0,1,2,...,complex.size()-1]
    std::sort(X.begin(), X.end(), comparator);
    // Setup "filtration_index_from_original_index_"
    auto & Y = filtration_index_from_original_index_();
    for ( uint64_t x = 0; x < X.size(); ++ x) Y[X[x]] = x; 
    // Setup "values_"
    auto & V = values_();
    V.resize(complex_().size());
    for ( uint64_t i = 0; i < V.size(); ++ i) V[i] = original_values[X[i]]; 
  }

  /// complex
  ///   Return the complex the filtration is associated to
  CubicalComplex const&
  complex ( void ) const {
    return data_ -> complex_;
  }

  /// cell
  ///   Return the ith cell in the filtration
  uint64_t
  original ( uint64_t filtered ) const {
    return original_index_from_filtration_index_()[filtered];
  }

  /// index
  ///   Given a cell, return its position in the filtration
  uint64_t
  filtered ( uint64_t original) const {
    return filtration_index_from_original_index_()[original];
  } 

  /// value
  ///   Return the value associated with an index in the filtration
  double 
  value ( uint64_t filtered ) const {
    return values_(filtered);
  }

private:

  std::shared_ptr<Filtration_> data_;

  /// complex_
  ///   accessor for data_ -> complex_
  CubicalComplex &
  complex_ ( void ) {
    return data_ -> complex_;
  }

  /// original_order_from_filtration_order_
  ///   accessor for data_ -> original_order_from_filtration_order_
  std::vector<uint64_t> &
  original_index_from_filtration_index_ ( void ) const {
    return data_ -> original_order_from_filtration_order_;
  }

  /// filtration_order_from_original_order_
  ///   accessor for data_ -> filtration_order_from_original_order_
  std::vector<uint64_t> &
  filtration_index_from_original_index_ ( void ) const {
    return data_ -> filtration_order_from_original_order_;
  }

  /// values_
  ///   accessor for data_ -> values_
  double &
  values_ ( void ) const {
    return data_ -> values_;
  }

};

#endif
