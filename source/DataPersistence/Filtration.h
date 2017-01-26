/// Filtration.h
/// Shaun Harker 2017-01-25-2318
/// MIT LICENSE

#ifndef DATAPERSISTENCE_FILTRATION_H
#define DATAPERSISTENCE_FILTRATION_H

#include "common.h"

#include "CubicalComplex.h"
#include "CubicalCell.h"

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
  std::vector<std::pair<CubicalCell, double>> filtration_;
  std::unordered_map<CubicalCell, uint64_t> cell_indexing_;
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
  ///     valuator  : a function which takes a CubicalCell and returns a value
  ///     direction : either "ascending" or "descending" (gives desired ordering of values)
  Filtration ( CubicalComplex complex,
               std::function<double(CubicalCell const&)> valuator,
               std::string const& direction ) {
    assign(complex, valuator, direction);
  }

  // Filtration
  //   constructor
  void
  assign ( CubicalComplex complex_in,
           std::function<double(CubicalCell const&)> valuator,
           std::string const& direction ) {
    data_ = std::make_shared<Filtration_>();
    complex_() = complex_in;
    // Initialize the filtration
    typedef std::pair<CubicalCell, double> CubicalCellValue;
    for ( auto cell : complex() ) filtration_() . push_back ( { cell, valuator(cell) } );
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
    std::sort(filtration_().begin(), filtration_().end(), comparator);
    // Index the cells
    for ( uint64_t i = 0; i < complex().size(); ++ i ) cell_indexing_()[cell(i)] = i; 
  }

  /// complex
  ///   Return the complex the filtration is associated to
  CubicalComplex const&
  complex ( void ) const {
    return data_ -> complex_;
  }

  /// cell
  ///   Return the ith cell in the filtration
  CubicalCell const&
  cell ( uint64_t i ) const {
    // if ( i >= complex().size() ) {
    //   std::cout << "Filtration::cell(" << i << ")\n";
    //   throw std::out_of_range("Filtration::cell index out of bounds ");
    // }
    return filtration_()[i].first;
  }

  /// index
  ///   Given a cell, return its position in the filtration
  uint64_t
  index ( CubicalCell const& cell) const {
    return cell_indexing_().find(cell) -> second;
  } 

  /// value
  ///   Return the value associated with a cell in the filtration
  double 
  value ( uint64_t i ) const {
    // if ( i >= complex().size() ) {
    //   std::cout << "Filtration::value(" << i << ")\n";
    //   throw std::out_of_range("Filtration::value index out of bounds ");
    // }
    return filtration_()[i].second;
  }

private:

  std::shared_ptr<Filtration_> data_;

  /// complex
  ///   Return the complex the filtration is associated to
  CubicalComplex &
  complex_ ( void ) {
    return data_ -> complex_;
  }

  /// filtration
  ///   accessor for data_ -> filtation_;
  std::vector<std::pair<CubicalCell, double>> &
  filtration_ ( void ) const {
    return data_ -> filtration_;
  }

  /// cell_indexing
  ///   accessor for data_ -> filtation_;
  std::unordered_map<CubicalCell, uint64_t> &
  cell_indexing_ ( void ) const {
    return data_ -> cell_indexing_;
  }

};

#endif
