/// Slice.h
/// Shaun Harker 2017-01-27-0419
/// MIT LICENSE

#ifndef DATAPERSISTENCE_SLICE_H
#define DATAPERSISTENCE_SLICE_H

#include "common.h"

// TODO: comments explaining what this is.
// main idea: A <= B <= X < C <= D. consider [A, D) a cartesian product. index it 0,1,2,...,N-1.
//            Now regard [B,C) as a subset of this cartesian product. Can you quickly iterate
//            through the indices?
// design goal:
//   std::vector<uint64_t> A,B,C,D;
//   ... define them ...
//   auto slice = Slice(A,B,C,D);
//   for ( auto ambient_index : slice ) { ... } 

class SliceIterator {
public:
  // default constructor gives "end"
  SliceIterator ( void ) {
    address_ = -1;
  }
  // parameter constructor gives "begin"
  SliceIterator ( std::vector<uint64_t> const& A,
            std::vector<uint64_t> const& B,
            std::vector<uint64_t> const& C,
            std::vector<uint64_t> const& D) : A_(A), B_(B), C_(C), D_(D) {
    dim_ = A_.size();
    address_ = 0;
    X_ = B_;
    G_.resize(dim_); // G for "gap"
    J_.resize(dim_, 1); // J for "jump"
    for ( uint64_t i = 0; i < dim_; ++ i ) G_[i] = D_[i] - A_[i];
    std::partial_sum (G_.begin(), G_.end()-1, J_.begin()+1, std::multiplies<uint64_t>());
    size_ = 1;
    for ( uint64_t i = 0; i < dim_; ++ i ) size_ *= C_[i] - B_[i];
    for ( uint64_t i = 0; i < dim_; ++ i ) G_[i] = J_[i]*(C_[i] - B_[i]);
    for ( uint64_t i = 0; i < dim_; ++ i ) address_ += J_[i] * (B_[i] - A_[i]);
    if ( size_ == 0 ) address_ = -1;
  }

  bool
  operator == ( SliceIterator const& rhs ) const {
    return address_ == rhs.address_;
  }

  bool
  operator != ( SliceIterator const& rhs ) const {
    return address_ != rhs.address_;
  }

  SliceIterator &
  operator ++ ( void ) {
    for ( int i = 0; i < dim_; ++ i ) {
      ++ X_[i];
      address_ += J_[i];
      if ( X_[i] == C_[i] ) { 
        X_[i] = B_[i];
        address_ -= G_[i];
      } else return *this;
    }
    address_ = -1;
    return *this;
  }

  uint64_t
  operator * ( void ) const {
    return address_;
  }

  uint64_t
  size ( void ) const {
    return size_;
  }
private:
  std::vector<uint64_t> A_;
  std::vector<uint64_t> B_;
  std::vector<uint64_t> C_;
  std::vector<uint64_t> D_;
  std::vector<uint64_t> X_;
  std::vector<uint64_t> G_;
  std::vector<uint64_t> J_;
  uint64_t dim_;
  uint64_t size_;
  uint64_t address_;
};

class Slice { 
public:
  Slice ( std::vector<uint64_t> const& A,
        std::vector<uint64_t> const& B,
        std::vector<uint64_t> const& C,
        std::vector<uint64_t> const& D) {
    begin_ = SliceIterator(A,B,C,D);
    end_ = SliceIterator();
  }

  SliceIterator
  begin ( void ) const {
    return begin_;
  }

  SliceIterator
  end ( void ) const {
    return end_;
  }

  uint64_t
  size ( void ) const {
    return begin_ . size ();
  }

private:
  SliceIterator begin_;
  SliceIterator end_;
};

#endif
