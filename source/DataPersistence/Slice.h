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
    address = -1;
  }
  // parameter constructor gives "begin"
  SliceIterator ( std::vector<uint64_t> const& A,
            std::vector<uint64_t> const& B,
            std::vector<uint64_t> const& C,
            std::vector<uint64_t> const& D) : A(A), B(B), C(C), D(D) {
    n = A.size();
    address = 0;
    X = B;
    G.resize(n); // G for "gap"
    J.resize(n, 1); // J for "jump"
    for ( uint64_t i = 0; i < n; ++ i ) G[i] = D[i] - A[i];
    std::partial_sum (G.begin(), G.end()-1, J.begin()+1, std::multiplies<uint64_t>());
    size_ = 1;
    for ( uint64_t i = 0; i < n; ++ i ) size_ *= C[i] - B[i];
    for ( uint64_t i = 0; i < n; ++ i ) G[i] = J[i]*(C[i] - B[i]);
    for ( uint64_t i = 0; i < n; ++ i ) address += J[i] * (B[i] - A[i]);
    // print_vector(A, "A");
    // print_vector(B, "B");
    // print_vector(C, "C");
    // print_vector(D, "D");
    // print_vector(X, "X");
    // print_vector(G, "G");
    // print_vector(J, "J");
    // std::cout << "address = " << address << "\n";
    // std::cout << "n = " << n << "\n";
    // std::cout << "size = " << size_ << "\n";
  }

  bool
  operator == ( SliceIterator const& rhs ) const {
    return address == rhs.address;
  }

  bool
  operator != ( SliceIterator const& rhs ) const {
    return address != rhs.address;
  }

  SliceIterator &
  operator ++ ( void ) {
    for ( int i = 0; i < n; ++ i ) {
      ++ X[i];
      address += J[i];
      if ( X[i] == C[i] ) { 
        X[i] = B[i];
        address -= G[i];
      } else return *this;
    }
    address = -1;
    return *this;
  }

  uint64_t
  operator * ( void ) const {
    return address;
  }

  uint64_t
  size ( void ) const {
    return size_;
  }
private:
  std::vector<uint64_t> A;
  std::vector<uint64_t> B;
  std::vector<uint64_t> C;
  std::vector<uint64_t> D;
  std::vector<uint64_t> X;
  std::vector<uint64_t> G;
  std::vector<uint64_t> J;
  uint64_t n;
  uint64_t size_;
  uint64_t address;
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
