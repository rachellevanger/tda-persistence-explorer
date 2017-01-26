/// SmartConstIterator.h
/// Shaun Harker 2017-01-25-1803
/// MIT LICENSE

#ifndef SMARTCONSTITERATOR_H
#define SMARTCONSTITERATOR_H

/// Synoposis.
///   This is a class which for a given type makes an iterator
///   out of an instance of that type and a "successor" function.

// This only works for types T supporting a "clone" function
// The clone method should produce a distinct copy of the object
// The copy need not be completely deep; only deep enough so that
// we can guarantee that "successor" does not mutate anything which
// if referred to elsewhere in the program.

template<class T>
struct SmartConstIterator_ {
  T t_;
  std::function<void(T &)> successor_; // does this need shared_ptr wrap?
};

template<class T>
class SmartConstIterator {
public:

  /// SmartConstIterator 
  ///   Default constructor
  SmartConstIterator () {}

  /// SmartConstIterator 
  ///   constructor
  SmartConstIterator ( T const& t, std::function<void(T &)> const& successor ) {
    assign(t,successor);
  }
  
  /// assign
  ///   constructor
  void
  assign ( T const& t, std::function<void(T &)> const& successor ) {
    data_ -> t_ = t;
    data_ -> successor_ = successor;
  }

  /// operator *
  ///   dereference iterator
  T const&
  operator * ( void ) const {
    return data_ -> t_;
  }

  /// operator ==
  ///   note: only compares data, not iteration pattern
  bool
  operator == ( SmartConstIterator<T> const& rhs ) const {
    return data_ -> t_ == rhs.data_ -> t_;
  }

  /// operator !=
  ///   note: only compares data, not iteration pattern
  bool
  operator != ( SmartConstIterator<T> const& rhs ) const {
    return not (*this == rhs);
  }

  /// operator ++
  SmartConstIterator<T>
  operator ++ ( void ) {
    if ( data_ . use_count() > 1 ) {
      // Must make a copy before in-place modification
      // or will invalidate other references
      auto newdata_ = std::make_shared<SmartConstIterator_<T>>();
      newdata_ -> t_ = data_ -> t_ . clone ();
      newdata_ -> successor_ = data_ -> successor_;
      data_ = newdata_;
    } 
    // In-place modification
    data_ -> successor_ ( data_ -> t_ );
    return *this;
  }

private:
  std::shared_ptr<SmartConstIterator_<T>> data_;
};

#endif
