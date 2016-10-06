/*
 * pointer_test.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: cyranka
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <list>
#include <algorithm>
#include <cstring>
#include <cmath>


class Simplex{
public:
  double value;
  int dim;
  int position; //index [0,...,N] //position in the list of all simplices of this dimension
  int globalPos; //index [0,...,N] //position in the whole list of simplices
  int listOfSimplices; //position in the list of all simplices

  std::list< Simplex* > boundary;
  std::list< Simplex* > coboundary;

  Simplex() : boundary(), coboundary() {
    value = HUGE_VAL; dim = 0; position = 0; globalPos= 0; listOfSimplices = -1;

  }

  Simplex(double v, int d, int p) : boundary(), coboundary(){
    value = v;
    dim = d;
    position = p;
    listOfSimplices = -1;
  }

  //sets value to min from coboundaries
  void evaluate(){
    value = HUGE_VAL;
    for (std::list<Simplex*>::iterator iterator = coboundary.begin(), end = coboundary.end(); iterator != end; ++iterator) {
      if( (**iterator).value < value ){
        value = (**iterator).value;
      }
    }
  }

  bool operator<(const Simplex& s2){
    const Simplex& s1 = *this;
    if(s1.value < s2.value) return true;
    //returns ​true if the first argument is less than (i.e. is ordered before) the second
    if(s1.value > s2.value) return false;
    if(s1.value == s2.value){
      if(s1.dim < s2.dim) return true;
      if(s1.dim > s2.dim) return false;
      if(s1.dim == s2.dim){
        if(s1.position < s2.position) return true;
        if(s1.position > s2.position) return false;

        std::cerr << "s1.value==s2.value, s1.dim==s2.dim, s1.position==s2.position, WHICH SHOULDN'T HAPPEN!\n";
        exit(1);
      }
    }
  }


};





std::ostream& operator<<(std::ostream& os, const Simplex& obj){
  //return os << "(" << obj.value << ", " << obj.dim << ", " << obj.position << ", " << obj.globalPos << ")";
  return os << "(" << obj.value << ", " << obj.dim << ", " << obj.position << ")";
}



int main(){
  std::unique_ptr<Simplex> array[3];

  array[0] = std::unique_ptr<Simplex>( new Simplex(1,0,0) );

  array[1] = std::unique_ptr<Simplex>( new Simplex(2,0,0) );

  array[0]->boundary.push_back(array[1].get());
  array[1]->boundary.push_back(array[0].get());

  std::cout << *array[0] << "\nbd:";
  std::cout << *array[0]->boundary.front() << "\n";

  std::cout << *array[1] << "\nbd:";
  std::cout << *array[1]->boundary.front() << "\n";

  array[0].swap(array[1]);

  std::cout << *array[0] << "\nbd:";
  std::cout << *array[0]->boundary.front() << "\n";

  std::cout << *array[1] << "\nbd:";
  std::cout << *array[1]->boundary.front() << "\n";


}









