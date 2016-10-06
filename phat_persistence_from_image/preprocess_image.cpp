#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cstring>
#include <cmath>

// PHAT include files
#include "phat/representations/vector_vector.h"
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/twist_reduction.h>
#include <phat/compute_persistence_pairs.h>
#include "CImg.h"

using namespace cimg_library;


bool preprocess(const char* file_name , const char* out, int order = 0 ){

  std::ofstream output_stream( out );
  std::ofstream coordinates( "coords.txt" );
  
  if( output_stream.fail() )
    return false;
  
  /* Load the image from file and save its dimensions */
  CImg<int> image( file_name );

  const int data_size_x_ = image . width( ),
            data_size_y_ = image . height( ) ;
  
   
  output_stream << 2 << "\n";
  output_stream << data_size_x_ << "\n";
  output_stream << data_size_y_ << "\n";

  /* Allocate memory for image data */
  int data[ data_size_x_ ][ data_size_y_ ];
  
  /* Load the data from the first channel of the  image */
  for( int x = 0; x < data_size_x_; ++ x){
    for( int y = 0; y < data_size_y_; ++ y){
      data[ x ][ y ] = image( x, y, 0, 1);      
    }
  }
  

  if(order == 0){

    for(int x  =  data_size_x_ - 1; x  >=  0; -- x){

      for(int y  =  0; y  <  data_size_y_; ++ y){
        output_stream << data[ x ][ y ] << "\n";
        coordinates << x << " " << y << "\n";
      }
    }

  }

  if(order == 1){

    for(int x  =  0; x  < data_size_x_; ++ x){

          for(int y  =  0; y  <  data_size_y_; ++ y){

            output_stream << data[ x ][ y ] << "\n";
            coordinates << x << " " << y << "\n";

          }
        }
    
  }
  
  output_stream.close();
  coordinates.close();

}


int main(int argc, char *argv[]){
  
  preprocess(argv[1], argv[2], atoi(argv[3]) );
  
}
