// Code to compute persistence homology using PHAT.
//Ordered simplices
//
// 2016
//
//

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

class Simplex{
public:
  double value;
  int dim;
  int position; //index [0,...,N] //position in the list of all simplices of this dimension
  int globalPos; //index [0,...,N] //position in the whole list of simplices
  int listOfSimplices; //position in the list of all simplices
  double xCoord;
  double yCoord;
  
  std::list< Simplex* > boundary;
  std::list< Simplex* > coboundary;
  
  Simplex() : boundary(), coboundary() {
    value = HUGE_VAL; dim = 0; position = 0; globalPos= 0; listOfSimplices = -1;
    xCoord = 0; yCoord = 0;
  }
  
  Simplex(double v, int d, int p) : boundary(), coboundary(){
    value = v;
    dim = d;
    position = p;
    listOfSimplices = -1; globalPos = -1; xCoord = -1; yCoord = -1;
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
  

  //sets value to max of boundaries
    void dual_evaluate(){
      value = -1;
      for (std::list<Simplex*>::iterator iterator = boundary.begin(), end = boundary.end(); iterator != end; ++iterator) {
        if( (**iterator).value > value ){
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




class Pixel{
public:
  int val;
  int xCoord;
  int yCoord;
  Pixel(int v, int xc, int yc){
    val = v; xCoord = xc; yCoord = yc;
  }
};



enum Ss{ sub , super };



typedef std::unique_ptr<Simplex> pSimplex;

bool compare_ptrs(const pSimplex &p1, const pSimplex &p2){
  return *p1 < *p2;
}

           
    
std::ostream& operator<<(std::ostream& os, const Simplex& obj){
  return os << "(" << obj.value << ", " << obj.dim << ", " << obj.position << ", " << obj.globalPos << ", " << obj.xCoord << ", " << obj.yCoord << ")";
  //return os << "(" << obj.value << ", " << obj.dim << ", " << obj.position << ")";
}


std::vector<int> read_data(std::vector<int>& xcoords, std::vector<int>& ycoords, std::string input_file_name, std::string output_file_name_prefix, int data_sizes[]){
  bool verbose = true;  
  
  std::cout << "Input file name: " << input_file_name << std::endl;
  std::cout << "Output file name prefix: " << output_file_name_prefix << std::endl;
  

  const char *input_fname = input_file_name . c_str ();

   std::ifstream sol_file(input_fname);

   FILE * coords_file = fopen("coords.txt", "r");

   // Make sure the file was opened correctly
   if (!sol_file . is_open()) {
     std::cout << "Unable to open file: " << input_fname << std::endl;
     exit (1);
   }

   int d;
   sol_file >> d;
   if(d != 2){
     std::cerr << "Working only for 2D!\n";
     exit(1);
   }
   sol_file >> data_sizes[0];
   std::cout << "Read x_dim=" << data_sizes[0] << "\n";
   sol_file >> data_sizes[1];
   std::cout << "Read y_dim=" << data_sizes[1] << "\n";
   
   
   // Create a vector of {value, index} pairs
   std::vector< int > v_sol_indx;
   double v = 0.0;

   if (verbose) {
     std::cout << std::endl;
     std::cout << "Reading data ...";
   }

   int indx = 0, d1, d2;

   // Read the data from file
   while (sol_file >> v) {
     fscanf( coords_file, "%d %d", &d1, &d2 );
     xcoords.push_back(d1);
     ycoords.push_back(d2);
     v_sol_indx . push_back( v );  // Add pair
     ++indx;                                             // Increment index
   }

   sol_file . close();
   fclose( coords_file );

   if (verbose) {
     std::cout << " done!" << std::endl;
   }

   // Number of values read from file
   int num_values = v_sol_indx . size();

   // Check if size of data matches data read from file
   int total_data_size = data_sizes [0] * data_sizes [1];

   if (num_values != total_data_size) {
     std::cout << "Data sizes do not match! " << num_values << " " << total_data_size << std::endl;
     exit (1);
   }

   if (verbose) {
     std::cout << "Grid size: " << data_sizes [0] << " x " << data_sizes [1] << std::endl;
   }
   
   return v_sol_indx;
}




//adds simplices vertices
void generateVertices(std::vector<pSimplex>& out, std::list<Pixel>& pixels, int data_sizes[], int suborsuper){

  int i, j;

  if(suborsuper == sub){
    const int N = (data_sizes[0] + 1) * (data_sizes[1] + 1) ;
    for( i = 0; i < N ; i ++){
      out[i] = pSimplex( new Simplex(HUGE_VAL, 0, i) );
    }
  }
  if(suborsuper == super){
    const int N = (data_sizes[0]) * (data_sizes[1]) ;
    std::list<Pixel>::reverse_iterator iterator, end;

    for( i = 0, iterator = pixels.rbegin(), end = pixels.rend() ; i < N, iterator != end ; i++, ++iterator){
      out[i] = pSimplex( new Simplex( 255-(*iterator).val, 0, i) );
      out[i]->xCoord = (*iterator).xCoord;
      out[i]->yCoord = (*iterator).yCoord;
    }
  }
}



//adds simplices edges horizontal
void generateEdgesH(std::vector<pSimplex>& simplices, int Nvert, int data_sizes[], int &index, int suborsuper){
  int i, j;
  int N, SIZE;
  if(suborsuper == sub){
    N = (data_sizes[0] + 1) * (data_sizes[1] + 1) - 1;
    SIZE = data_sizes[0];
  }
  if(suborsuper == super){
    N = (data_sizes[0]) * (data_sizes[1]) - 1;
    SIZE = data_sizes[0] - 1;
  }
    
  int array_index = Nvert; // indexes simplices in the corresponding array of simplices 
  
  //adds horizontal edges first  
  for( i = 0; i < N ; i ++){
    if( i == 0 || i % (SIZE + 1) != SIZE ){
      Simplex* s = new Simplex(HUGE_VAL, 1, index++) ;
  
      s->boundary .push_back( simplices[i].get() );
      s->boundary .push_back( simplices[i+1].get() );

      
      simplices[array_index] = pSimplex(s);           
            
      simplices[i]->coboundary.push_back( simplices[array_index].get() );
      simplices[i+1]->coboundary.push_back( simplices[array_index].get() );
      
      array_index++;
    }
  }  
}


//adds simplices edges vertical
void generateEdgesV(std::vector<pSimplex>& simplices, int NvertNedgesH, int data_sizes[], int &index, int suborsuper){
  int i, j;
  int N, SIZE;
  if(suborsuper == sub){
    N = (data_sizes[0] ) * (data_sizes[1] + 1);
    SIZE = data_sizes[0];
  }
  if(suborsuper == super){
    N = (data_sizes[0] - 1) * (data_sizes[1] );
    SIZE = data_sizes[0] - 1;
  }


  int array_index = NvertNedgesH;
 
  //adds horizontal edges first  
  for( i = 0; i < N ; i ++){
    Simplex* s = new Simplex(HUGE_VAL, 1, index++);
    
    s->boundary .push_back( simplices[i].get() );
    s->boundary .push_back( simplices[ i + SIZE + 1 ].get() );


    simplices[array_index] = pSimplex(s);        
    
    simplices[ i ]->coboundary.push_back( simplices[array_index].get() );
    simplices[ i + SIZE + 1 ]->coboundary.push_back( simplices[array_index].get() );
    
    array_index++;
  }
}

//adds cubes 
void generateCubes(std::vector<pSimplex>& simplices, int Nvert, int NedgesH, int NedgesV, std::list<Pixel>& pixels, int data_sizes[], int suborsuper){
  int i,j;
  int N, SIZE;
  if( suborsuper == sub){
    N = pixels.size();
    SIZE = data_sizes[0];
  }
  if(suborsuper == super){
    N = ( data_sizes[0] - 1 ) * ( data_sizes[1] - 1 );
    SIZE = data_sizes[0] - 1;
  }
    
  int index = 0;
  
  std::list<Pixel>::iterator iterator, end;

  for( i = 0, iterator = pixels.begin(), end = pixels.end(); i < N && iterator != end; i++, ++iterator){

    Simplex* s;
    if( suborsuper == sub){
      s = new Simplex( (*iterator).val, 2, i) ;
      s->xCoord = (*iterator).xCoord;
      s->yCoord = (*iterator).yCoord;
    }
    if( suborsuper == super ){
      s = new Simplex( HUGE_VAL, 2, i) ;
    }

    int HedgesSTART, VedgesSTART;
    if(suborsuper == sub){
      HedgesSTART = Nvert;
      VedgesSTART = Nvert + NedgesH;
    }
    if(suborsuper == super){
      HedgesSTART = Nvert + NedgesV;
      VedgesSTART = Nvert;
    }

    s->boundary.push_back( simplices[ HedgesSTART + i ].get() );
    s->boundary.push_back( simplices[ HedgesSTART + i + SIZE ].get() );

    if(i != 0 && i % SIZE == 0) index++;

    s->boundary.push_back( simplices[ VedgesSTART + i + index ].get() );
    s->boundary.push_back( simplices[ VedgesSTART + i + 1 + index ].get() );

    simplices[Nvert + 2*NedgesH + i] = pSimplex( s );
    
    simplices[ HedgesSTART + i ]->coboundary.push_back( simplices[ Nvert + 2*NedgesH + i ].get() );
    simplices[ HedgesSTART + i + SIZE ]->coboundary.push_back( simplices[ Nvert + 2*NedgesH + i ].get() );
    
    simplices[ VedgesSTART + i + index ]->coboundary.push_back( simplices[ Nvert + 2*NedgesH + i ].get() );
    simplices[ VedgesSTART + i + 1 + index ]->coboundary.push_back( simplices[ Nvert + 2*NedgesH + i ].get() );
  }
  
}


//evaluates values of simplices
void evaluate(std::vector<pSimplex>& simplices, int Nvert, int Ncubes, int data_sizes[], int suborsuper){
  int N = simplices.size();
  int i;
  if(suborsuper == sub){
    for(i = Nvert; i < N - Ncubes; i++){
      simplices[i]->evaluate();
    }

    for(i = 0; i < Nvert; i++){
      simplices[i]->evaluate();
    }
  }
  if(suborsuper == super){
    for(i = Nvert; i < N; i++){
      simplices[i]->dual_evaluate();
    }

  }

}

//this function sets the x and y coordinates of the cells,
//it should be updated such that coordinates of eges, and vertices have not integer values, but halves
void updateCoords(int Nvert, int NedgesH, int NedgesV, std::vector<pSimplex>& simplices , int suborsuper){

  int i, j;
  const int startC = Nvert + NedgesH + NedgesV;

  if( suborsuper == sub ){
    //update coordinates for edges (boundaries of cubes)
    for(i = startC ; i < simplices.size() ; i++){
      j = 0;
      for (std::list<Simplex*>::iterator iterator = simplices[i]->boundary.begin(), end = simplices[i]->boundary.end(); iterator != end; ++iterator, ++j) {
        (*iterator)->xCoord = simplices[i]->xCoord;
        (*iterator)->yCoord = simplices[i]->yCoord;

      }
    }

    //update coordinates for vertices (boundaries of horizontal edges)
    int startVH = Nvert;
    int startVV = Nvert + NedgesH;

    for( i = startVH ; i < startVV ; i++ ){
      j = 0;
      for (std::list<Simplex*>::iterator iterator = simplices[i]->boundary.begin(), end = simplices[i]->boundary.end(); iterator != end; ++iterator, ++j) {
        (*iterator)->xCoord = simplices[i]->xCoord;
        (*iterator)->yCoord = simplices[i]->yCoord;

      }
    }

    //update coordinates for vertices (boundaries of vertical edges)
    for(i = startVV ; i < startC ; i++){
      j = 0;
      for (std::list<Simplex*>::iterator iterator = simplices[i]->boundary.begin(), end = simplices[i]->boundary.end(); iterator != end; ++iterator, ++j) {
        (*iterator)->xCoord = simplices[i]->xCoord;
        (*iterator)->yCoord = simplices[i]->yCoord;

      }
    }
  }

  if( suborsuper == super ){
    double x, y;
    for( i = Nvert ; i < simplices.size(); i++ ){
      x = 0; y = 0;
      for (std::list<Simplex*>::iterator iterator = simplices[i]->boundary.begin(), end = simplices[i]->boundary.end(); iterator != end; ++iterator, ++j) {
        x += (*iterator)->xCoord; y += (*iterator)->yCoord;
      }
      simplices[i]->xCoord = x / simplices[i]->boundary.size();
      simplices[i]->yCoord = y / simplices[i]->boundary.size();
    }
  }

}


int main(int argc, char *argv[])
{
  std::string in(argv[1], strlen(argv [1]));
  std::string out(argv[2], strlen(argv [2]));

  std::string outAll(argv[2], strlen(argv [2]));
  std::string outH0(argv[2], strlen(argv [2]));
  std::string outH1(argv[2], strlen(argv [2]));

  std::string subsuper(argv[3], strlen(argv[3]));
    
  // int output_type = atoi( argv[4] );



  int suborsuper;

  if( std::strcmp( subsuper.c_str(), "sub" ) == 0){
    suborsuper = sub;
  }else{
    if( std::strcmp( subsuper.c_str(), "super" ) == 0){
      suborsuper = super;
    }else{
      std::cerr << "LAST PARAMETER IS EITHER 'sub' OR 'super'\n";
      throw std::runtime_error("LAST PARAMETER IS EITHER 'sub' OR 'super'\n");
    }
  }

  
  outH0.append("_").append(subsuper).append("_H0");
  outH1.append("_").append(subsuper).append("_H1");
  outAll.append("_").append(subsuper).append("_all.csv");

  int data_sizes[2];

	/* Load the image from file and save its dimensions */
  CImg<int> image( argv[1] );

  const int data_size_x_ = image . width( ),
            data_size_y_ = image . height( ) ;

  data_sizes[0] = data_size_x_;
  data_sizes[1] = data_size_y_;

  //build a double-linked list with pixel values
  std::list<Pixel> pixels;

  //for( int x = data_size_x_ - 1; x >= 0 ; -- x){
  for( int x = 0; x < data_size_x_ ; ++ x){
    for( int y = data_size_y_ - 1; y >= 0; --y){
      pixels.push_back( Pixel( image( x, y, 0, 1), x, y ) );
    }
  }

  std::vector<pSimplex> simplices;

  if( suborsuper == sub ){
    int Nvert = (data_sizes[0] + 1) * (data_sizes[1] + 1);
    int NedgesH = ( (data_sizes[0] ) * (data_sizes[1] + 1) );
    int NedgesV = ( (data_sizes[0] ) * (data_sizes[1] + 1) );
    int Ncubes = pixels.size();

    simplices.resize( Nvert + NedgesH + NedgesV + Ncubes );
  
    generateVertices(simplices, pixels, data_sizes, suborsuper);
    int index = 0;

    generateEdgesH(simplices, Nvert, data_sizes, index, suborsuper);
    generateEdgesV(simplices, Nvert + NedgesH, data_sizes, index, suborsuper);
    generateCubes( simplices, Nvert, NedgesH, NedgesV, pixels, data_sizes, suborsuper);

    evaluate(simplices, Nvert, Ncubes, data_sizes, suborsuper);
  
    updateCoords(Nvert, NedgesH, NedgesV, simplices, suborsuper);
  }
  
  if( suborsuper == super ){

    int Nvert = (data_sizes[0]) * (data_sizes[1]);
    int NedgesH = ( (data_sizes[0] - 1 ) * ( data_sizes[1] ) );
    int NedgesV = ( (data_sizes[0] - 1 ) * ( data_sizes[1] ) );
    int Ncubes = ( data_sizes[0] - 1 ) * ( data_sizes[1] - 1 );

    simplices.resize( Nvert + NedgesH + NedgesV + Ncubes );

    generateVertices(simplices, pixels, data_sizes, suborsuper);
    int index = 0;

    generateEdgesV(simplices, Nvert , data_sizes, index, suborsuper);
    generateEdgesH(simplices, Nvert + NedgesV, data_sizes, index, suborsuper);
    generateCubes( simplices, Nvert, NedgesH, NedgesV, pixels, data_sizes, suborsuper);

    evaluate(simplices, Nvert, Ncubes, data_sizes, suborsuper);

    updateCoords(Nvert, NedgesH, NedgesV, simplices, suborsuper);
  }


  std::sort( simplices.begin(), simplices.end(), [](const pSimplex &p1, const pSimplex &p2){
    return *p1 < *p2;
  } );


  for(int i=0; i < simplices.size(); i++){
    simplices[i]->globalPos = i;
  }
  

  //DEBUGING CODE START

  /*std::cout << "simplices:\n";
  for(int i=0; i < simplices.size(); i++)
  {
    std::cout << *(simplices[i]) ;

    std::list<Simplex*>::iterator iterator, end;
    std::cout << "\n bd: ";
    for (iterator = simplices[i]->boundary.begin(), end = simplices[i]->boundary.end(); iterator != end; ++iterator) {
        std::cout << **iterator << ", ";
    }

    std::cout << "\n cobd:";
    for (iterator = simplices[i]->coboundary.begin(), end = simplices[i]->coboundary.end(); iterator != end; ++iterator) {
        std::cout << **iterator << ", ";
    }
    std::cout << "\n\n";
  }*/

  //DEBUGING CODE END



  //PHAT CODE START

  //build a boundary matrix
  phat::boundary_matrix< phat::vector_vector > boundary_matrix;
  const int num_cols = simplices.size();
  boundary_matrix . set_num_cols(num_cols);
  std::vector< phat::index > col_vector;
  int col_indx = 0;
  
  // Add cells to the boundary matrix (filtration)
  for(int k0 = 0; k0 < num_cols; ++k0) {
    std::vector< phat::index > col_vector;
    col_vector . clear();    
    
    for (std::list<Simplex*>::iterator iterator = simplices[k0]->boundary.begin(), end = simplices[k0]->boundary.end(); iterator != end; ++iterator) {
      col_vector.push_back((*iterator)->globalPos);
    }        
    
    std::sort( col_vector.begin(), col_vector.end() );
    
    // Set current column (add empty column)
    const int col_dim = simplices[col_indx]->dim;
    boundary_matrix . set_dim( col_indx, col_dim);    
    boundary_matrix . set_col( col_indx++, col_vector );

  }  
  
  //boundary_matrix.save_ascii( out );
  
  phat::persistence_pairs pairs;
  
  phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );
  
  int NUM_PAIRS = pairs.get_num_pairs();

    std::ofstream outputH0 ( outH0 );
    std::ofstream outputH1 ( outH1 );

    std::ofstream outputAll ( outAll );
    outputAll << "dim,birth,b_x,b_y,b_z,death,d_x,d_y,d_z\n" ;


    for( int i =0; i < NUM_PAIRS; i++){
      if( simplices[pairs.get_pair(i).first]->value != simplices[pairs.get_pair(i).second]->value ){

        if( (int)boundary_matrix . get_dim( pairs.get_pair(i).first ) == 0 ) {
          //differientiaties the output depending if sub / super
          if(suborsuper == sub)
            outputH0 <<  simplices[pairs.get_pair(i).first]->value << " " << simplices[pairs.get_pair(i).second]->value << "\n";
          if(suborsuper == super)
            outputH0 <<  255 - simplices[pairs.get_pair(i).first]->value << " " << 255 - simplices[pairs.get_pair(i).second]->value << "\n";
  
        }

        if( (int)boundary_matrix . get_dim( pairs.get_pair(i).first ) == 1 ){
          //differientiaties the output depending if sub / super
          if(suborsuper == sub)
            outputH1 <<  simplices[pairs.get_pair(i).first]->value << " " << simplices[pairs.get_pair(i).second]->value << "\n";
          if(suborsuper == super)
            outputH1 <<  255 - simplices[pairs.get_pair(i).first]->value << " " << 255 - simplices[pairs.get_pair(i).second]->value << "\n";

        }

        if(suborsuper == sub)
          outputAll << simplices[pairs.get_pair(i).first]->dim << ", " << simplices[pairs.get_pair(i).first]->value << ", " << simplices[pairs.get_pair(i).first]->xCoord << ", " << simplices[pairs.get_pair(i).first]->yCoord << ", 0, "
          << simplices[pairs.get_pair(i).second]->value << ", " << simplices[pairs.get_pair(i).second]->xCoord << ", " << simplices[pairs.get_pair(i).second]->yCoord << ", 0" << "\n";
        if(suborsuper == super)
          outputAll << simplices[pairs.get_pair(i).first]->dim << ", " <<  255 - simplices[pairs.get_pair(i).first]->value << ", " << simplices[pairs.get_pair(i).first]->xCoord << ", " << simplices[pairs.get_pair(i).first]->yCoord << ", 0, "
          << 255 - simplices[pairs.get_pair(i).second]->value << ", " << simplices[pairs.get_pair(i).second]->xCoord << ", " << simplices[pairs.get_pair(i).second]->yCoord << ", 0" << "\n";

      }
    }



  //PHAT CODE END

    outputH0.close();
    outputH1.close();
    outputAll.close();


}
