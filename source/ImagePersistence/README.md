# README #


# SOURCE CODE ORIGINALLY FROM: https://bitbucket.org/dzako/phat_persistence_from_image

Original code by Jacek Cyranka 

The code assumes that phat library is located in `phat` folder in the directory above, where this program is compiled, i.e.
this is parameter `../../phat` in the makefile.

1. PROGRAM '2d_pic_to_bd_matrix.cpp'

computes a cubical complex using PHAT library for persistence, it builds a complex with a prescribed ordering of cells,
computes the boundary matrix of the resulting complex, and computes its persistence using PHAT, and ouputs H0 and H1 
**nondiagonal** generators to files. The parameters are 

`./2d_pic_to_bd_matrix {data_in file name} {output file name} {sub|super}`, 
for example `./2d_pic_to_bd_matrix data.in out sub`, where 

 * `data.in` - the input file with vector of pixel values.
 * `out` - the output filename prefix. `out_{sub|super}_H0` will contain H0 generators, `out_{sub|super}_H1` will contain H1 generators, and `out_{sub|super}_all.csv` will contain both dimensions, including geometric information about the location of the paired critical cells.
 * `{sub|super}` - select if *sub level* filtered complex or *super level* filtered complex is used for homology computation. The complex built for sub level is the *'regular' complex*, whereas the complex built for super level is the *dualized complex*.



## Details

We have put together an interface to the persistent homology computing code "Phat"

We are given an N x M image.

## Sublevel sets.

We create a cubical complex of dimension 2 with N x M 2-cells.

Each cell is assigned a value:

* 2-cells have value from image.
* 1-cells are assigned a value equal to the min of the values on their coboundary
* 0-cells are assigned a value equal to the min of the values on their coboundary

We index the cells so that 
  value(cell1) < value(cells2) -> index(cell1) < index(cell2)

## Superlevel sets.

For superlevel sets we do not only negative the image, but also consider a dualization
of the cubical complex.

We create a cubical complex of dimension 2 with (N-1) x (M-1) 2-cells.
Note that such a cubical complex has N x M vertices, which can be put in 1-1
correspondence with the N x M image.

Each cell is assigned a value:

* 0-cells are given the value corresponding to the _negated_ image pixels
* 1-cells are assigned a value equal to the max of the values on their boundary
* 2-cells are assigned a value equal to the max of the values on their boundary


