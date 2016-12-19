# README #

DataPersistence

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

We are given an N x M image. We want to do "sublevel" and "superlevel" persistence. In both cases we assign pixel data to cells in a particular way and sort the cells according to their values. The filtration is then the nested sequence of subcomplexes given by adding each cell in turn.

The guiding principle is that what we want about the filtration is:

(H) For any k, the first k cells in the filtration give a closed subcomplex. 

In the case of sublevel sets, we start by assigning pixel data to two-cells of a N x M complex and then assign values to lower dimensional cells. Cells are then sorted in ascending order of value. Ties are broken arbitrarily. 

In the superlevel case we start by assigning the pixel data to the vertices of a (N-1) x (M-1) complex and then assign values to higher dimensional cells. Cells are then sorted in descending order of value. Ties are broken arbitrarily. 

In both cases we need a means of propagating values onto cells which give the (H) property.

When making the filtration in order of ascending values, to get (H) we need to ensure the following two equivalent conditions:
(1) the value of a cell is greater than or equal to the maximum value on its boundary cells.
(2) the value of a cell is less than or equal to the minimum value on its coboundary cells.

When making the filtration in order of descending values, to get (H) we need to ensure the following two equivalent conditions:
(3) the value of a cell is greater than or equal to the minimum value on its boundary cells.
(4) the value of a cell is less than or equal to the maximum value on its coboundary cells.

The implementation is intended to do the following:

* Sub-level. Two-cells (squares) are initialized with pixel values. By construction we ensure the value of any other cell is the min of the values on its coboundary cells. The cells are then sorted in order of ascending values to obtain a filtration. By (2), the filtration obeys (H). 
* Super-level. Zero-cells (vertices) are initialized with pixel values. By construction we ensure that the value of any other cell is the min of the values on its boundary. The cells are then sorted in order of descending values to obtain a filtration. By (3), the filtration obeys (H).



