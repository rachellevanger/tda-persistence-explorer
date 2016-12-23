# ImagePersistence

Jacek Cyranka, Shaun Harker

ImagePersistence.cpp is a C++ program which uses the image library CImg to load 2D images, create sublevel of superlevel filtered complexes, and compute persistence using PHAT. It saves the output to file. 

Installation is accomplished by typing

    ./install.sh


This creates a command-line program `ImagePersistence` with the following usage:

```
  Usage: ImagePersistence <image_filename> <output_filename> <mode>
   <image_filename>  : input image filename
   <output_filename> : output filename (a CSV file)
   <mode>            : A mode of execution (either "sub" or "super") 
                       to indicate whether to compute
                       sublevel or superlevel persistence.
```

## Implementation Details

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



