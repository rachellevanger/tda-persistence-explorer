# README #

SOURCE CODE ORIGINALLY FROM: https://bitbucket.org/dzako/phat_persistence_from_image

Original code by Jacek Cyranka 

The code assumes that phat library is located in `phat` folder in the directory above, where this program is compiled, i.e.
this is parameter `../phat` in the makefile.

1. PROGRAM '2d_pic_to_bd_matrix.cpp'

computes a cubical complex using PHAT library for persistence, it builds a complex with a prescribed ordering of cells,
computes the boundary matrix of the resulting complex, and computes its persistence using PHAT, and ouputs H0 and H1 
**nondiagonal** generators to files. The parameters are 

`./2d_pic_to_bd_matrix {data_in file name} {output file name} {sub|super} {0,1,2}`, 
for example `./2d_pic_to_bd_matrix data.in out sub 0`, where 

 * `data.in` - the input file with vector of pixel values
 * `out` - the output filename , `out_H0` will contain H0 generators, `out_H1` will contain H1 generators
 * `{sub|super}` - select if *sub level* filtered complex or *super lever* filtered complex is used for homology computation. The complex built for sub level is the *'regular' complex*, whereas the complex built for super level is the *dualized complex*
 * `{0,1,2}` - either `0`, `1` or `2` integer value defining the output format: for `0` only values of birth/death are outputted,
for `1` only x and y coordinates of the generating cells are outputted, having format `birth_x birth_y death_x death_y`,
for `2` full information is outputted in format `value(dimension,birth_x,birth_y) value(dimension,death_x,death_y)`	.
