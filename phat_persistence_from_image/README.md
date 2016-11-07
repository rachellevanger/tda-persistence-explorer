# README #

SOURCE CODE ORIGINALLY FROM: https://bitbucket.org/dzako/phat_persistence_from_image

Original code by Jacek Cyranka 

The code assumes that phat library is located in `phat` folder in the directory above, where this program is compiled, i.e.
this is parameter `../../phat` in the makefile.

1. PROGRAM '2d_pic_to_bd_matrix.cpp'

computes a cubical complex using PHAT library for persistence, it builds a complex with a prescribed ordering of cells,
computes the boundary matrix of the resulting complex, and computes its persistence using PHAT, and ouputs H0 and H1 
**nondiagonal** generators to files. The parameters are 

`./2d_pic_to_bd_matrix {data_in file name} {output file name} {sub|super}`, 
for example `./2d_pic_to_bd_matrix data.in out sub`, where 

 * `data.bmp` - the input image in bmp format.
 * `out` - the output filename prefix. `out_{sub|super}_H0` will contain H0 generators, `out_{sub|super}_H1` will contain H1 generators, and `out_{sub|super}_all.csv` will contain both dimensions, including geometric information about the location of the paired critical cells.
 * `{sub|super}` - select if *sub level* filtered complex or *super level* filtered complex is used for homology computation. The complex built for sub level is the *'regular' complex*, whereas the complex built for super level is the *dualized complex*.

