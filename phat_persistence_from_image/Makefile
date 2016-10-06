PROGRAM = 2d_pic_to_bd_matrix
#PROGRAM = preprocess_image

SRC =

export PREFIX = ../phat
export CXX = g++-4.9
export CXXOPTS = -O2 -std=c++11  

CXXINC = -I$(PREFIX)/include -I../
LIBPATHS = -L$(PREFIX) -L/usr/X11R6/lib 
CXXFLAGS = $(CXXINC) $(CXXOPTS)
LIBRARIES = -lm -lpthread -lX11

.SUFFIXES: .c .cpp .hpp

.c.o:
	$(CXX) $(CXXFLAGS) -c $<

OBJ = $(addsuffix .o, $(basename $(SRC))) $(addsuffix .o, $(PROGRAM))

default: $(PROGRAM)

$(PROGRAM): $(OBJ) $(PROGRAM).cpp
	$(CXX) $(LIBPATHS) -o $(PROGRAM) $(OBJ) $(LIBRARIES)
	chmod 700 $(PROGRAM)

clean:
	rm *.o *~

	
