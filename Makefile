#
#	Makefile
#

DBG := 0
CXX 			:= /usr/bin/g++ 

CFLAGS 			:= -std=c++11 -m64 -Wall -Wno-reorder -fPIC -DIL_STD

ifeq ($(DBG), 1)
	CFLAGS 		+= -g -O0
else
	CFLAGS 		+= -O3 -DNDEBUG
endif

ifeq ($(MORE), 1)
	CFLAGS 		+= -DDBG
endif

LDFLAGS         = -no-pie

SRC 			= $(wildcard *.cpp)
HEADER          = $(wildcard *.h)

OBJ 			= $(SRC:.cpp=.o)

EXECUTABLE		= split-gen

all: $(SRC) $(EXECUTABLE) $(HEADER) $(HEADER_LEFDEF) Makefile

$(EXECUTABLE) : $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) $(LIB) $(LDFLAGS) -o $@

%.o: %.cpp $(HEADER) $(HEADER_LEFDEF) Makefile
	$(CXX) -c $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)
