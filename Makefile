TARGETNAME = bbhsim
SRCS     = $(wildcard *.cpp)
OBJS     = $(SRCS:%.cpp=%.o)
OPTIMIZATION = -O3
CXXFLAGS = -std=c++17 -g $(OPTIMIZATION)
LIBS     = -lm -lgsl -lcblas -lconfig++

ifeq ($(HAS_HDF5), yes)
	CXXFLAGS += -DHAS_HDF5
	LIBS += -lhdf5
endif


all: $(OBJS)
	$(CXX) $^ $(LIBS) -o $(TARGETNAME)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ -c

clean:
	rm -f *.o $(TARGETNAME)

.PHONY: all clean