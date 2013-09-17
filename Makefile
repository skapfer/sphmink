all: sphmink

clean:
	rm -f sphmink *.o

sphmink: sphmink.o readpoly.o
	$(CXX) -D_DEBUG -DDEBUG -Wall -o $@ $^ -lgsl -lgslcblas

.PHONY: all clean
