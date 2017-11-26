
compiler=g++  #if we want to change it anytime
flags = -c -Wall  #compiler flags



all: main.cpp Matrix2017.cpp
	$(compiler) main.cpp Matrix2017.cpp -o Matrix

clean:
	rm -rf *.o

debug: main.cpp Matrix2017.cpp
	$(compiler) -g main.cpp Matrix2017.cpp -o debug
	gdb ./debug

profile: main.cpp Matrix2017.cpp
	$(compiler) -pg main.cpp Matrix2017.cpp -o profile
	gprof ./profile
