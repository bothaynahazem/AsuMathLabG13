
compiler=g++  #if we want to change it anytime
flags = -c -Wall  #compiler flags



Matrix2017: main.cpp Matrix2017.cpp
	$(compiler) main.cpp Matrix2017.cpp -o Matrix

clean:
	rm -rf *o
