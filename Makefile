
compiler=g++  #if we want to change it anytime
flags = -c -Wall  #compiler flags



Matrix2017: main.test class.test
	$(compiler) main.o Matrix2017.o -o Matrix2017

main.test: main.cpp
	$(compiler) $(flags) main.cpp

class.test: Matrix2017.cpp
	$(compiler) $(flags) Matrix2017.cpp

clean:
	rm -rf *o 
	
