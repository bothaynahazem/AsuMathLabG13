#include <iostream>
#include <string>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string.h>
#include <stdio.h>
#include "stdarg.h"
#include <cmath>
#include <complex>
using namespace std;

int LUPDecompose(double **A, int N, double Tol, int *P);
double LUPDeterminant(double **A, int *P, int N);
int LUPDecompose(complex<double> **A, int N, double Tol, int *P);
complex<double> LUPDeterminant(complex<double> **A, int *P, int N);

class Matrix
{
	int nRows, nColumns;
	double** values;
	complex<double>** cvalues;
	string type;


public:
	Matrix(); //default constructor
	~Matrix(); //default destructor

	enum MI { MI_ZEROS, MI_ONES, MI_EYE, MI_RAND, MI_VALUE };

	/*constructors*/
	Matrix(int nRows, int nColumns, int initialization = MI_ZEROS, double initializationValue = 0.0);
	Matrix(string type,int nRows, int nColumns, int initialization=MI_ZEROS, complex<double> initializationValue=0);
	Matrix(int nRows, int nColumns, double first, ...);
	Matrix(Matrix& m);
	Matrix(const Matrix& m);
	Matrix(double d);
	Matrix(complex<double> d);
	Matrix(const string s);
  Matrix(int nRows, int nColumns, complex<double> first, ...);
  Matrix(const string type,const string s);
	/*copy functions*/
	void copy(const Matrix& m);
	void copy(const double d);
  void copy(complex<double> d);
	void copy(const string s);
  void copy(const string type,const string s);
	void reset();

	string getString();
  string getAltString();
  string getAltStringNoB();
	/*operations*/
	void add(const Matrix& m);
	void sub(const Matrix& m);
	/*
	it is used to:
	subtract a matrix elements' values from another matrix
	example:
	Matrix a,b;
	a.sub(b);
	*/
	void mul(Matrix& m);
	/*
	it is used to:
	multiply (this) Matrix to another Matrix (m) then saves the result in (this) Matrix
	so it doesn't return anything
	example:
	Matrix a,b;
	a.mul(b);
	*/

	/*operators*/
	Matrix operator=(const Matrix& m);
	Matrix operator=(const double d);
	Matrix operator=(const complex<double> d);
	Matrix operator=(const string s);

	void operator+=(Matrix& m);
	void operator+=(double d);
        void operator+=(complex<double> d);
	Matrix operator+();
	Matrix operator+(Matrix& m);
	/*
	it is used to:
       add two matrices and return a matrix
       example:
       Matrix a,b,c;
       a=b+c;
       //or a=a+b;
	*/

	Matrix operator+(double d);
	Matrix operator+(complex<double> d);
	/*
    it used to:
    add a value to each element of the matrix the called the operator then outputs it
    example:
    double d;
    Matrix a,b;
    a=a+d;
    //or
    a=b+d;
    */

	Matrix operator-();
	/*
	it is used to:
	to multiply each element of (this) Matrix with -1
	it returns a Matrix just to allow this line (a=-b;)
	but you can use it like this (-a;)

	example:
	Matrix a,b;
	a=-b;
	//or
	-a;
	*/
	void operator-=(Matrix& m);
	/*
	it is used to:
	subtract matrix elements (m) from another matrix elements (this matrix)//the obj that called the operator
	example:
	Matrix a,b;
	a-=b;
	*/
	void operator-=(double d);
	void operator-=(complex<double> d);
	/*
	it is used to:
	subtract a double (d) from the matrix elements (this)
	example:
	Matrix a;
	double d;
	a-=d;
	*/
	Matrix operator-(Matrix& m);
	/*
	it is used to:
	subtract two matrices,, and return the result matrix
	example:
	Matrix a,b,c;
	a=b-c;
	//or
	a=a-b;
	*/
	Matrix operator-(double d);
	Matrix operator-(complex<double> d);
	/*
	it is used to:
	subtract double from matrix,, and return the result matrix
	example:
	double d;
	Matrix a,b;
	b=a-d;
	*/

	void operator*=(Matrix& m);
	/*
	it is used to:
	multiply (this) Matrix to a Matrix (m)
	and then saves the result in (this) Matrix
	so it doesn't return anything
	example:
	Matrix a,b;
	a*=b;
	*/
	void operator*=(double d);
	void operator*=(complex<double> d);
	/*
	it is used to:
	multiply each element of (this) Matrix by double d then saves the result in (this) Matrix
	so it doesn't return anything
example:
	double d;
	Matrix a;
	a*=d;
	*/
	Matrix operator*(Matrix& m);
	/*
	it is used to:
	multiply (this) Matrix by another Matrix (m) and returns the result
example:
	Matrix a,b,c;
	a=b*c;
	//or
	a=a*b;
	*/
	Matrix operator*(double d);
	Matrix operator*(complex<double> d);
	/*
	it is used to:
	multiply elements of (this) Matrix by double (d) and return the result Matrix


	example:
	double d;
	Matrix a,b;
	a=b*d;
	//or
	a=a*d;
	*/

	void operator/=(Matrix& m);
	void operator/=(double d);
        void operator/=(complex<double> d);
	Matrix operator/(Matrix& m);
	Matrix operator/(double d);
        Matrix operator/(complex<double> d);


	Matrix operator++(); //Pre Increment
	/*
	it is used to:
	increase the value of elements of (this) Matrix by 1
	it returns Matrix so it can allow this sentence (a=b++;)
	but you can use it like this b++;
	and the returned matrix will be deleted and won't give any runtime error
	example:
	Matrix a,b;
	a++;
	//or
	b=a++;
	*/
	Matrix operator++(int);//Post Increment, int is not used
	/*
		it is used to:
	increase the value of elements of (this) Matrix by int 1
	but it returns the old value
	the int is used to differentiate between it and +++()
	the old value is saved in C
	example:
	int x
	Matrix a,b;
	a++x;
	//or
	b=a++x;
	*/
	Matrix operator--(); //Pre Increment
	/*
	it is used to:
	decrease the value of elements of (this) Matrix by 1
	and return the new value
	example:
	Matrix a,b;
	a--;
	//or
	b=a--;
	*/

	Matrix operator--(int); //Post Increment, int is not used
	/*
	it is used to:
	decrease the value of elements of (this) Matrix by 1
	and returns the old value
	the old value is saved in r
	example:
	Matrix a,b;
	int x;
	a--x;
	//or
	b=a--x;
	*/

	/*ADVANCED OPERATIONS*/
	Matrix rdivide(double d, const Matrix& m);
    /*simple right array division
     --> 1./A <-- for example
     dividing 1 by every element in the array A
    */
	Matrix rdivide(const Matrix& m1, const Matrix& m2);
    /*more complex right array division
     --> A./B <-- for example
     dividing every element in A by every corresponding element in the array B

     from MatLab:
     The numerator input a can be complex,
     but the denominator b requires a real-valued input.
     If a is complex, the real and imaginary parts of a are independently divided by b.
    */

	friend istream& operator >> (istream &is, Matrix& C); //input stream
	friend ostream& operator << (ostream &os, Matrix& C); //output stream

	void setSubMatrix(int iR, int iC, Matrix& m);
	Matrix getSubMatrix(int r, int c, int nRows, int nColumns);
	Matrix getCofactor(int r, int c);

	void addColumn(Matrix& m);
	void addRow(Matrix& m);

	double& operator[](int i)
	{
		return values[i / nColumns][i%nColumns];
	}
	double& operator()(int i)
	{
		return values[i / nColumns][i%nColumns];
	}
	double& operator()(int r, int c)
	{
		return values[r][c];
	}
	int getn()
	{
		return nRows*nColumns;
	};
	int getnRows()
	{
		return nRows;
	};

	int getnColumns()
	{
		return nColumns;
	};

	double getDeterminant();
	double getdDeterminant();
	Matrix getTranspose();
	Matrix getInverse();
	Matrix div(Matrix &m);
    static Matrix sin(Matrix&s);
    static Matrix cos(Matrix&s);
    static Matrix tan(Matrix&s);
    static Matrix sinh(Matrix&s);
    static Matrix cosh(Matrix&s);
    static Matrix tanh(Matrix&s);
    static Matrix asin(Matrix&s);
    static Matrix acos(Matrix&s);
    static Matrix atan(Matrix&s);
    static Matrix a_sinh(Matrix&s);
    static Matrix a_cosh(Matrix&s);
    static Matrix a_tanh(Matrix&s);
    static Matrix log(Matrix&s);
    static Matrix log10(Matrix&s);
    static Matrix pow(Matrix&s,double x);
    static Matrix rpow(Matrix&s,double x);
    static Matrix exp(Matrix&s);
    static Matrix sqrt(Matrix&s);
    static Matrix cubicrt(Matrix&s);
    static Matrix ceil(Matrix&s);
    static Matrix floor(Matrix&s);
    static complex<double> complex_parser(const string cs);
	complex<double> getcDeterminant();
};
