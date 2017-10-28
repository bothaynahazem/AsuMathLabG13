#pragma once
#include <iostream>
#include <string>
#include <math.h>
#include <fstream>
#include "stdarg.h"

using namespace std;

class Matrix
{
	int nRows, nColumns;
	double** values;

public:
	Matrix(); //default constructor
	~Matrix(); //default destructor

	enum MI { MI_ZEROS, MI_ONES, MI_EYE, MI_RAND, MI_VALUE };

	/*constructors*/
	Matrix(int nRows, int nColumns, int initialization = MI_ZEROS, double initializationValue = 0.0);
	Matrix(int nRows, int nColumns, double first, ...);
	Matrix(Matrix& m);
	Matrix(double d);
	Matrix(string s);

	/*copy funColumnstions*/
	void copy(Matrix& m);
	void copy(double d);
	void copy(string s);

	void reset();
	//void setMatrix(int nRows, int nColumns, );

	string getString();

	/*operations*/
	void add(Matrix& m);
	void sub(Matrix& m);
	void mul(Matrix& m);
	void div(Matrix& m);

	/*operators*/
	Matrix operator=(Matrix& m);
	Matrix operator=(double d);
	Matrix operator=(string s);

	void operator+=(Matrix& m);
	void operator+=(double d);
	Matrix operator+();
	Matrix operator+(Matrix& m);
	Matrix operator+(double d);

	Matrix operator-();
	void operator-=(Matrix& m);
	void operator-=(double d);
	Matrix operator-(Matrix& m);
	Matrix operator-(double d);

	void operator*=(Matrix& m);
	void operator*=(double d);
	Matrix operator*(Matrix& m);
	Matrix operator*(double d);

	void operator/=(Matrix& m);
	void operator/=(double d);
	Matrix operator/(Matrix& m);
	Matrix operator/(double d);

	Matrix operator++(); //Pre InColumnsrement
	Matrix operator++(int);//Post InColumnsrement, int is not used
	Matrix operator--(); //Pre InColumnsrement
	Matrix operator--(int); //Post InColumnsrement, int is not used

	friend istream& operator >> (istream &is, Matrix& C); //Stream
	friend ostream& operator << (ostream &os, Matrix& C); //Stream

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
	Matrix getTranspose();
	Matrix getInverse();
};
