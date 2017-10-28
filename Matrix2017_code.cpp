#include "Matrix2017.h"

int max(int a, int b)
{
	return (a > b ? a : b);
}

/*Matrix class implementation*/
Matrix::Matrix() //default constructor
{
	nRows = nColumns = 0;
	values = NULL;
}
Matrix::~Matrix() //default destructor
{
	reset();
}

Matrix::Matrix(int nRows, int nColumns, int initialization, double initializationValue)
{
	this->nRows = nRows;
	this->nColumns = nColumns;
	if ((nRows*nColumns) == 0)
	{
		values = NULL;
		return;
	}
	values = new double*[nRows];
	for (int iR = 0;iR<nRows;iR++)
	{
		values[iR] = new double[nColumns];
		for (int iC = 0;iC<nColumns;iC++)
		{
			switch (initialization)
			{
			case MI_ZEROS:
				values[iR][iC] = 0;
				break;
			case MI_ONES:
				values[iR][iC] = 1;
				break;
			case MI_EYE:
				values[iR][iC] = (iR == iC) ? 1 : 0;
				break;
			case MI_RAND: values[iR][iC] = (rand() % 1000000) / 1000000.0;
				break;
			case MI_VALUE:
				values[iR][iC] = initializationValue;
				break;
			}
		}
	}
}

Matrix::Matrix(int nRows, int nColumns, double first, ...)
{
	this->nRows = nRows;
	this->nColumns = nColumns;
	if ((nRows*nColumns) == 0)
	{
		values = NULL;
		return;
	}
	values = new double*[nRows];
	va_list va;
	va_start(va, first);
	for (int iR = 0;iR<nRows;iR++)
	{
		values[iR] = new double[nColumns];
		for (int iC = 0;iC<nColumns;iC++)
		{
			values[iR][iC] = (iC == 0 && iR == 0) ? first : va_arg(va, double);
		}
	}
	va_end(va);
}

Matrix::Matrix(Matrix& m)
{
	nRows = nColumns = 0;
	values = NULL;
	copy(m);
}

Matrix::Matrix(string s)
{
	nRows = nColumns = 0;
	values = NULL;
	copy(s);
}

Matrix::Matrix(double d)
{
	nRows = nColumns = 0;
	values = NULL;
	copy(d);
}

void Matrix::copy(Matrix& m)
{
	reset();
	this->nRows = m.nRows;
	this->nColumns = m.nColumns;
	if ((nRows*nColumns) == 0)
	{
		values = NULL;
		return;
	}
	values = new double*[nRows];
	for (int iR = 0;iR<nRows;iR++)
	{
		values[iR] = new double[nColumns];
		for (int iC = 0;iC<nColumns;iC++)
		{
			values[iR][iC] = m.values[iR][iC];
		}
	}
}

void Matrix::copy(double d)
{
	reset();
	this->nRows = 1;
	this->nColumns = 1;
	values = new double*[1];
	values[0] = new double[1];
	values[0][0] = d;
}

void Matrix::copy(string s)
{
	reset();
	char* buffer = new char[s.length() + 1];
	strcpy_s(buffer, s.length()+1, s.c_str());
	char* lineContext;
	char* lineSeparators = ";\r\n";
	char* line = strtok_s(buffer, lineSeparators, &lineContext);
	while (line)
	{
		Matrix row; char* context;char* separators = " []";
		char* token = strtok_s(line, separators, &context);
		while (token)
		{
			Matrix item = atof(token);
			row.addColumn(item);
			token = strtok_s(NULL, separators, &context);
		}
		if (row.nColumns>0 && (row.nColumns == nColumns || nRows == 0))
			addRow(row);line = strtok_s(NULL, lineSeparators, &lineContext);
	}
	delete[] buffer;
}

void Matrix::reset()
{
	if (values)
	{
		for (int i = 0;i<nRows;i++)
			delete[] values[i];
		delete[] values;
	}
	nRows = nColumns = 0;
	values = NULL;
}

string Matrix::getString()
{
	string s;
	for (int iR = 0; iR<nRows; iR++)
	{
		for (int iC = 0; iC<nColumns; iC++)
		{
			cout << values[iR][iC] << " ";

			char buffer[50];
			//sprintf_s(buffer, 50, "%g\t", values[iR][iC]);
			s += buffer;
		}
		cout<<endl;
		s += "\n";
	}
	return s;
}

Matrix Matrix::operator=(Matrix& m)
{
	copy(m);
	return *this;
}

Matrix Matrix::operator=(double d) { copy(d);return *this; }
Matrix Matrix::operator=(string s) { copy(s);return *this; }

void Matrix::add(Matrix& m)
{
	if (nRows != m.nRows || nColumns != m.nColumns)
		throw("Invalid matrix dimension");
	for (int iR = 0; iR < nRows; iR++)
	{
		for (int iC = 0;iC < nColumns;iC++)
			values[iR][iC] += m.values[iR][iC];
	}
}
void Matrix::operator+=(Matrix& m) { add(m); }
void Matrix::operator+=(double d) { add(Matrix(nRows, nColumns, MI_VALUE, d)); }
Matrix Matrix::operator+(Matrix& m) { Matrix r = *this;r += m;return r; }
Matrix Matrix::operator+(double d) { Matrix r = *this;r += d;return r; }

void Matrix::sub(Matrix& m)
{
	if (nRows != m.nRows || nColumns != m.nColumns)
		throw("Invalid matrix dimension");
	for (int iR = 0;iR < nRows;iR++)
	{
		for (int iC = 0;iC < nColumns;iC++)
			values[iR][iC] -= m.values[iR][iC];
	}
}

void Matrix::operator-=(Matrix& m) { sub(m); }
void Matrix::operator-=(double d) { sub(Matrix(nRows, nColumns, MI_VALUE, d)); }
Matrix Matrix::operator-(Matrix& m) { Matrix r = *this;r -= m;return r; }
Matrix Matrix::operator-(double d) { Matrix r = *this;r -= d;return r; }

void Matrix::mul(Matrix& m)
{
	if (nRows != m.nRows || nColumns != m.nColumns)
		throw("Invalid matrix dimension");
	Matrix r(nRows, m.nColumns);
	for (int iR = 0;iR<r.nRows;iR++)
	{
		for (int iC = 0;iC < r.nColumns;iC++)
		{
			r.values[iR][iC] = 0;
			for (int k = 0;k < m.nColumns;k++)
				r.values[iR][iC] += values[iR][k] * m.values[k][iC];
		}
	}
	copy(r);
}

void Matrix::operator*=(Matrix& m) { mul(m); }
void Matrix::operator*=(double d)
{
	for (int iR = 0;iR<nRows;iR++)
		for (int iC = 0;iC<nColumns;iC++)
			values[iR][iC] *= d;
}
Matrix Matrix::operator*(Matrix& m) { Matrix r = *this;r *= m;return r; }
Matrix Matrix::operator*(double d) { Matrix r = *this;r *= d;return r; }

Matrix Matrix::operator++() { add(Matrix(nRows, nColumns, MI_VALUE, 1.0));return *this; }

Matrix Matrix::operator++(int) { Matrix C = *this;add(Matrix(nRows, nColumns, MI_VALUE, 1.0));return C; }
Matrix Matrix::operator--() { add(Matrix(nRows, nColumns, MI_VALUE, -1.0));return *this; }
Matrix Matrix::operator--(int) { Matrix r = *this;add(Matrix(nRows, nColumns, MI_VALUE, -1.0));return r; }
Matrix Matrix::operator-()
{
	for (int iR = 0;iR < nRows;iR++)
	{
		for (int iC = 0;iC < nColumns;iC++)
			values[iR][iC] = -values[iR][iC];
	}
	return *this;
}

Matrix Matrix::operator+() { return *this; }

void Matrix::setSubMatrix(int r, int c, Matrix& m)
{
	if ((r + m.nRows)>nRows || (c + m.nColumns)>nColumns)
		throw("Invalid matrix dimension");
	for (int iR = 0;iR<m.nRows;iR++)
		for (int iC = 0;iC<m.nColumns;iC++)
			values[r + iR][c + iC] = m.values[iR][iC];
}

Matrix Matrix::getSubMatrix(int r, int c, int nRows, int nColumns)
{
	if ((r + nRows)>nRows || (c + nColumns)>nColumns)
		throw("Invalid matrix dimension");
	Matrix m(nRows, nColumns);
	for (int iR = 0;iR<m.nRows;iR++)
		for (int iC = 0;iC<m.nColumns;iC++)
			m.values[iR][iC] = values[r + iR][c + iC];
	return m;
}

/*void Matrix::setMatrix()
{

}*/
void Matrix::addColumn(Matrix& m)
{
	Matrix n(max(nRows, m.nRows), nColumns + m.nColumns);
	n.setSubMatrix(0, 0, *this);
	n.setSubMatrix(0, nColumns, m);
	*this = n;
}

void Matrix::addRow(Matrix& m)
{
	Matrix n(nRows + m.nRows, max(nColumns, m.nColumns));
	n.setSubMatrix(0, 0, *this);
	n.setSubMatrix(nRows, 0, m);
	*this = n;
}

Matrix Matrix::getCofactor(int r, int c)
{
	if (nRows <= 1 && nColumns <= 1)
		throw("Invalid matrix dimension");
	Matrix m(nRows - 1, nColumns - 1);
	for (int iR = 0;iR<m.nRows;iR++)
		for (int iC = 0;iC<m.nColumns;iC++)
		{
			int sR = (iR<r) ? iR : iR + 1;
			int sC = (iC<c) ? iC : iC + 1;
			m.values[iR][iC] = values[sR][sC];
		}
	return m;
}

double Matrix::getDeterminant()
{
	if (nRows != nColumns)
		throw("Invalid matrix dimension");
	if (nRows == 1 && nColumns == 1)return values[0][0];
	double value = 0, m = 1;
	for (int iR = 0;iR<nRows;iR++)
	{
		value += m * values[0][iR] * getCofactor(0, iR).getDeterminant();
	m *= -1;
	}
	return value;
}

istream& operator >> (istream &is, Matrix& m)
{
	string s;getline(is, s, ']');s += "]";
	m = Matrix(s);return is;
}
ostream& operator << (ostream &os, Matrix& m)
{
	os << m.getString();
	return os;
}
Matrix Matrix::getInverse()//inverse=(1/determinant)*transpose of cofactor matrix
{
    if(nRows!=nColumns)//inverse can only be done on square matrices
        throw("Invalid Matrix Dimension");
  Matrix n(nRows,nColumns);// copy matrix
  for (int iR = 0;iR<n.nRows;iR++)
		for (int iC = 0;iC<n.nColumns;iC++)
		{
		    n.values[iR][iC]=values[iR][iC];
		}
  double det_value=n.getDeterminant();//determinant value of the matrix
  Matrix m(nRows,nColumns);//cofactor matrix
  int sign=1;
  for (int iR = 0;iR<m.nRows;iR++)
    for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=sign*n.getCofactor(iR,iC).getDeterminant();//getting detreminant values of cofactor matrix
		    sign*=-1;//following sign rule in matrices
		}
     m.getTranspose();//transpose of cofactor matrix
     m*=(1/det_value);
     return m;
}
