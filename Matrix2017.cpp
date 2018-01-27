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

void Matrix::copy(const Matrix& m)
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
	strncpy(buffer, s.c_str(), s.length() + 1); //copying the input string to a buffer because strtok will destruct it.

	const char* lineSeparators = ";\r\n"; //separators used to indicate a ln has been terminated
	char* line = strtok(buffer, lineSeparators); //tokenizes the first ln
     	char* Remainlines = strtok(NULL, ""); //tokenizes the remaining lns

	while (line) //line here is my token
	{
		Matrix row; //empty matrix
	        const char* separators = " []"; //row separator is space
						//while [] are used for the first and last rows only (as they have to be removed)
		char* token = strtok(line, separators); //tokenizes the line into numbers (still in a string form)

		while (token) //the token here is the ln I'm extracting numbers from
		{
			const double token_value=atof(token); //converts the tokens into doubles

			Matrix item;
			item = (const double)token_value; //filling the matrix with numbers
			row.addColumn(item); //add each item in its correct column in the (row) matrix

			token = strtok(NULL, separators); //gets the next token (or number in our case)

		}

		if ((row.nColumns>0) && (row.nColumns == nColumns || nRows == 0)) //if there were no rows before in the "this" matrix
			addRow(row); //add this row to "this"

	        line = strtok(Remainlines, lineSeparators); //tokenizing the next ln
        	Remainlines = strtok(NULL, ""); //tokenizing the ln next to it
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
	string s="";
	for (int iR = 0; iR<nRows; iR++)
	{
		for (int iC = 0; iC<nColumns; iC++)
		{
			//cout << values[iR][iC] << " ";

			char buffer[50]="";
			snprintf(buffer, 50, "%g\t", values[iR][iC]);
			s += buffer;
		}
		//cout << endl;
		s += "\n";
	}
	return s;
}

Matrix Matrix::operator=(const Matrix& m)
{
	copy(m);
	return *this;
}

Matrix Matrix::operator=(const double d) { copy(d);return *this; }
Matrix Matrix::operator=(string s) { copy(s);return *this; }

void Matrix::add(const Matrix& m)
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

void Matrix::sub(const Matrix& m)
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
	if (nColumns != m.nRows) //that's how matrices are multiplied
		throw("Invalid matrix dimension for multiplication");

	Matrix r(nRows, m.nColumns); //the dim of the product matrix

	for (int iR = 0; iR<r.nRows; iR++)
	{
		for (int iC = 0; iC < r.nColumns; iC++)
		{
			r.values[iR][iC] = 0; //initializing this particular element of the matrix with zero

			for (int k = 0; k < m.nColumns; k++)
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
Matrix Matrix::operator*(Matrix& m) { Matrix r = *this;	r *= m;	return r; }
Matrix Matrix::operator*(double d) { Matrix r = *this;	r *= d;	return r; }


Matrix Matrix::operator/(Matrix& m) //C = A/B where C, A and B are all matrices
{
	Matrix r = *this;
	Matrix ret;
	ret = r.div(m);
    return ret;
}

Matrix Matrix::operator/(double d) //C = A/B where C, A are matrices and B is a double
{
	Matrix r = *this;
	for (int iR = 0; iR<r.nRows; iR++)
		for (int iC = 0; iC<r.nColumns; iC++)
			r.values[iR][iC] /= d;
	return r;
}

void Matrix::operator/=(Matrix& m) // Divides by m and stores the result in the calling function
{
	*this = div(m);
}

void Matrix::operator/=(double d) // Divides by d (element wise) and stores the result in the calling function
{
	for (int iR = 0; iR<nRows; iR++)
		for (int iC = 0; iC<nColumns; iC++)
			values[iR][iC] /= d;
}

Matrix Matrix::div(Matrix& m)//div C = A/B = A * B.getInverse();
{
	Matrix r = *this;
	if (nColumns != m.nRows)
		throw("First matrix must have the same number of columns as the rows in the second matrix.");
	Matrix t ;
	t= m.getInverse(); // get the inverse of B
	r *= t;
	return r;
}

Matrix Matrix::operator++() { add(Matrix(nRows, nColumns, MI_VALUE, 1.0));return *this; }

Matrix Matrix::operator++(int) { Matrix C = *this;add(Matrix(nRows, nColumns, MI_VALUE, 1.0));return C; }
Matrix Matrix::operator--() { add(Matrix(nRows, nColumns, MI_VALUE, -1.0));return *this; }

//add (-1) to each element of the matrix, then return the matrix.
Matrix Matrix::operator--(int) //int is not used.
{
	Matrix r = *this;
	add(Matrix(nRows, nColumns, MI_VALUE, -1.0));
	return r;
}

//return the same matrix multiplied by a negative sign for each element.
Matrix Matrix::operator-()
{
	for (int iR = 0;iR < nRows;iR++)
	{
		for (int iC = 0;iC < nColumns;iC++)
			values[iR][iC] = -values[iR][iC];
	}
	return *this;
}

//return the same matrix.
Matrix Matrix::operator+() { return *this; }

//copy a submatrix (m) into a matrix, r & c are row & columns where we want to copy.
void Matrix::setSubMatrix(int r, int c, Matrix& m)
{
	if ((r + m.nRows)>nRows || (c + m.nColumns)>nColumns)
		throw("Invalid matrix dimension");
	for (int iR = 0;iR<m.nRows;iR++)
		for (int iC = 0;iC<m.nColumns;iC++)
			values[r + iR][c + iC] = m.values[iR][iC];
}

//extract a submatrix from matrix, r & c are row & column where we want to extract. nRows & nColumns are the rows & columns of the submatrix.
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

//add column to matrix (m).
void Matrix::addColumn(Matrix& m)
{
	Matrix n(max(nRows, m.nRows), nColumns + m.nColumns);
	n.setSubMatrix(0, 0, *this);
	n.setSubMatrix(0, nColumns, m);
	*this = n;
}

//add row to matrix (m).
void Matrix::addRow(Matrix& m)
{
	Matrix n(nRows + m.nRows, max(nColumns, m.nColumns));
	n.setSubMatrix(0, 0, *this);
	n.setSubMatrix(nRows, 0, m);
	*this = n;
}

//return cofactor matrix, r & c are element's row & column which we want to get its cofactor.
Matrix Matrix::getCofactor(int r, int c)
{
	//valid only for (2*2) matrices.
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

//return the determinant of the matrix.
double Matrix::getDeterminant()
{
	double result=0;
	int *p = new int [nRows+1];
		Matrix copy = *this;
		 if(LUPDecompose(copy.values,nRows,0.001,p)){
			 result= LUPDeterminant(copy.values,p,nRows);
			 if(result>0 && result<0.01)result=0;
			 return result;
		 }


}
/* INPUT: A - array of pointers to rows of a square matrix having dimension N
 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
 */
int LUPDecompose(double **A, int N, double Tol, int *P) {

    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0; //failure, matrix is degenerate

        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return 1;  //decomposition done
}

/* INPUT: A,P filled in LUPDecompose; N - dimension.
 * OUTPUT: Function returns the determinant of the initial matrix
 */
double LUPDeterminant(double **A, int *P, int N) {

    double det = A[0][0];

    for (int i = 1; i < N; i++)
        det *= A[i][i];

    if ((P[N] - N) % 2 == 0)
        return det;
    else
        return -det;
}


Matrix Matrix::rdivide(double d, const Matrix& m)
{
    Matrix result(m.nRows,m.nColumns);

    for (int i=0; i<m.nRows; i++)
    {
        for (int j=0; j<m.nColumns; j++)
        {
            result.values[i][j] = ( d / m.values[i][j] );
        }

    }

    return result;
}

istream& operator >> (istream &is, Matrix& m) //inputs the matrix through "cin>>" example: Matrix myMatrix; cin>>myMatrix;
{
	string s;
	getline(is, s, ']'); //] is the delimiter at which the getline knows this is the end of the string
	s += "]"; //because it wasn't saved in the actual string and the constructor uses it
	m = Matrix( s); //uses the constructor which takes an input string
	return is;
}

ostream& operator << (ostream &os, Matrix& m) //prints out the matrix elements using "cout<<"
{
	os << m.getString();
	return os;
}

Matrix Matrix::getInverse() //inverse=(1/determinant)*transpose of cofactor matrix
{
	if (nRows != nColumns) //inverse can only be done on square matrices
		throw("Invalid Matrix Dimension");
	Matrix n=*this;
	double det_value = n.getDeterminant(); //determinant value of the matrix
	Matrix m(nRows, nColumns); //cofactor matrix
	if (det_value==0)
        {
		int x=NAN;
		m=x;
		return m;
	}


	int sign_c =1;
	int sign_r=1;

	for (int iR = 0;iR<m.nRows;iR++)
    {
        sign_c=sign_r;
		for (int iC = 0;iC<m.nColumns;iC++)
		{
			m.values[iR][iC] = sign_c * n.getCofactor(iR, iC).getDeterminant();//getting detreminant values of cofactor matrix
			sign_c *=-1;//following sign rule in matrices
		}
        sign_r*=-1;
    }
	m=m.getTranspose();//transpose of cofactor matrix
	m *= (1 / det_value);
	return m;
}


Matrix Matrix::getTranspose() {
	Matrix x(nColumns, nRows);
	for (int ir = 0; ir<x.nRows;ir++) {
		for (int ic = 0; ic<x.nColumns;ic++) {
			x.values[ir][ic] = values[ic][ir];
		}
	}
	return x;
}
Matrix Matrix::sin(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
            m.values[iR][iC]=std::sin(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::cos(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=std::cos(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::tan(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
            m.values[iR][iC]=std::tan(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::acos(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=std::acos(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::asin(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=std::asin(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::atan(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=std::atan(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::sinh(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=std::sinh(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::cosh(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=std::cosh(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::tanh(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=std::tanh(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::a_sinh(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=asinh(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::a_cosh(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=acosh(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::a_tanh(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
		    m.values[iR][iC]=atanh(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::log(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
            m.values[iR][iC]=std::log(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::log10(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
            m.values[iR][iC]=std::log10(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::exp(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
            m.values[iR][iC]=std::exp(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::sqrt(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
            m.values[iR][iC]=std::sqrt(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::cbrt(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
            m.values[iR][iC]=std::cbrt(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::rpow(Matrix&s,double x)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
            m.values[iR][iC]=std::pow(s.values[iR][iC],x);
		}
    }
    return m;
}
Matrix Matrix::pow(Matrix&s,double x)
{
    if(s.nRows!=s.nColumns)throw("Invalid matrix dimensions");
    if(x<-1)throw("Invalid Power");
    if(x==0)
    {
        Matrix m(s.nRows,s.nColumns,MI_EYE);
        return m;
    }
    else if(x==1)
    {
        Matrix m(s.nRows,s.nColumns);
         m.copy(s);
         return m;
    }
    else if(x==-1)
    {
        Matrix m(s.nRows,s.nColumns);
        m=s.getInverse();
        return m;
    }
    else
    {
      Matrix m(s.nRows,s.nColumns);
      m.copy(s);
      for(int i=1;i<x;i++)
        {
           m*=s;
        }
     return m;
    }
}
Matrix Matrix::floor(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
            m.values[iR][iC]=std::floor(s.values[iR][iC]);
		}
    }
    return m;
}
Matrix Matrix::ceil(Matrix&s)
{
    Matrix m(s.nRows,s.nColumns);
    for (int iR = 0;iR<m.nRows;iR++)
    {
		for (int iC = 0;iC<m.nColumns;iC++)
		{
            m.values[iR][iC]=std::ceil(s.values[iR][iC]);
		}
    }
    return m;
}

//Matrix Matrix::rdivide(const Matrix& m1, const Matrix& m2) //not yet implemented
//{
//    if (m1.nRows != m2.nRows || m1.nColumns != m2.nColumns)
//		throw("Invalid matrix dimension");
//
//    Matrix result;
//
//
//    return result;
//
//}
