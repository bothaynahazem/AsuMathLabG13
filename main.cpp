#include "Matrix2017.h"

int main()
{
	Matrix mymatrix1("[25 25; 25 25]");
	Matrix mymatrix2("[5 5; 1 25]");
	Matrix mymatrix3;
	mymatrix3=mymatrix1/mymatrix2;

	cout<<mymatrix3;

	int x;
//	cin >> x;
	return 0;
}
