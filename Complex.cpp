#include "Complex.h"

Complex::Complex()
{
  R = I = 0.0;
}

Complex::Complex(double R, double I)
{this->R = R; this->I = I;}

Complex::Complex(Complex& C)
{
  copy(C);
}

void Complex::copy(Complex& C)
{
  R = C.R; I = C.I;
}

string Complex::getString()
{
  char text[100];
  if(I==0)sprintf(text, "%g", R);
  else if(R==0)sprintf(text, "%g * i", I);
  else if(I>0)sprintf(text, "%g + %g * i", R, I);
  else if(I<0)sprintf(text, "%g - %g * i", R, -I);
  return string(text);
}

double Complex::magnitude()
{
  return sqrt(R*R+I*I);
}
double Complex::angle()
{
  return atan2(I, R);
}
void Complex::negative()
{
  R*=-1; I*=-1;
}
double Complex::real()
{
  return R;
}
double Complex::imaginary()
{
  return I;
}

Complex Complex::addComplex(Complex& A, Complex& B)
{
  Complex C;
  C.R = A.R + B.R;
  C.I = A.I + B.I;
  return C;
}

void Complex::add(Complex& C)
{
  R += C.R;
  I += C.I;
}

Complex Complex::operator=(Complex& C)
{
  copy(C);
  return *this;
}
Complex Complex::operator=(double D)
{
  R = D;
  I = 0;
  return *this;
}
void Complex::operator+=(Complex& C)
{
  add(C);
}
void Complex::operator+=(double D)
{
  R += D;
}
Complex Complex::operator+(Complex& C)
{
  return addComplex(*this, C);
}
Complex Complex::operator+(double D)
{
  return addComplex(*this, Complex(D, 0));
}

Complex Complex::operator-()
{
    R=-R;
    I=-I;
  return *this;
}
Complex::operator const string()
{
  return getString();
}
istream& operator >> (istream& is, Complex& C)
{
  is>>C.R;
  is>>C.I;
  return is;
}
ostream& operator << (ostream& os, Complex& C)
{
  os<<C.getString();
  return os;
}

Complex Complex::operator++()
{
  R++;
  return *this;
}
Complex Complex::operator++(int)
{
  Complex C = *this;
  R+=1;
  return C;
}
double Complex::operator[](string name)
{
  if(name=="magnitude")
    return magnitude();
  if(name=="angle")
    return angle();
  if(name=="real")
    return real();
  if(name=="imaginary")
    return imaginary();
  return 0.0;
}
double Complex::operator()(string name, string info)
{
  if(name=="angle")
  {
    if(info=="degree")
      return angle()*180.0/PI;
      if(info=="radian" || info.length()==0)
      return angle();
  }
  return (*this)[name];
}

//wrong function khalis khalis
Complex Complex::operator*(Complex& A, Complex& B)
{
    Complex Product;
    double RealV = A.real()*B.real() - A.imaginary()*B.imaginary();
    double ImagV = A.real()*B.imaginary() + A.imaginary()*B.real();

    Product(RealV,ImagV);

    return Product;
}
