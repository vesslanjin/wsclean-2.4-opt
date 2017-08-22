#include <iostream>

#include "../../polynomialfitter.h"
#include "../../polynomialchannelfitter.h"

using namespace std;

double y(double x)
{
  return (-1.0/29700.0)*x*x + 82.0/33.0;
}

double inty(double x)
{
  return (-1.0/29700.0/3.0) * x*x*x + (82.0/33.0) * x;
}

int main(int argc, char* argv[])
{
  // function:
  // x=120, y=2
  // x=210, y=1
  // a x*x + b = y
  // a 14400 + b = 2, b = 2 - 14400a
  // a 44100 + b = 1
  // a 44100 + 2 - 14400a = 1
  // 29700a = -1, a = -1/29700
  // b = 2+14400/29700 = 82/33
  // y = -1/29700 x*x + 82/33

  // int_l^m y = [ -1/9900 x*x*x +  82/33 x ]_l^m
  
  // Divide into 120-150, 150-180, 180-210

  double y1 = y(135);
  double y2 = y(165);
  double y3 = y(195);
  std::cout.precision(14);
  std::cout << " 120-150: " << y1 << '\n';
  std::cout << " 150-180: " << y2 << '\n';
  std::cout << " 180-210: " << y3 << '\n';
  
  double iy1 = (inty(150) - inty(120)) / 30.0;
  double iy2 = (inty(180) - inty(150)) / 30.0;
  double iy3 = (inty(210) - inty(180)) / 30.0;
  std::cout << "int 120-150: " << iy1 << '\n';
  std::cout << "int 150-180: " << iy2 << '\n';
  std::cout << "int 180-210: " << iy3 << '\n';

  ao::uvector<double> aTerms{(82.0/33.0), 0.0, (-1.0/29700.0) };
  std::cout << "Real: y = " << aTerms[2] << "*x^2 + " << aTerms[0] << '\n';

  PolynomialFitter pFitter;
  pFitter.AddDataPoint(135, iy1, 1.0);
  pFitter.AddDataPoint(165, iy2, 1.0);
  pFitter.AddDataPoint(195, iy3, 1.0);
  ao::uvector<double> pTerms;
  pFitter.Fit(pTerms, 3);
  std::cout << "Polynomial regression: y = " << pTerms[2] << "*x^2 + " << pTerms[1] << " + " << pTerms[0] << '\n';
  std::cout << "Errors: " << pTerms[2]-aTerms[2] << ", " << pTerms[1]-aTerms[1] << ", " << pTerms[0]-aTerms[0] << '\n';
  std::cout << pTerms[2]*29700.0 << "*x*x/29700 + " << pTerms[1]*1000.0 << " x/1000 + " << pTerms[0]*33.0 << "/33\n";

  
  PolynomialChannelFitter cFitter;
  cFitter.AddChannel(120,150);
  cFitter.AddChannel(150,180);
  cFitter.AddChannel(180,210);
  cFitter.AddDataPoint(0, iy1);
  cFitter.AddDataPoint(1, iy2);
  cFitter.AddDataPoint(2, iy3);
  ao::uvector<double> cTerms;
  cFitter.Fit(cTerms, 3);
  std::cout << "Channel-aware fit: y = " << cTerms[2] << "*x^2 + " << cTerms[1] << " + " << cTerms[0] << '\n';
  std::cout << "Errors: " << cTerms[2]-aTerms[2] << ", " << cTerms[1]-aTerms[1] << ", " << cTerms[0]-aTerms[0] << '\n';
  
  std::cout << cTerms[2]*29700.0 << "*x*x/29700 + " << cTerms[1]*1000.0 << " x/1000 + " << cTerms[0]*33.0 << "/33\n";
}
