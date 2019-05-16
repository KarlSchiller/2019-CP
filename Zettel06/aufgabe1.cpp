#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <iomanip>
// #include <math.h>  // sqrt()

using namespace std;
using namespace Eigen;


// function to test minimization algorithms
template<typename T>
T quadr(T x)
{
  return x*x-2;
}


// Computes the first derivative of @f at @x
double first(double (*funptr)(double), double x)
{
  // double eps = 1e-8;
  double h = 1e-2;
  // if (abs(x) < eps)  // x near zero
  // {
    // h = sqrt(eps);
  // } else {
    // h = sqrt(eps)*x;
  // }
  return 0.5*(funptr(x+h)-funptr(x-h))/h;
}


// Computes the second derivative of @f at @x
double second(double (*funptr)(double), double x)
{
  // double eps = pow(10, -8);
  double h = 1e-2;
  // if (abs(x) < eps)  // x near zero
  // {
    // h = sqrt(eps);
  // } else {
    // h = sqrt(eps)*x;
  // }
  return (funptr(x+h)-2*funptr(x)+funptr(x-h))/(h*h);
}


/* bisection method to minimize function
 * @funptr      pointer to function
 * @lower       lower boundary of interval
 * @middle      point in interval with
 *                funptr(middle) < funpter(lower) and
 *                funptr(middle) < funpter(upper)
 * @upper       upper boundary of interval
 * @tol         error tolerance
 * @stream      output stream to print iteration steps
 */
void bisection(double (*funptr)(double),
                 double &lower,
                 double &middle,
                 double &upper,
                 double tol,
                 ostream &stream)
{
  double temp;
  int N = 0;
  while((upper-lower)>tol)
  {
    if(abs(upper-middle)>abs(middle-lower))
    {
      temp = (upper+middle)/2;
      if(funptr(temp) < funptr(middle))
      {
        lower = middle;
        middle = temp;
      } else {
        upper = temp;
      }
    } else {
      temp = (middle+lower)/2;
      if(funptr(temp) < funptr(middle))
      {
        upper = middle;
        middle = temp;
      } else {
        lower = temp;
      }
    }
    N++;
    stream << setprecision(10) << N << ";" << lower << ";" << middle << ";" << upper << endl;
  }
  cout << "\tStop Bisection after " << N << " steps" << endl;
}


/* newton method to minimize function
 * @funptr    pointer to function
 * @x         starting point
 */
void newton(double (*funptr)(double), double &x, double tol, ostream &stream)
{
  double old = x+2*tol;
  int i = 1;
  // do {
  while(abs(x-old) > tol) {
    old = x;
    x = old - first(funptr, old)/second(funptr, old);

    // stop iteration when there is no convergence
    if(first(funptr, old) <= first(funptr, x) && i != 1)
    // if(funptr(old) <= funptr(x))
    {
      cout << "\tStop Newton after " << i-1 << " steps" << endl;
      break;
    }

    stream  << setprecision(10) << i << ";" << x << ";" << first(funptr, x) << endl;
    i++;
  // } while (first(funptr, x) > tol);
  // } while (abs(x - old) > tol);
  }
}


int main()
{
  cout << "Aufgabe 1" << endl;
  ofstream file;

  double accuracy = 1e-9;
  double a = -0.5;
  double b = -0.1;
  double c = 2;
  file.open("build/aufg1-bisection.txt");
  file << "i;lower;middle;upper" << endl;
  file << "0;" << a << ";" << b << ";" << c << endl;
  bisection(quadr, a, b, c, accuracy, file);
  file.close();

  double x = 400;
  file.open("build/aufg1-newton.txt");
  file << "i;x;deriv" << endl;
  file << "0;" << x << ";" << first(quadr, x) << endl;
  newton(quadr, x, accuracy, file);
  file.close();
  return 0;
}
