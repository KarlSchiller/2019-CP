#include <iostream>
// #include <fstream>
#include <Eigen/Dense>
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
  double eps = 1e-8;
  double h;
  if (abs(x) < eps)  // x near zero
  {
    h = sqrt(eps);
  } else {
    h = sqrt(eps)*x;
  }
  return 0.5*(funptr(x+h)-funptr(x-h))/h;
}


// Computes the second derivative of @f at @x
double second(double (*funptr)(double), double x)
{
  double eps = pow(10, -8);
  double h;
  if (abs(x) < eps)  // x near zero
  {
    h = sqrt(eps);
  } else {
    h = sqrt(eps)*x;
  }
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
 */
void bisection(double (*funptr)(double),
                 double &lower,
                 double &middle,
                 double &upper,
                 double tol)
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
    cout << "\t" << N << ": " << lower << "  " << middle << "  " << upper << endl;
  }
}


// TODO: bisection with golden ratio
// double gr = 0.5*(3-sqrt(5));
// just for fun...


/* newton method to minimize function
 * @funptr    pointer to function
 * @x         starting point
 */
void newton(double (*funptr)(double), double &x)
{
  double old;
  int i = 0;
  do {
    cout << "\t" << i << ": " << x << endl;
    old = x;
    x = old - first(funptr, old)/second(funptr, old);
    i++;
  } while (funptr(x) < funptr(old));
  // TODO: Use bisection to find better point and start again
}


int main()
{
  cout << "Aufgabe 1" << endl;

  double accuracy = 1e-9;
  double a = -0.5;
  double b = -0.1;
  double c = 2;
  bisection(quadr, a, b, c, accuracy);
  cout << endl;
  // TODO: print interation steps to file to compare with newton

  double x = 400;
  newton(quadr, x);
  return 0;
  // TODO: print interation steps to file
}
