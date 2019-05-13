#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>  // sqrt()
#include <tuple>

using namespace std;
using namespace Eigen;

// Rosenbrock-Funktion
double rosen(double x1, double x2){
  return (1-x1)*(1-x1) + 100*(x2-x1*x1)*(x2-x1*x1);
}

/* Computes the first derivative of @f at @x, k ist nur zur Festlegung,
welche Variable grade abgeleitet wird*/
double first(double (*funptr)(double, double), double x1, double x2, int k)
{
  double eps = 1e-8;
  double h;
  if ((abs(x1) < eps) || (abs(x2) < eps))  // x near zero
  {
    h = sqrt(eps);
  } else if (k == 0){
    h = sqrt(eps)*x1;
  }
  else {
    h = sqrt(eps)*x2;
  }
  // Für beide Variablen auswerten
  if (k == 0){
    return 0.5*(funptr(x1+h, x2)-funptr(x1-h, x2))/h;
  }
  else{
    return 0.5*(funptr(x1, x2+h)-funptr(x1, x2-h))/h;
  }
}

VectorXd steepest(double (*funptr)(double, double), VectorXd x0){
  // Gradienten bestimmen
  VectorXd g(x0.size());
  VectorXd x_i = x0;
  g(0) = first(funptr, x_i(0), x_i(1), 0);
  g(1) = first(funptr, x_i(0), x_i(1), 1);
  return g;
}

int main()
{
    cout << "Beginn des Programms!" << endl;
    VectorXd x0(2);
    x0 << -1, -1;
    cout << steepest(rosen, x0) << endl;
    cout << "Ende des Programms!" << endl;
    return 0;
}
