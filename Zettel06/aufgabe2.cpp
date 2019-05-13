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

// Eindimensionale, zu minimierende Funktion
double minimize(double lam, double x1, double x2, double g1, double g2,
                double (*funptr)(double, double))
{
  return funptr(x1 + lam*g1, x2 + lam*g2);
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

// Intervalhalbierungsverfahren zur eindimensionalen Minimierung
void bisection(double (*funptr)(double, double, double, double, double, double (*funptr)(double, double)),
                 double x1, double x2, double g1, double g2, double (*funptr1)(double, double),
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
      if(funptr(temp, x1, x2, g1, g2, funptr1) < funptr(middle, x1, x2, g1, g2, funptr1))
      {
        lower = middle;
        middle = temp;
      } else {
        upper = temp;
      }
    } else {
      temp = (middle+lower)/2;
      if(funptr(temp, x1, x2, g1, g2, funptr1) < funptr(middle, x1, x2, g1, g2, funptr1))
      {
        upper = middle;
        middle = temp;
      } else {
        lower = temp;
      }
    }
    N++;
  }
}

VectorXd steepest(double (*funptr)(double, double), VectorXd x0, ofstream &stream){
  // Gradienten bestimmen
  VectorXd g(x0.size());
  VectorXd x_i = x0;
  double upper = 100, lower=0, middle = -3, lam;
  //do{
    g(0) = first(funptr, x_i(0), x_i(1), 0)*(-1);
    g(1) = first(funptr, x_i(0), x_i(1), 1)*(-1);


    stream << x_i(0) << ";" << g(0) << ";" << x_i(1) << ";" << g(1);

    stream << endl;
    bisection(minimize, x_i(0), x_i(1), g(0), g(1), rosen, lower, middle, upper, 1e-10);
    //cout << upper << " " << lower << endl;
    lam = (upper-lower)/2;
    cout << "Minimale Schrittweite: " << lam << endl;
    x_i = x_i + lam * g;
  //}while(g.norm() >= 1);
  return g;
}

// VectorXd conjugate(double (*funptr)(double, double), VectorXd x0, ofstream &stream){
//   // Start
//   // Gradienten bestimmen
//   VectorXd g_0(x0.size()), g_i(x0.size());
//   VectorXd x_i = x0;
//
//   g_0(0) = first(funptr, x0(0), x0(1), 0)*(-1);
//   g_0(1) = first(funptr, x0(0), x0(1), 1)*(-1);
//
//   VectorXd p = g_0;
//   // Hier kommt die eindimensionale Minimierung hin
//   /* l_min = minimierung, aber mit p
//   x_i = x_i + l_min*p
//   g_i(0) = first(funptr, x_i(0), x_i(1), 0)*(-1);
//   g_i(1) = first(funptr, x_i(0), x_i(1), 1)*(-1);
//   double m = (g_i.dot(g_i))/(g_0.dot(g_0));
//   p = g + m*p
//     */
// }

int main()
{
    cout << "Beginn des Programms!" << endl;
    VectorXd x0(2);
    x0 << -1, -1;
    ofstream file;
    file.open("build/gradient.txt", ios::trunc);
    file << "# x1 g(x1) x2 g(x2)" << endl;
    file << "x1; g1; x2; g2" << endl;
    cout << steepest(rosen, x0, file) << endl;
    file.close();
    cout << "Ende des Programms!" << endl;
    return 0;
}
