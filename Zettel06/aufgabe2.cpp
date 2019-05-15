#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>  // sqrt()
#include <tuple>
#include <iomanip>

using namespace std;
using namespace Eigen;

// Rosenbrock-Funktion
double rosen(double x1, double x2){
  return (1-x1)*(1-x1) + 100*(x2-x1*x1)*(x2-x1*x1);
}

// Funktion aus der b)
double funk_b(double x1, double x2){
  return 1/(1+exp(-10*(x1*x2 -3)*(x1*x2 -3))/(x1*x1 + x2*x2));
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
  double upper = 50, lower=-50, middle = 10, lam;
  //do{
  for(int i = 0; i<11123; i++){
    upper = 50, lower=-50, middle = 0;
    //cout << "iteration: " << i << endl;
    //cout << "upper: " << upper << "lower: " << lower << endl;
    g(0) = first(funptr, x_i(0), x_i(1), 0)*(-1);
    g(1) = first(funptr, x_i(0), x_i(1), 1)*(-1);


    stream << x_i(0) << ";" << g(0) << ";" << x_i(1) << ";" << g(1);

    stream << endl;

    if(minimize(upper, x_i(0), x_i(1), g(0), g(1), rosen) < minimize(middle, x_i(0), x_i(1), g(0), g(1), rosen))
    {
      cout << "ACHTUNG: upper " << setprecision(10) << minimize(upper, x_i(0), x_i(1), g(0), g(1), rosen) << " < middle " << minimize(middle, x_i(0), x_i(1), g(0), g(1), rosen) << endl;
    }
    if(minimize(lower, x_i(0), x_i(1), g(0), g(1), rosen) < minimize(middle, x_i(0), x_i(1), g(0), g(1), rosen))
    {
      cout << "ACHTUNG: lower < middle" << endl;
    }

    bisection(minimize, x_i(0), x_i(1), g(0), g(1), rosen, lower, middle, upper, 1e-8);
    lam = (upper-lower)/2;
    //cout << "Minimale Schrittweite einfach: " << lam << endl;
    x_i = x_i + lam * g;
    //cout << g.norm() << endl;
  //}while(g.norm() >= 1.63);
  }
  return x_i;
}

VectorXd conjugate(double (*funptr)(double, double), VectorXd x0, ofstream &stream){
  // Start
  // Gradienten bestimmen und Variablen initialisieren
  VectorXd g_0(x0.size()), g_i(x0.size());
  VectorXd x_i = x0, p;
  double upper = 100, lower=-100, middle = 0, lam, m;
  int i = 0;

  // Bestimmung des ersten Gradienten
  g_0(0) = first(funptr, x0(0), x0(1), 0)*(-1);
  g_0(1) = first(funptr, x0(0), x0(1), 1)*(-1);

  p = g_0;
  do{
  upper = 50, lower=-50, middle = 0;
  if(minimize(upper, x_i(0), x_i(1), p(0), p(1), funptr) < minimize(middle, x_i(0), x_i(1), p(0), p(1), funptr))
  {
    cout << "ACHTUNG: upper " << setprecision(10) << minimize(upper, x_i(0), x_i(1), p(0), p(1), funptr) << " < middle " << minimize(middle, x_i(0), x_i(1), p(0), p(1), funptr) << endl;
  }
  if(minimize(lower, x_i(0), x_i(1), p(0), p(1), funptr) < minimize(middle, x_i(0), x_i(1), p(0), p(1), funptr))
  {
    cout << "ACHTUNG: lower < middle" << endl;
  }
  // Minimierung eindimensional bezgl. lambda
  bisection(minimize, x_i(0), x_i(1), p(0), p(1), funptr, lower, middle, upper, 1e-8);

  // lambda gewählt als Mitte zwischen upper und lower
  lam = (upper-lower)/2;
  //cout << "minimale Schrittweite konjugiert: " << lam << endl;

  // Update auf nächstes x
  x_i = x_i + lam*p;
  stream << x_i(0) << ";" << g_i(0) << ";" << x_i(1) << ";" << g_i(1);

  stream << endl;

  //cout << x_i << endl;
  g_i(0) = first(funptr, x_i(0), x_i(1), 0)*(-1);
  g_i(1) = first(funptr, x_i(0), x_i(1), 1)*(-1);
  //cout << "g_0: " << g_0 << endl;

  m = (g_i.dot(g_i))/(g_0.dot(g_0));
  g_0 = g_i;
  // Update auf neue Richtung p
  p = g_i + m*p;
  i++;
  //cout << g_i.norm() <<  endl;
}while(g_i.norm() > 0.04); // bei 1.63 beginnts wieder zu steigen
  //cout << i << endl;
  return x_i;
}

int main()
{
    cout << "Beginn des Programms!" << endl;
    VectorXd x0(2);
    x0 << -1, -1;
    ofstream file;
    file.open("build/gradient.txt", ios::trunc);
    file << "# x1 g(x1) x2 g(x2)" << endl;
    file << "x1;g1;x2;g2" << endl;
    cout << "Steepest: " << steepest(rosen, x0, file) << endl;
    cout << endl;
    file.close();

    file.open("build/conjugate.txt", ios::trunc);
    file << "# x1 g(x1) x2 g(x2)" << endl;
    file << "x1;g1;x2;g2" << endl;
    cout << "Konjugiert: " << conjugate(rosen, x0, file) << endl;
    cout << endl;
    file.close();

    // Aufgabenteil b)
    VectorXd x1(2), x2(2), x3(2);
    x1 << 1.5, 2.3;
    x2 << -1.7, -1.9;
    x3 << 0.5, 0.6;

    file.open("build/b1.txt", ios::trunc);
    file << "# x1 g(x1) x2 g(x2)" << endl;
    file << "x1;g1;x2;g2" << endl;
    cout << "Erster Startwert: " << conjugate(funk_b, x1, file) << endl;
    cout << endl;
    file.close();

    file.open("build/b2.txt", ios::trunc);
    file << "# x1 g(x1) x2 g(x2)" << endl;
    file << "x1;g1;x2;g2" << endl;
    cout << "Zweiter Startwert: " << conjugate(funk_b, x2, file) << endl;
    cout << endl;
    file.close();

    file.open("build/b3.txt", ios::trunc);
    file << "# x1 g(x1) x2 g(x2)" << endl;
    file << "x1;g1;x2;g2" << endl;
    cout << "Dritter Startwert: " << conjugate(funk_b, x3, file) << endl;
    cout << endl;
    file.close();
    cout << "Ende des Programms!" << endl;
    return 0;
}
