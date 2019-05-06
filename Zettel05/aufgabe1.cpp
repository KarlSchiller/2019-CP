#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>  // sqrt()

using namespace std;
using namespace Eigen;

// Mittelpunktsregel
double mittel(double (*funptr)(double, double, double, double), double a, double b, double N,
double x_strich, double y_strich, double x){
  double h = (b-a)/N;
  double res=0;
  for(int i = 1; i<N+1; i++){
    res = res + funptr(x, x_strich ,y_strich, a-h/2+i*h);
  }
  res *= h;
  return res;
}

// Mittelpunktsregel Nummer 2
double mittel2(double (*funptr)(double, double, double, double), double a, double b, double N,
double x_strich, double y_strich, double x){
  double h = (b-a)/N;
  double res=0;
  for(int i = 1; i<N+1; i++){
    res = res + mittel(funptr, a, b, N, x_strich, y_strich, x);
  }
  res *= h;
  return res;
}

// Mittelpunktsregel Nummer 3
double mittel3(double (*funptr)(double, double, double, double), double a, double b, double N,
double x_strich, double y_strich, double x){
  double h = (b-a)/N;
  double res=0;
  for(int i = 1; i<N+1; i++){
    res = res + mittel2(funptr, a, b, N, x_strich, y_strich, x);
  }
  res *= h;
  return res;
}

double funk_z(double x, double x_strich, double y_strich, double z_strich){
  return 1/sqrt((x-x_strich)*(x-x_strich) + y_strich*y_strich + z_strich*z_strich);
}


int main() {
  cout << "Beginn des Programms!" << endl;
  // Aufgabenteil a)
  // Initialisieren der Größen
  VectorXd n = VectorXd::LinSpaced(70, 11, 80);
  VectorXd x = 0.1*n;
  double a = 1.0;

  // Erzeugen
  VectorXd y_strich = VectorXd::LinSpaced(x.size(), -10, 10);
  VectorXd x_strich = y_strich;
  double res = 0;
  VectorXd pot(x.size());
  double N = 11.0;

  for(int i = 0; i<x.size(); i++){
    for(int j=0; j<y_strich.size(); j++){
      res += mittel3(funk_z, -a, a, N, x_strich(j), y_strich(j), x(i));
    }
    pot(i) = res;
    res = 0;
  }

  ofstream file;
  file.open("build/aufg_a.txt", ios::trunc);
  file << "# x pot(x)" << endl;
  file << "x;pot" << endl;
  for(int i = 0; i<x.size(); i++){
    file << x(i) << ";";
    file << pot(i) << endl;
  }
  file.close();

  //cout << pot << endl;
  cout << "Ende des Programms!" << endl;
  return 0;
}
