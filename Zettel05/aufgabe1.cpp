#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>  // sqrt()

using namespace std;
using namespace Eigen;

// Simpsonregel
double simpson(double (*funptr)(double, double, double, double), double a,
double b, int N, double x, double x_strich, double y_strich){
  double h = (b-a)/N;
  double n = N/2;
  double res=0, x2k=0, x2k1 = 0;
  if(fmod(N, 2) ==0){
    for(int i = 1; i<n; i++){
      x2k = a + i*2*h;
      res += 2*funptr(x, x_strich, y_strich,x2k);
    }
    for(int i = 1; i<n+1; i++){
        x2k1 = a + (2*i-1)*h;
        res += 4*funptr(x, x_strich, y_strich, x2k1);
    }
    res += funptr(x, x_strich, y_strich, a) + funptr(x, x_strich, y_strich, b);
    res *= h/3;
  }
  return res;
}

double funk_z(double x, double x_strich, double y_strich, double z_strich){
  return 1/sqrt((x-x_strich)*(x-x_strich) + y_strich*y_strich + z_strich*z_strich);
}


int main() {
  cout << "Beginn des Programms!" << endl;
  // Initialisieren der Größen
  VectorXd n = VectorXd::LinSpaced(70, 11, 80);
  VectorXd x = VectorXd::LinSpaced(70, -10, 10);
  VectorXd y_strich = VectorXd::LinSpaced(70, -10, 10);
  VectorXd x_strich = y_strich;
  double res = 0;
  double zwischen, N = 2.0;

  res = simpson(funk_z, -x(0)/n(0), x(0)/n(0), N, x(0), x_strich(0), y_strich(0));
  zwischen = res;
  // Aufrufen der Funktion
  for(int i = 0; i<70; i++){
    return 0;
  }

  cout << "Ende des Programms!" << endl;
  return 0;
}
