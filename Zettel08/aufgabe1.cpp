#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>  // sqrt()

using namespace std;
using namespace Eigen;


// Fuktion des harmonischen Oszillators
VectorXd func(VectorXd r, double m){
  return -m*r;
}
/*
Funktion zur Implementierung von Runge Kutta
* f     auszuwertende Funktion
* T     t Element [0, T]
* h     Schrittweite
*/
VectorXd runge_kutta(VectorXd (*f)(VectorXd, double), double T, double h){
  int N = T/h;
  VectorXd tn(N+1), k1, k2, k3, k4;
  for (int i = 0; i<=N; i++){
    tn(i) = i*h;
  }

  //k1 = h*f()
  return tn;
}
int main() {
  cout << "Beginn des Programms!" << endl;
  double T = 20.0;
  double h = 2.0;
  cout << runge_kutta(func, T, h) << endl;
  cout << "\nEnde des Programms!" << endl;
  return 0;
}
