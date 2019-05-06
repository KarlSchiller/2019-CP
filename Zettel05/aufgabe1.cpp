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
    res = res + funptr(x, x_strich, y_strich, a-h/2+i*h);
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

// Im Zähler x_strich oder x?
double funk_c(double x, double x_strich, double y_strich, double z_strich){
  return x_strich/sqrt((x-x_strich)*(x-x_strich) + y_strich*y_strich + z_strich*z_strich);
}


int main() {
  cout << "Beginn des Programms!" << endl;
  // Aufgabenteil a) und b)
  // Initialisieren der Größen
  VectorXd n = VectorXd::LinSpaced(70, 11, 80);
  VectorXd n_b = VectorXd::LinSpaced(11, 0, 10);

  VectorXd x = 0.1*n;
  VectorXd x_b = 0.1*n_b;
  double a = 1.0;

  // Erzeugen des Gitters
  VectorXd y_strich = VectorXd::LinSpaced(100, -10, 10);
  VectorXd x_strich = y_strich;

  /* Problematisch, wenn das Gitter und y nahe 0 werden, dann wird durch
  0 geteilt und das Integral divergiert.
  Noch zu klären: Ist massiv von der Wahl des Gitters abhängig
  Macht einene Unterschied, was ich als Gitter wähle */
  VectorXd y_strich_b = VectorXd::LinSpaced(100, -20, 20);
  VectorXd x_strich_b = y_strich_b;

  double res = 0, res2=0;
  VectorXd pot(x.size()), pot_b(x_b.size()), pot_c_1(x.size()), pot_c_2(x_b.size());
  double N = 11.0;

  for(int i = 0; i<x.size(); i++){
    for(int j=0; j<y_strich.size(); j++){
      res += mittel3(funk_z, -a, a, N, x_strich(j), y_strich(j), x(i));
      res2 += mittel3(funk_c, -a, a, N, x_strich(j), y_strich(j), x(i));
    }
    pot(i) = res;
    pot_c_1(i) = res2;
    res = 0;
    res2 = 0;
  }

  for(int i = 0; i<x_b.size(); i++){
    for(int j=0; j<y_strich_b.size(); j++){
      res += mittel3(funk_z, -a, a, N, x_strich_b(j), y_strich_b(j), x_b(i));
      res2 += mittel3(funk_c, -a, a, N, x_strich_b(j), y_strich_b(j), x_b(i));
    }
    pot_b(i) = res;
    pot_c_2(i) = res2;
    res = 0;
    res2 = 0;
  }

  // Speichern Aufgabenteil a)
  ofstream file;
  file.open("build/aufg_a.txt", ios::trunc);
  file << "# x pot(x)" << endl;
  file << "x;pot" << endl;
  for(int i = 0; i<x.size(); i++){
    file << x(i) << ";";
    file << pot(i) << endl;
  }
  file.close();

  // Speichern Aufgabenteil b)
  file.open("build/aufg_b.txt", ios::trunc);
  file << "# x pot(x)" << endl;
  file << "x;pot" << endl;
  for(int i = 0; i<x_b.size(); i++){
    file << x_b(i) << ";";
    file << pot_b(i) << endl;
  }
  file.close();

  // Speichern Aufgabenteil c)
  file.open("build/aufg_c_1.txt", ios::trunc);
  file << "# x pot(x)" << endl;
  file << "x;pot" << endl;
  for(int i = 0; i<x.size(); i++){
    file << x(i) << ";";
    file << pot_c_1(i) << endl;
  }
  file.close();

  file.open("build/aufg_c_2.txt", ios::trunc);
  file << "# x pot(x)" << endl;
  file << "x;pot" << endl;
  for(int i = 0; i<x_b.size(); i++){
    file << x_b(i) << ";";
    file << pot_c_2(i) << endl;
  }
  file.close();

  //cout << pot << endl;
  cout << "Ende des Programms!" << endl;
  return 0;
}
