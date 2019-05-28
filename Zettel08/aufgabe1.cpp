#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>  // sqrt()

using namespace std;
using namespace Eigen;


// // Funktion des harmonischen Oszillators
double feld(double r, double m){
  return -m*r;
}

VectorXd funktion(VectorXd y, double m){
  unsigned int d = y.size()/2;
  VectorXd temp(2*d);

  /* Schreibe die letzten d Einträge von y in die ersten d Einträge von y'
  und die ersten d Einträge von y in die letzten d Einträge von y'*/
  for(int i = 0; i<2*d; i++){
    if (i < d){
      temp(i) = y(i+d);
    }
    else{
      temp(i) = 1/m*feld(y(i-d), m);
    }
  }
  return temp;
}
/*
Funktion zur Implementierung von Runge Kutta
* f     auszuwertende Funktion
* T     t Element [0, T]
* h     Schrittweite
* m     Masse
* r     Anfangswert für r
* v     Anfangswert für v
*/
void runge_kutta(VectorXd (*f)(VectorXd, double), double T, int N, double m, VectorXd r, VectorXd v, ofstream &file){
  //int N = T/h;
  double h = T/N;
  unsigned int d = r.size();
  VectorXd tn(N+1), k1, k2, k3, k4, y(2*d), y_next(2*d);
  MatrixXd ergebnis(2*d, N+1);
  // Die Startwerte in die Matrix schreiben
  ergebnis.col(0).segment(0,d) = r;
  ergebnis.col(0).segment(d,d) = v;

  // Erstellen der Zeiten tn und schreiben in ein file
  for (int i = 0; i<=N; i++){
    tn(i) = i*h;
    file << tn(i) << " ";
  }
  file << endl;
  cout << tn.size() << endl;
  // Initialisieren des y-Vektors
  for (int i = 0; i < 2*d; i++){
    if(i < d){
      y(i) = r(i);
    }
    else{
      y(i) = v(i-d);
    }
  }

  // Implementierung des Runge-Kutta-Verfahrens
  // Schritt 0 ist schon gemacht
  for (int i = 1; i < N+1; i++){
    k1 = h*f(y, m);
    k2 = h*f(y+0.5*k1, m);
    k3 = h*f(y+0.5*k2, m);
    k4 = h*f(y+k3, m);
    y_next = y + 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
    y = y_next;
    ergebnis.col(i) = y;
  }

  // Schreiben der Ergebnisse in ein File
  for(int i = 0; i<ergebnis.rows()/2; i++){
    for(int j = 0; j<ergebnis.cols(); j++){
      file << ergebnis(i, j) << " ";
    }
    file << endl;
  }
}
int main() {
  cout << "Beginn des Programms!\n" << endl;
  // Initialisierung der benötigten Größen
  double T = 20.0;
  int N = 100;
  double m = 2.0;
  unsigned int d = 3;
  VectorXd r(d), v(d);
  MatrixXd ergebnis(2*d, N+1);
  ofstream file;
  // Anfangswerte
  r << 1, 0, 3;
  v << 0, 1, 0;
  
  file.open("build/ergebnis.txt", ios::trunc);
  runge_kutta(funktion, T, N, m, r, v, file);
  file.close();

  cout << "\nEnde des Programms!" << endl;
  return 0;
}
