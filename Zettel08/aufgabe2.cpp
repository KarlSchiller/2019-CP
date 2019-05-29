#include <iostream>
#include <fstream>
#include <Eigen/Dense>
// #include <math.h>  // sqrt, round
// #include <tuple>
// #include <complex>
// #include <iomanip> // setprecision

using namespace std;
using namespace Eigen;


// Kraftfeld des Problems
VectorXd feld(VectorXd r, double m){
  return -m*r;
}


// Potential des Problems zur Energieüberprüfung
double pot(VectorXd r, double m){
  return 0.5*m*r.dot(r);
}


/* Hilfsfunktion für runge_kutta
 * Berechnet aus dem Vektor y den Vektor y':
 * Schreibe die letzten d Einträge von y in die ersten d Einträge von y'
 * und die ersten d Einträge von y in die letzten d Einträge von y'
 */
VectorXd next_step(VectorXd y, double m){
  unsigned int d = y.size()/2;
  VectorXd temp(2*d);

  temp.segment(0,d) = y.segment(d,d);
  temp.segment(d,d) = 1/m*feld(y.segment(0,d),m);

  return temp;
}


// Aus Aufgabe 1 kopiert
/*
Funktion zur Implementierung von Runge Kutta
* T         Obere Grenze des Zeitintervalls
* h         Schrittweite
* m         Masse des Oszillators
* r0        Anfangswert für r
* v0        Anfangswert für v
* file      File, in das geschrieben wird
* energie   Gesamtenergie des Oszillators
*/
void runge_kutta(
        double T,
        int N,
        double m,
        VectorXd r0,
        VectorXd v0,
        ofstream &file,
        VectorXd &energie
        )
{
  double h = T/N;
  unsigned int d = r0.size();
  VectorXd tn(N+1), k1, k2, k3, k4, y(2*d), y_next(2*d);
  MatrixXd ergebnis(2*d, N+1);
  // Die Startwerte in die Matrix schreiben
  ergebnis.col(0).segment(0,d) = r0;
  ergebnis.col(0).segment(d,d) = v0;

  // Initialisieren des y-Vektors
  for (int i = 0; i < 2*d; i++){
    if(i < d){
      y(i) = r0(i);
    }
    else{
      y(i) = v0(i-d);
    }
  }
  energie(0) = 0.5*m*v0.dot(v0) + pot(r0, m);
  // Implementierung des Runge-Kutta-Verfahrens
  // Schritt 0 ist schon gemacht
  for (int i = 1; i < N+1; i++){
    k1 = h*next_step(y, m);
    k2 = h*next_step(y+0.5*k1, m);
    k3 = h*next_step(y+0.5*k2, m);
    k4 = h*next_step(y+k3, m);
    y_next = y + 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
    y = y_next;
    ergebnis.col(i) = y;
    // Gesamtenergie berechnen
    energie(i) = 0.5*m*y.segment(d,d).dot(y.segment(d,d)) + pot(y.segment(0,d), m);
  }
  // cout << "Energie" << endl << energie << endl;

  // Speichern der Zeiten
  for (int i = 0; i<=N; i++){
    tn(i) = i*h;
    file << tn(i) << " ";
  }
  file << endl;

  // Speichern der Ergebnisse
  for(int i = 0; i<ergebnis.rows()/2; i++){
    for(int j = 0; j<ergebnis.cols(); j++){
      file << ergebnis(i, j) << " ";
    }
    file << endl;
  }
}

int main() {
  cout << "Beginn Aufgabe 2!\n" << endl;

  // Initialisierung der benötigten Größen
  double T = 20.0;      // obere Grenze des Zeitintervalls
  int N = 20;          // Anzahl Schritte
  double m = 2.0;       // Masse
  unsigned int d = 3;   // Dimension
  VectorXd r(d), v(d), energie(N+1);
  ofstream file;
  // Aufgabenteil a) mit r und v senkrecht
  r << 1, 0, 3;
  v << 0, 1, 0;

  file.open("build/aufg2_a_unharm.txt", ios::trunc);
  runge_kutta(T, N, m, r, v, file, energie);
  file.close();

  cout << "\nEnde Aufgabe 2!" << endl;
  return 0;
}
