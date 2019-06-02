#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <math.h>  // sqrt()
#include <time.h> // für Zeitmessungen

using namespace std;
using namespace Eigen;

// Potential zur DGL in Arbeit
// double pot(VectorXd x, double alpha){
//   return 0.5*x.squaredNorm() + alpha*x;
// }

// DGL wie angegeben
VectorXd dgl(VectorXd x, VectorXd x_punkt, double alpha){
  return -x-alpha*x_punkt;
}

/* Hilfsfunktion für adams-bashfort
 * Berechnet aus dem Vektor y den Vektor y':
 * Schreibe die letzten d Einträge von y in die ersten d Einträge von y'
 * und die ersten d Einträge von y in die letzten d Einträge von y'
 */
VectorXd next_step(VectorXd y, double alpha){
  unsigned int d = y.size()/2;
  VectorXd temp(2*d);

  temp.segment(0,d) = y.segment(d,d);
  temp.segment(d,d) = dgl(y.segment(0,d), y.segment(d,d),alpha);
  return temp;
}


MatrixXd runge_kutta(VectorXd (*f)(VectorXd, double), double T, int N, double alpha,
VectorXd r, VectorXd v){
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
    //file << tn(i) << " ";
  }
  //file << endl;

  // Initialisieren des y-Vektors
  y.segment(0,d) = r;
  y.segment(d,d) = v;

  // Berechnung der Energie aus kinetischer und potentieller Energie
  //energie(0) = 0.5*m*v.dot(v) + pot(r, m);

  // Implementierung des Runge-Kutta-Verfahrens
  // Schritt 0 ist schon gemacht
  for (int i = 1; i < N+1; i++){
    k1 = h*f(y, alpha);
    k2 = h*f(y+0.5*k1, alpha);
    k3 = h*f(y+0.5*k2, alpha);
    k4 = h*f(y+k3, alpha);
    y_next = y + 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
    y = y_next;
    ergebnis.col(i) = y;
    // Gesamtenergie berechnen
    //energie(i) = 0.5*m*y.segment(d,d).dot(y.segment(d,d)) + pot(y.segment(0,d), m);
  }

  // Schreiben der Ergebnisse in ein File
  for(int i = 0; i<ergebnis.rows()/2; i++){
    for(int j = 0; j<ergebnis.cols(); j++){
      //file << ergebnis(i, j) << " ";
    }
    //file << endl;
  }
  return ergebnis;
}
/*
Funktion zur Implementierung von Adams-Bashfort
* f         auszuwertende Funktion
* T         obere Grenze des Zeitintervalls
* h         Schrittweite
* m         Masse
* r         Anfangswert für r
* v         Anfangswert für v
* file      File, in das geschrieben wird
* energie   Gesamtenergie des Oszillators
*/
void adams_bashfort(VectorXd (*f)(VectorXd, double), double T, int N, double alpha,
VectorXd x, VectorXd x_punkt, ofstream &file, VectorXd &energie){
  double h = T/N;
  unsigned int d = x.size();
  VectorXd y(2*d), y_next(2*d), tn(N+1);
  MatrixXd ergebnis(2*d, N+1), runge(2*d, 4);

  y.segment(0,d) = x;
  y.segment(d,d) = x_punkt;
  runge = runge_kutta(next_step, T*3.0/N, 3, alpha, x, x_punkt);


  // Schritte, in denen ich noch die Ergebnisse aus Runge-Kutta verwenden muss
  ergebnis.col(0) = runge.col(3) + h/24.0*(55*f(runge.col(3), alpha) - 59*f(runge.col(2), alpha)
                    +37*f(runge.col(1), alpha) -9*f(runge.col(0), alpha));

  ergebnis.col(1) = ergebnis.col(0) + h/24.0*(55*f(ergebnis.col(0), alpha) - 59*f(runge.col(3), alpha)
                    +37*f(runge.col(2), alpha) -9*f(runge.col(1), alpha));

  ergebnis.col(2) = ergebnis.col(1) + h/24.0*(55*f(ergebnis.col(1), alpha) - 59*f(ergebnis.col(0), alpha)
                    +37*f(runge.col(3), alpha) -9*f(runge.col(2), alpha));

  ergebnis.col(3) = ergebnis.col(2) + h/24.0*(55*f(ergebnis.col(2), alpha) - 59*f(ergebnis.col(1), alpha)
                    +37*f(ergebnis.col(0), alpha) -9*f(runge.col(3), alpha));

  //ToDO: Potential bestimmen
  energie(0) = 0.5*ergebnis.col(0).segment(d,d).squaredNorm();
  energie(1) = 0.5*ergebnis.col(1).segment(d,d).squaredNorm();
  energie(2) = 0.5*ergebnis.col(2).segment(d,d).squaredNorm();
  energie(3) = 0.5*ergebnis.col(3).segment(d,d).squaredNorm();
  // Ab hier Automatisierung mit ergebnis
  for(int i = 4; i < N+1; i++){
    y_next = ergebnis.col(i-1) + h/24.0*(55*f(ergebnis.col(i-1), alpha) - 59*f(ergebnis.col(i-2), alpha)
    + 37*f(ergebnis.col(i-3), alpha) - 9*f(ergebnis.col(i-4), alpha));
    ergebnis.col(i) = y_next;
    energie(i) = 0.5*y.segment(d,d).dot(y.segment(d,d));
  }

  // Erstellen der Zeiten tn und schreiben in ein file
  for (int i = 0; i<=N; i++){
    tn(i) = i*h;
    file << tn(i) << " ";
  }
  file << endl;

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

  unsigned int d = 3;
  double alpha = 0.1, T = 20.0;
  int N = 300;
  VectorXd x(d), x_punkt(d), y(2*d), energie(N+1);
  ofstream file;
  MatrixXd test;

  x << 1, 2, 3;
  x_punkt << 0, 1, 0;
  y.segment(0,d) = x;
  y.segment(d,d) = x_punkt;

  // Aufgabenteil a)
  file.open("build/aufg1_a1.txt", ios::trunc);
  adams_bashfort(next_step, T, N, alpha, x, x_punkt, file, energie);
  file.close();

  // Kurz Aufgabenteil b) eingeschoben :)
  file.open("build/aufg1_b.txt", ios::trunc);
  file << "zeit energie" << endl;
  for(int i=0; i<=N; i++)
  {
      file << setprecision(10) << i*T/N << " " << energie(i) << endl;
  }
  file.close();

  alpha = -0.1;
  file.open("build/aufg1_a2.txt", ios::trunc);
  adams_bashfort(next_step, T, N, alpha, x, x_punkt, file, energie);
  file.close();

  alpha = 0.0;
  file.open("build/aufg1_a3.txt", ios::trunc);
  adams_bashfort(next_step, T, N, alpha, x, x_punkt, file, energie);
  file.close();

  // Aufgabenteil c)
  double time_a = 0.0, time_r = 0.0, tstart;
  T = 2000.0, alpha = 0.1;
  N = 30000;
  VectorXd energie2(N+1);

  tstart = clock();
  adams_bashfort(next_step, T, N, alpha, x, x_punkt, file, energie2);
  time_a = clock() - tstart;
  time_a = time_a/CLOCKS_PER_SEC;
  cout << "Zeit Adams: " << time_a << ".sec" << endl;

  tstart = clock();
  runge_kutta(next_step, T, N, alpha, x, x_punkt);
  time_r = clock() - tstart;
  time_r = time_r/CLOCKS_PER_SEC;
  cout << "Zeit Runge: " << time_r << ".sec" << endl;

  cout << "\nEnde des Programms!" << endl;
  return 0;
}
