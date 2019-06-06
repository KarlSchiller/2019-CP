#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>  // sqrt, round
// #include <complex>
#include <iomanip> // setprecision

using namespace std;
using namespace Eigen;

// Höhenabhängige Dichte
double rho(double h){
  return 1.2*1e3*exp(-1.17*1e-4*h);
}

/* DGL wie angegeben
Cw    Luftwiderstandsbeiwert
A     Querschnittsfläche
m_R   Masse Raumschiff
R     Radius Erde
m_E   Masse Erde
G     Gravitationskonstante

Die beiden Minusse vor den Kräften sind dadurch motiviert, dass der Ursprung
meines Koordinatensystems im Kern der Erde liegt
*/
VectorXd dgl(VectorXd x, VectorXd v){
  double Cw = 0.2, A = 10, m_R = 1e4, R = 6370*1e3, m_E = 5.79*1e24, G = 6.67*1e-11;
  return 1/m_R*(-(Cw*A*rho(x.norm()-R))/2*v.norm()*v - G*m_E*m_R/(x.norm()*x.norm()*x.norm())*x);
}

/* Hilfsfunktion für RK_Schritt
 * Berechnet aus dem Vektor y den Vektor y':
 * Schreibe die letzten d Einträge von y in die ersten d Einträge von y'
 * und die ersten d Einträge von y in die letzten d Einträge von y'
 */
VectorXd next_step(VectorXd y){
  unsigned int d = y.size()/2;
  VectorXd temp(2*d);

  temp.segment(0,d) = y.segment(d,d);
  temp.segment(d,d) = dgl(y.segment(0,d),y.segment(d,d));
  return temp;
}

// Definition für einen einzelnen Runge-Kutta Schritt
VectorXd RK_Schritt(VectorXd y, double h){
  VectorXd k1, k2, k3, k4;
  k1 = h*next_step(y);
  k2 = h*next_step(y+0.5*k1);
  k3 = h*next_step(y+0.5*k2);
  k4 = h*next_step(y+k3);
  return y + 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
}
/* Funktion zum Lösen der DGL mit runge kutta mit Schrittweitenanpassung,
solange, bis die Kapsel die Erde erreicht hat
n = 5 laut Kierfeld Skript für RK 4. Ordnung
*/
VectorXd dgl_loesen(VectorXd y){
  unsigned int d = y.size()/2, i = 0;
  double h = 1, R = 6370*1e3, eps = 1e-4;
  VectorXd y_0 = y, y_1, y_alt;
  double delta_y;

  // Ab hier wird geprüft, ob die Kapsel schon die Erdoberfläche erreicht hat
  while(y.segment(0,d).norm() > R){
    // 2 Schritte mit Schrittweite h
    y = RK_Schritt(y_0, h);
    y_1 = RK_Schritt(y, h);

    // 1 Schritt mir Schrittweite 2h
    y_alt = RK_Schritt(y_0, 2*h);

    // Berechnen des Fehlers
    delta_y = (y_1 - y_alt).norm();

    // Abschätzung mit selbstgewähltem Fehlermaß, an der liegts immerhin nicht
    // Bei diesem Fall wird mit obiges mit dem gleichen y wiederholt, die Schrittweite wurde aber reduziert
    if (abs(delta_y)/y.norm() > eps){
      h = h*pow(eps*y.norm()/abs(delta_y), 1/5);
      // Schnellere Alternativ-Methode zur Implementierung der Formel, läuft leider nicht durch
      // h = h*exp(1.0/5.0*log(eps*y.norm()/delta_y.norm()));
    }
    // In diesem Fall wird wiederholt, und zwar mit einem neuen y_0 und vergrößerter Schrittweite
    else{
      h = h*pow(eps*y.norm()/abs(delta_y), 1/5);
      // Schnellere Alternativ-Methode zur Implementierung der Formel, läuft leider nicht durch
      // h = h*exp(1.0/5.0*log(eps*y.norm()/delta_y.norm()));
      y_0 = y;
    }
    i++;
  }
  cout << "Durchläufe: " << i << endl;
  return y;
}

// Funktion, die zur Minimierung verwendet werden sollte
// Motivation steht in der PDF, v_0 wäre der zu minimierende Parameter
double mini_a(VectorXd y, double v_0){
  unsigned int d = y.size()/2;
  VectorXd y_0 = y, r, v;
  y = dgl_loesen(y_0);
  r = y.segment(0,d);
  v = y.segment(d,d);
  // cout << v.dot(r) << endl;
  return abs(v.dot(r)/r.norm());
}

int main() {
    cout << "Beginn Aufgabe 2!\n" << endl;
    ofstream file;
    unsigned int d = 2, N=1000;
    double R = 6370*1e3, h = 130*1e3;
    VectorXd x(d), x_punkt(d), y(2*d), werte;
    MatrixXd ergebnis(2*d, N);

    // initialisierung des Ortsvektors
    x << 0, R+h;
    y.segment(0,d) = x;

    // Variation der Anfangswerte
    werte = VectorXd::LinSpaced(N, 0, 7.8*1e3);

    file.open("build/aufg2_a.txt", ios::trunc);
    for(int i = 0; i < N; i++){
      x_punkt << werte(i), 0;
      y.segment(d,d) = x_punkt;
      // ergebnis.col(i) = dgl_loesen(y);
      // Aufruf der Funktionen mit Übertrag in ein file
      file << werte(i) << " ";
      file << mini_a(y, 0) << endl;
      // cout << ergebnis.col(i) << endl;
    }
    file.close();
    // file.open("build/aufg2_test.txt", ios::trunc);
    // // Schreiben der Ergebnisse in ein File
    // for(int i = 0; i<ergebnis.rows()/2; i++){
    //   for(int j = 0; j<ergebnis.cols(); j++){
    //     file << ergebnis(i, j) << " ";
    //   }
    //   file << endl;
    // }
    // file.close();

    cout << "\nEnde Aufgabe 2!" << endl;
    return 0;
}
