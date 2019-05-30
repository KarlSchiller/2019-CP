#include <iostream>
#include <fstream>
#include <Eigen/Dense>
// #include <math.h>  // sqrt, round
// #include <tuple>
// #include <complex>
#include <iomanip> // setprecision

using namespace std;
using namespace Eigen;


// Kraftfeld des Problems
VectorXd feld(VectorXd r, double m, double G, int Fall){
  double r_norm = r.norm();
  switch(Fall){
    case 1: return -m*G*r/r_norm; break;
    case 2: return -m*G*r*1.1/(pow(r.norm(), 3.1)); break;
    case 3: return -m*G*r*0.9/(pow(r.norm(), 2.9)); break;
    default: return -m*G*r/r_norm; // default mäßig den einfachen Fall zurück geben
  }
    //return -m*G*r/r_norm;

    // Fall für alpha = 1.1
    //return -m*G*r*1.1/(pow(r.norm(), 3.1));

    // Fall für alpha = 0.9
    //return -m*G*r*0.9/(pow(r.norm(), 2.9));
}


// Potential des Problems zur Energieüberprüfung
double pot(VectorXd r, double m, double G, int Fall){
  switch(Fall){
    case 1: return -m*G/r.norm(); break;
    case 2: return -m*G/(pow(r.norm(), 1.1)); break;
    case 3: return -m*G/(pow(r.norm(), 0.9)); break;
    default: return -m*G/r.norm(); // default mäßig den einfachen Fall zurück geben
  }
    // return 0.5*m*r.dot(r);
    //return -m*G/r.norm();

    // Fall für alpha =
    //return -m*G/(pow(r.norm(), 0.9));
    // return -m*G/(pow(r.norm(), 1.1));
}


/* Hilfsfunktion für runge_kutta
 * Berechnet aus dem Vektor y den Vektor y':
 * Schreibe die letzten d Einträge von y in die ersten d Einträge von y'
 * und die ersten d Einträge von y in die letzten d Einträge von y'
 */
VectorXd next_step(VectorXd y, double m, double G, int Fall){
    unsigned int d = y.size()/2;
    VectorXd temp(2*d);

    temp.segment(0,d) = y.segment(d,d);
    temp.segment(d,d) = 1/m*feld(y.segment(0,d),m,G, Fall);

    return temp;
}


// Aus Aufgabe 1 kopiert
/*
Funktion zur Implementierung von Runge Kutta
* T           Obere Grenze des Zeitintervalls
* h           Schrittweite
* m           Masse des Planeten
* G           Gravitationskonstante
* r0          Anfangswert für r
* v0          Anfangswert für v
* file        File, in das geschrieben wird
* energie     Gesamtenergie des Oszillators
* drehimpuls  Matrix mit den Drehimpulsvektoren
* Fall        1: alpha = 1, 2: alpha = 1.1, 3: alpha = 0.9
*/
void runge_kutta(
        double T,
        int N,
        double m,
        double G,
        VectorXd r0,
        VectorXd v0,
        ofstream &file,
        VectorXd &energie,
        MatrixXd &drehimpuls,
        int Fall
        )
{
    double h = T/N;
    unsigned int d = r0.size();
    VectorXd tn(N+1), k1, k2, k3, k4, y(2*d), y_next(2*d);
    MatrixXd ergebnis(2*d, N+1);
    Vector3d r_l, v_l;  // Hilfsgrößen zur Berechnung des Drehimpulses

    // Startwerte
    ergebnis.col(0).segment(0,d) = r0;
    ergebnis.col(0).segment(d,d) = v0;
    y.segment(0,d) = r0;
    y.segment(d,d) = v0;
    energie(0) = 0.5*m*v0.dot(v0) + pot(r0, m, G, Fall);
    r_l = y.segment(0,d);
    v_l = y.segment(d,d);
    drehimpuls.col(0) = m*r_l.cross(v_l);

    // Implementierung des Runge-Kutta-Verfahrens
    // Schritt 0 ist schon gemacht
    for (int i = 1; i < N+1; i++){
        k1 = h*next_step(y, m, G, Fall);
        k2 = h*next_step(y+0.5*k1, m, G, Fall);
        k3 = h*next_step(y+0.5*k2, m, G, Fall);
        k4 = h*next_step(y+k3, m, G, Fall);
        y_next = y + 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
        y = y_next;
        ergebnis.col(i) = y;

        // Gesamtenergie berechnen
        energie(i) = 0.5*m*y.segment(d,d).squaredNorm() + pot(y.segment(0,d), m, G, Fall);

        // Drehimpuls berechnen
        r_l = y.segment(0,d);
        v_l = y.segment(d,d);
        drehimpuls.col(i) = m*r_l.cross(v_l);
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
    double T = 200.0;      // obere Grenze des Zeitintervalls
    int N = 3000;         // Anzahl Schritte
    double m = 2.0;       // Masse
    double G = 1;         // Gravitationskonstante
    unsigned int d = 3;   // Dimension
    int Fall = 1;         // 1: alpha = 1, 2: alpha = 1.1, 3: alpha = 0.9
    VectorXd r(d), v(d), energie(N+1);
    MatrixXd drehimpuls(d, N+1);
    ofstream file;

    // Aufgabenteil a) mit r und v senkrecht
    r << 1, 0, 0;
    v << 0, 0.5, 0.8;

    file.open("build/aufg2_a_ellipse.txt", ios::trunc);
    runge_kutta(T, N, m, G, r, v, file, energie, drehimpuls, Fall);
    file.close();
    //cout << energie << endl;

    // Aufgabenteil b): Prüfe Energieerhaltung
    file.open("build/aufg2_b_energie.txt", ios::trunc);
    file << "zeit energie" << endl;
    for(int i=0; i<=N; i++)
    {
        file << setprecision(10) << i*T/N << " " << energie(i) << endl;
    }
    file.close();

    // Aufgabenteil b): Prüfe Drehimpulserhaltung
    file.open("build/aufg2_b_drehimpuls.txt", ios::trunc);
    file << "zeit Lx Ly Lz" << endl;
    for(int i=0; i<=N; i++)
    {
        file << setprecision(10) << i*T/N;
        for(int dim=0; dim<d; dim++)
        {
            file << " " << drehimpuls(dim, i);
        }
        file << endl;
    }
    file.close();

    // Problem bei kleinen Anfangsgeschwindigkeiten
    v << 0, 0.005, 0.008;
    file.open("build/aufg2_a_schmetterling.txt", ios::trunc);
    runge_kutta(T, N, m, G, r, v, file, energie, drehimpuls, Fall);
    file.close();

    // Aufgabenteil c)
    // alpha = 1.1
    v << 0, 0.5, 0.8; // nicht vergessen den Vektor wieder zurück zu setzen :D
    Fall = 2;
    file.open("build/aufg2_c_11.txt", ios::trunc);
    runge_kutta(T, N, m, G, r, v, file, energie, drehimpuls, Fall);
    file.close();

    // alpha = 0.9
    Fall = 3;
    file.open("build/aufg2_c_09.txt", ios::trunc);
    runge_kutta(T, N, m, G, r, v, file, energie, drehimpuls, Fall);
    file.close();

    cout << "\nEnde Aufgabe 2!" << endl;
    return 0;
}
