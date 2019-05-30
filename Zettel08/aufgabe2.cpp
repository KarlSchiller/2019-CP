#include <iostream>
#include <fstream>
#include <Eigen/Dense>
// #include <math.h>  // sqrt, round
// #include <complex>
#include <iomanip> // setprecision

using namespace std;
using namespace Eigen;


// Kraftfeld des Problems
VectorXd feld(VectorXd r, double m, double G, double alpha){
    return -alpha*m*G*r/pow(r.norm(), alpha+2);
}


// Potential des Problems zur Energieüberprüfung
double pot(VectorXd r, double m, double G, double alpha){
    return -m*G/(pow(r.norm(), alpha));
}


/* Hilfsfunktion für runge_kutta
 * Berechnet aus dem Vektor y den Vektor y':
 * Schreibe die letzten d Einträge von y in die ersten d Einträge von y'
 * und die ersten d Einträge von y in die letzten d Einträge von y'
 */
VectorXd next_step(VectorXd y, double m, double G, double alpha){
    unsigned int d = y.size()/2;
    VectorXd temp(2*d);

    temp.segment(0,d) = y.segment(d,d);
    temp.segment(d,d) = 1/m*feld(y.segment(0,d),m,G,alpha);

    return temp;
}


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
* alpha       Anpassung des Potentials
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
        double alpha
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
    energie(0) = 0.5*m*v0.dot(v0) + pot(r0, m, G,alpha);
    r_l = y.segment(0,d);
    v_l = y.segment(d,d);
    drehimpuls.col(0) = m*r_l.cross(v_l);

    // Implementierung des Runge-Kutta-Verfahrens
    // Schritt 0 ist schon gemacht
    for (int i = 1; i < N+1; i++){
        k1 = h*next_step(y, m, G, alpha);
        k2 = h*next_step(y+0.5*k1, m, G, alpha);
        k3 = h*next_step(y+0.5*k2, m, G, alpha);
        k4 = h*next_step(y+k3, m, G, alpha);
        y_next = y + 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
        y = y_next;
        ergebnis.col(i) = y;

        // Gesamtenergie berechnen
        energie(i) = 0.5*m*y.segment(d,d).squaredNorm() + pot(y.segment(0,d), m, G, alpha);

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


/* Kraftfeld Keppler Planet und Mond
 * r0       Vektor des betrachteten Körpers
 * r1       Vektor des anderen Körpers
 * m0       Masse des betrachteten Körpers
 * m1       Masse des anderen Körpers
 * G        Gravitationskonstante
 */
VectorXd feld_2(VectorXd r0, VectorXd r1, double m0, double m1, double G){
    unsigned int d = r0.size();
    VectorXd feld;
    feld = -m0*G*r0/pow(r0.norm(),3);
    feld = feld -m0*m1*G*(r1-r0)/pow(r1.norm()-r0.norm(),3);
    return feld;
}


/* Hilfsfunktion für runge_kuttta Planet und Mond
 * Berechnet aus dem Vektor y den Vektor y':
 * Schreibe die letzten d Einträge von y in die ersten d Einträge von y'
 * und die ersten d Einträge von y in die letzten d Einträge von y'
 * mp       Masse Planet
 * mm       Masse Mond
 * G        Gravitationskonstante
 */
VectorXd next_step_2(VectorXd y, double mp, double mm, double G){
    unsigned int d = y.size()/4;
    VectorXd temp(4*d);

    // Ortsvektoren
    temp.segment(0,d) = y.segment(2*d,d);
    temp.segment(d,d) = y.segment(3*d,d);
    // Geschwindigkeitsvektoren
    temp.segment(2*d,d) = 1/mp*feld_2(y.segment(0,d),y.segment(d,d),mp,mm,G);
    temp.segment(3*d,d) = 1/mm*feld_2(y.segment(d,d),y.segment(0,d),mm,mp,G);

    return temp;
}


/*
Runge Kutte für das Keppler-Problem mit Planet und Mond
* T           Obere Grenze des Zeitintervalls
* h           Schrittweite
* m0          Masse des Planeten
* m1          Masse des Mondes
* G           Gravitationskonstante
* r0          Anfangswert für r des Planeten
* v0          Anfangswert für v des Planeten
* r1          Anfangswert für r des Mondes
* v1          Anfangswert für v des Mondes
* file        File, in das geschrieben wird
* alpha       Anpassung des Potentials
*/
void planet_und_mond(
        double T,
        int N,
        double m0,
        double m1,
        double G,
        VectorXd r0,
        VectorXd v0,
        VectorXd r1,
        VectorXd v1,
        ofstream &file
        )
{
    double h = T/N;
    unsigned int d = r0.size();
    VectorXd tn(N+1), k1, k2, k3, k4, y(4*d), y_next(4*d);
    MatrixXd ergebnis(4*d, N+1);

    // Startwerte
    ergebnis.col(0).segment(0,d) = r0;
    ergebnis.col(0).segment(d,d) = r1;
    ergebnis.col(0).segment(2*d,d) = v0;
    ergebnis.col(0).segment(3*d,d) = v1;
    y.segment(0,d) = r0;
    y.segment(d,d) = r1;
    y.segment(2*d,d) = v0;
    y.segment(3*d,d) = v1;

    // Implementierung des Runge-Kutta-Verfahrens
    // Schritt 0 ist schon gemacht
    for (int i = 1; i < N+1; i++){
        k1 = h*next_step_2(y, m0, m1, G);
        k2 = h*next_step_2(y+0.5*k1, m0, m1, G);
        k3 = h*next_step_2(y+0.5*k2, m0, m1, G);
        k4 = h*next_step_2(y+k3, m0, m1, G);
        y = y + 1.0/6.0*(k1 + 2*k2 + 2*k3 + k4);
        ergebnis.col(i) = y;
    }

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
    double alpha = 1.;         // 1: alpha = 1, 2: alpha = 1.1, 3: alpha = 0.9
    VectorXd r(d), v(d), energie(N+1);
    MatrixXd drehimpuls(d, N+1);
    ofstream file;

    // Aufgabenteil a) mit r und v senkrecht
    r << 1, 0, 0;
    v << 0, 0.5, 0.8;

    file.open("build/aufg2_a_ellipse.txt", ios::trunc);
    runge_kutta(T, N, m, G, r, v, file, energie, drehimpuls, alpha);
    file.close();

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

    // Aufgabenteil a) Problem bei kleinen Anfangsgeschwindigkeiten
    v << 0, 0.005, 0.008;
    file.open("build/aufg2_a_schmetterling.txt", ios::trunc);
    runge_kutta(T, N, m, G, r, v, file, energie, drehimpuls, alpha);
    file.close();
    // nicht vergessen den Vektor wieder zurück zu setzen :D
    v << 0, 0.5, 0.8;

    // Aufgabenteil c) Verändere das Potential
    alpha = 1.1;
    file.open("build/aufg2_c_11.txt", ios::trunc);
    runge_kutta(T, N, m, G, r, v, file, energie, drehimpuls, alpha);
    file.close();

    alpha = 0.9;
    file.open("build/aufg2_c_09.txt", ios::trunc);
    runge_kutta(T, N, m, G, r, v, file, energie, drehimpuls, alpha);
    file.close();

    // Aufgabenteil d) und e) Planet mit Mond
    N = 5000;
    T = 24;
    d = 3;                      // Dimension
    double mp = 2;              // Masse des Planeten
    double mm = 0.01*mp;        // Masse des Mondes
    G = 1;                      // Gravitationskonstante
    VectorXd rp(d), rm(d), vp(d), vm(d);
    rp << 1, 0, 0;
    rm << 1.1, 0, 0;
    vp << 0, sqrt(3/2), 0;
    vm << 0, 4.7, 0;
    file.open("build/aufg2_d.txt", ios::trunc);
    planet_und_mond(T, N, mp, mm, G, rp, vp, rm, vm, file);
    file.close();

    cout << "\nEnde Aufgabe 2!" << endl;
    return 0;
}
