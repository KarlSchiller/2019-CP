#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <math.h>  // sqrt()
#include <random>
#include <string>


using namespace std;
using namespace Eigen;


/*  Berechnet die Kraft des Lennard-Jones-Potentials
 *  Keine Kraft bei Abstaenden groesser des Cutoff
 *  INPUT       r       2D-Ortsvektor zwischen zwei Teilchen
 *              rc      Cutoff
 *  OUTPUT      ljkraft Kraftvektor nach dem LJ-Potential
 */
VectorXd lj_kraft(VectorXd r, double rc)
{
    VectorXd ljk = VectorXd::Zero(2,1);
    if (r.norm()>rc)
    {
        ljk << 0, 0;
    }
    else
    {
        ljk = r*24*(2*pow(r.norm(),-14)-pow(r.norm(),-8));
    }
    return ljk;
}


/*  Findet Naechsten Nachbar von Teilchen r2 zu Teilchen r1
 *  Bei einem Cutoff von L/2 wechselwirkt Teilchen 1 bei r1 nur mit dem nächstgelegenen
 *  Teilchen 2, da alle anderen Bildteilchen von j weiter als L/2 enfternt liegen.
 *  INPUT       r1      Ortsvektor
 *              r2      Ortsvektor
 *              L       Kantenlaenge des Simulationsvolumens
 */
VectorXd naechsterNachbar(VectorXd r1, VectorXd r2, double L)
{
    double rc = L/2;
    VectorXd shift = VectorXd::Zero(2, 1);     // Verschiebungsvektor
    VectorXd diff = VectorXd::Zero(2, 1);   // Abstandsvektor zw. Teilchen und Bildteilchen

    // Pruefe zunaechst die Bildraeume

    // Bildraum rechts
    shift << L, 0;
    diff = r1-(r2+shift);
    if (diff.norm()<rc) {
        return r2+shift;
    }
    // Bildraum links
    diff = r1-(r2-shift);
    if (diff.norm()<rc) {
        return r2-shift;
    }

    // Bildraum oben
    shift << 0, L;
    diff = r1-(r2+shift);
    if (diff.norm()<rc) {
        return r2+shift;
    }
    // Bildraum unten
    diff = r1-(r2-shift);
    if (diff.norm()<rc) {
        return r2-shift;
    }

    // Bildraum oben rechts
    shift << L, L;
    diff = r1-(r2+shift);
    if (diff.norm()<rc) {
        return r2+shift;
    }
    // Bildraum unten links
    diff = r1-(r2-shift);
    if (diff.norm()<rc) {
        return r2-shift;
    }

    // Bildraum unten rechts
    shift << L, -L;
    diff = r1-(r2+shift);
    if (diff.norm()<rc) {
        return r2+shift;
    }
    // Bildraum oben links
    diff = r1-(r2-shift);
    if (diff.norm()<rc) {
        return r2-shift;
    }

    // Kein Bildteilchen nah genug, also muss das reale Teilchen am naechsten sein.
    return r2;
}


/*  Berechne wirkende Kraft auf jedes Teilchen
 *  INPUT       r       Ortsvektoren der Teilchen
 *              N       Anzahl der Teilchen
 *              F       (alte) Kraftvektoren der Teilchen
 *              L       Kantenlaenge des Gebiets [0,L]x[0,L]
 *  OUTPUT      F       wird ueberschrieben
 */
void update_kraft(MatrixXd &r, MatrixXd &F, unsigned int N, double L)
{
    double rc = L/2;
    F.setZero();
    // Folgende Groessen als Zwischenspeicher benoetigt
    VectorXd delta_r = VectorXd::Zero(2, 1);     // Abstandsvektor
    VectorXd r_next = VectorXd::Zero(2, 1);     // nächster Nachbar
    VectorXd F_temp = VectorXd::Zero(2, 1);     // Kraft zwischen j und i

    for (int i=0; i<(N-1); i++)
    {
        for (int j=i+1; j<N; j++)   // Keine Selbstwecheelwirkung oder Doppelzaehlung
        {
            r_next = naechsterNachbar(r.col(i), r.col(j), L);
            delta_r = r.col(i)-r_next;
            if(delta_r.norm()<rc)
            {
                F_temp = lj_kraft(delta_r, L);
                F.col(i) = F.col(i) + F_temp;
                F.col(j) = F.col(j) - F_temp;
            }
        }
    }
}


/*  Periodische Randbedingungen: Verhindere, dass Teilchen aus dem Simulationsgebiet
 *  heraus laufen
 *  EIN SIMULATIONSGEBIET, SIE ZU KNECHTEN, SIE ALLE ZU FINDEN,
 *  INS DUNKEL ZU TREIBEN UND EWIG ZU BINDEN.
 *  INPUT       r       2xN Matrix mit Ortsvektoren der Teilchen
 *              N       Anzahl der Teilchen
 *              L       Kantenlaenge des Gebiets [0,L]x[0,L]
 */
void periodische_RB(MatrixXd &r, unsigned int N, double L)
{
    // Der Modulo Operator ist besser als eine Addition, falls sich
    // Teilchen in einem Zeitschritt mehr als L entfernt haben
    // mit grossen Geschwindigkeiten durch hohe Temperaturen.
    // Leider macht er bei T=100 numerische Probleme :(
    for (int n=0; n<N; n++)     // Durchlaeuft alle Teilchen
    {
        // Teilchen nach links abgedriftet
        if(r(0,n) < 0) {
            // r(0,n) = fmod(r(0,n), L) + L;
            r(0,n) += L;
        }
        // Teilchen nach rechts abgedriftet
        if(r(0,n) > L) {
            // r(0,n) = fmod(r(0,n), L);
            r(0,n) -= L;
        }
        // Teilchen nach unten abgedriftet
        if(r(1,n) < 0) {
            // r(1,n) = fmod(r(1,n), L) + L;
            r(1,n) += L;
        }
        // Teilchen nach oben abgedriftet
        if(r(1,n) > L) {
            // r(1,n) = fmod(r(1,n), L);
            r(1,n) -= L;
        }
    }
}


/*  Initialisiere Molekulardynamik Simulation
 *  INPUT       N       Anzahl Teilchen
 *              L       Kantenlaenge der Box
 *              r       2xN-Matrix fuer Orte
 *              v       2xN-Matrix fuer Geschwindigkeiten
 *              Tinit   Starttemperatur
 *  OUTPUT      r und v werden ueberschrieben
 */
void md_init(unsigned int N, double L, MatrixXd &r, MatrixXd &v, double Tinit)
{
    double sqrtN = sqrt(N);
    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<> dis(0.0, 1.0);      // Zufallszahlen zwischen 0 und 1

    for (int n=0; n<sqrtN; n++)
    {
        for (int m=0; m<sqrtN; m++)
        {
            // Ordne Punkte auf einem Gitter an
            r(0, sqrtN*n+m) = 1./(2*sqrtN)*(1.+2.*n)*L;
            r(1, sqrtN*n+m) = 1./(2*sqrtN)*(1.+2.*m)*L;
            // Zufaellige Geschwindigkeiten
            v(0, sqrtN*n+m) = dis(gen);
            v(1, sqrtN*n+m) = dis(gen);
        }
    }

    // Nun Schwerpunktsbewegung auf Null setzen
    VectorXd vmean = 1./N * v.rowwise().sum();
    for (int n=0; n<N; n++)
    {
        v.col(n) = (v.col(n) - vmean);
    }

    // Geschwindigkeiten werden mit @skal umskaliert, um Temperatur auf Tinit zu setzen.
    double Nf = 2*N-2;                              // Anzahl Freiheitsgrade
    double skal = Tinit*Nf/v.colwise().squaredNorm().sum();
    for (int n=0; n<N; n++)
    {
        v.col(n) = sqrt(skal)*v.col(n);
    }
}


/*  Berechne potentielle Energie
 *  INPUT       r       2xN Matrix der Ortsvektoren
 *              N       Anzahl der Teilchen
 *              L       Kantenlaenge des Simulationsvolumens
 *  OUTPUT      Epot    Potentielle Energie der Konfiguration
 */
double potEnergie(MatrixXd &r, unsigned int N, double L)
{
    double rc = L/2;
    // Folgende Groessen als Zwischenspeicher benoetigt
    VectorXd delta_r = VectorXd::Zero(2, 1);     // Abstandsvektor
    VectorXd r_next = VectorXd::Zero(2, 1);
    double Epot = 0.;

    for (int i=0; i<(N-1); i++)
    {
        for (int j=i+1; j<N; j++)   // Keine Selbstwecheelwirkung oder Doppelzaehlung
        {
            r_next = naechsterNachbar(r.col(i), r.col(j), L);
            delta_r = r.col(i)-r_next;
            if(delta_r.norm()<rc)
            {
                // Lennard-Jones-Potential
                Epot += 4*(pow(delta_r.norm(), -12)-pow(delta_r.norm(), -6));
            }
        }
    }
    return Epot;
}


/*  Molekulardynamik-Simulation unter Verwendung des Verlet-Algorithmus
 *  INPUT       r       2xN Matrix der Ortsvektoren
 *              v       2xN Matrix der Geschwindigkeitsvektoren
 *              N       Anzahl der Teilchen
 *              L       Kantenlaenge des Simulationsvolumens
 *              h       Schrittweite in der Zeit
 *              Tinit   Temperatur zum Startzeitpunkt t=0
 *              tequi   Zeit, bis zu welcher Aequilibriert werden soll
 *              file_equi   file, in welchen Analyse der Aequilibrierung geschrieben wird
 *              equi_analyse    Energie und Temperaturbetrachtung der Aequilibrierung ja/nein
 *              tmax    Zeit, bis zu welcher Simuliert werden soll
 *              messung file, in welchem eigentliche Messung geschrieben wird
 *              thermostat  Thermostat angeschlossen -> Temperatur konstant?
 */
void md_simulation(MatrixXd &r,
                   MatrixXd &v,
                   unsigned int N,
                   double L,
                   double h,
                   double Tinit,
                   double tequi,
                   string file_equi,
                   bool equi_analyse,
                   double tmax,
                   string file_messung,
                   bool thermostat)
{
    double t = 0;                       // Aktuelle Zeit der Simulation

    cout << "Schritte Aequilibrierung:  " << tequi/h << endl;
    cout << "Schritte Simulation:       " << tmax/h << endl;

    ofstream stream;
    unsigned int Nf = 2*N-2;            // Anzahl Freiheitsgrade (in 2D)

    md_init(N, L, r, v, Tinit);

    MatrixXd F = MatrixXd::Zero(2, N);  // Matrix der Kraftvektoren
    update_kraft(r, F, N, L);
    MatrixXd Fold;  // vorherige Kraftmatrix zur Berechnung der neuen Geschwindigkeiten

    // AEQUILIBRIERUNG
    if (equi_analyse)
    {
        stream.open(file_equi, ios::trunc);
        stream.precision(10);
        stream << "t vsx vsy Ekin Epot Temp" << endl;
    }
    while (t < tequi)
    {
        if (equi_analyse)  // Analysiere Energien, Schwerpunktsbewegung und Temperatur
        {
            // Berechne Schwerpunktsgeschwindigkeit
            VectorXd vS = 1./N * v.rowwise().sum();
            stream << t << " " << vS(0) << " " << vS(1);
            // Berechne kinetische Energie
            double Ekin = 0.5*v.colwise().squaredNorm().sum();
            stream << " " << Ekin;
            // berechne potentielle Energie
            double Epot = potEnergie(r, N, L);
            stream << " " << Epot;
            // berechne Temperatur
            double T = 2*Ekin/Nf;
            stream << " " << T;
            stream << endl;
        }

        Fold = F;
        // Neue Ortsvektoren
        for (int i=0; i<N; i++)
        {
            r.col(i) = r.col(i)+v.col(i)*h+0.5*h*h*F.col(i);
        }
        // periodische RB: verhindere Herauslaufen aus Box
        periodische_RB(r, N, L);
        // Neue Kraefte
        update_kraft(r, F, N, L);
        // neue Geschwindigkeiten
        for (int i=0; i<N; i++)
        {
            v.col(i) = v.col(i)+0.5*h*(F.col(i)+Fold.col(i));
        }

        if (thermostat)
        {
            // Temperatur wird konstant gehalten,
            // dazu umskalierung der Geschwindigkeiten auf Tinit
            double skal = Tinit*Nf/v.colwise().squaredNorm().sum();
            for (int n=0; n<N; n++)
            {
                v.col(n) = sqrt(skal)*v.col(n);
            }
        }

        t += h;
    }
    if (equi_analyse)
    {
        stream.close();
    }


    // SIMULATION
    stream.open(file_messung, ios::trunc);
    stream.precision(10);
    stream << "t Temp" << endl;
    while (t < tmax+tequi)
    {
        // Messe Observablen des Systems
        stream << t-tequi;
        // berechne Temperatur
        double T = v.colwise().squaredNorm().sum()/Nf;
        stream << " " << T;
        // TODO: Paarkorrelationsfunktion implementieren und hier einfuegen
        stream << endl;

        Fold = F;
        // Neue Ortsvektoren
        for (int i=0; i<N; i++)
        {
            r.col(i) = r.col(i)+v.col(i)*h+0.5*h*h*F.col(i);
        }
        // periodische RB: verhindere Herauslaufen aus Box
        periodische_RB(r, N, L);
        // Neue Kraefte
        update_kraft(r, F, N, L);
        // neue Geschwindigkeiten
        for (int i=0; i<N; i++)
        {
            v.col(i) = v.col(i)+0.5*h*(F.col(i)+Fold.col(i));
        }

        if (thermostat)
        {
            // Temperatur wird konstant gehalten,
            // dazu umskalierung der Geschwindigkeiten auf Tinit
            double skal = Tinit*Nf/v.colwise().squaredNorm().sum();
            for (int n=0; n<N; n++)
            {
                v.col(n) = sqrt(skal)*v.col(n);
            }
        }

        t += h;
    }
    stream.close();
}


void save_data(MatrixXd &r, MatrixXd &v, ofstream &file)
{
    file << "x y vx vy" << endl;
    file.precision(10);
    for(int n=0; n<r.cols(); n++)
    {
        file << r(0,n) << " " << r(1,n) << " ";
        file << v(0,n) << " " << v(1,n) << endl;
    }
}


int main()
{
    cout << "Beginn des Programms!" << endl << endl;

    // Einstellbare Parameter
    unsigned int N = 16;    // Anzahl Teilchen
    double L = 8;           // Kantenlaenge der Box
    double h = 0.01;        // Zeitschritt des Verlet-Algorithmus
    double Tinit = 1;       // Temperatur bei Initialisierung t=0
    double tequi = 50;     // Äquilibrierungszeit
    double tmax = 500;      // Maximale Zeit der Simulation
    bool thermostat = false; // kein Thermostat angeschlossen, Temperatur variabel

    ofstream stream;
    MatrixXd r = MatrixXd::Zero(2, N);
    MatrixXd v = MatrixXd::Zero(2, N);

    // Teste LJ-Kraft
    // cout << lj_kraft(r.col(0)-r.col(1), L/2) << endl;

    // Teste Berechnung des Kraftfeldes
    // MatrixXd F = MatrixXd::Zero(2,N);
    // update_kraft(r, F, N, L);
    // cout << F << endl;

    // Verifikation der Init-Funktion
    md_init(N, L, r, v, Tinit);
    stream.open("build/init.txt", ios::trunc);
    save_data(r, v, stream);
    stream.close();

    // Aufgabenteil b) Betrachte Äquilibrierungsphase
    cout << "MD-Simulation 1 - T=1" << endl;
    md_simulation(r, v, N, L, h, Tinit,
            tequi, "build/equilibration.txt", true,
            tmax, "build/messung_T1.txt",
            thermostat);
    cout << endl;

    // Aufgabenteil c) Betrachte T und Paarkorrelationsfunktion
    Tinit = 100;
    // passe bei so hohen Temperaturen die Schrittweite an
    h = 0.001;
    tequi = 5;
    tmax = 50;
    cout << "MD-Simulation 2 - T=100" << endl;
    md_simulation(r, v, N, L, h, Tinit,
            tequi, "build/dummy.txt", false,
            tmax, "build/messung_T1e2.txt",
            thermostat);
    // cout << r << endl << v << endl;
    cout << endl;
    h = 0.01;
    tequi = 50;
    tmax = 500;

    Tinit = 0.01;
    cout << "MD-Simulation 3 - T=0.01" << endl;
    md_simulation(r, v, N, L, h, Tinit,
            tequi, "build/dummy.txt", false,
            tmax, "build/messung_T1e-2.txt",
            thermostat);
    cout << endl;


    // Aufgabenteil d) Thermostat
    thermostat = true;
    cout << "MD-Simulation 4 - Thermostat eingebaut" << endl;
    md_simulation(r, v, N, L, h, Tinit,
            tequi, "build/thermostat_equi.txt", true,
            tmax, "build/thermostat_mess.txt",
            thermostat);
    cout << endl;
    // Konfiguration am letzten Zeitschritt
    stream.open("build/endkonfig.txt", ios::trunc);
    save_data(r, v, stream);
    stream.close();

    cout << "Ende des Programms!" << endl;
    return 0;
}
