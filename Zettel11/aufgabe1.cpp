#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <math.h>  // sqrt()
#include <random>
#include <chrono>  // get elapsed time

using namespace std;
using namespace Eigen;


/*  Hilfsmethode zur Testumgebung Simulated Annealing
 *  INPUT       delta   Abstand der Ortsvektoren
 *  OUTPUT      N       Anzahl der Ortsvektoren
 */
unsigned int sa_number_of_pos_vecs(double delta)
{
    unsigned int N = 8;     // Eckpunkte
    N += 2*(3/delta -1);    // aeussere vertikale Begrenzungen
    N += 2*(2/delta -1);    // innere vertikale Begrenzungen
    N += 4/delta -1;        // oberste horizontale Begrenzung
    N += 2*(1/delta -1);    // unterste vertikale Begrenzungen
    N += 2/delta -1;        // mittlere vertikale Begrenzung
    return N;
}


/*  Berechne Weglaenge einer Permutation
 *  INPUT       r       Matrix der Ortsvektoren
 *              perm    Permutation
 *  OUTPUT      L       Weglaenge
 */
inline double weglaenge(MatrixXd &r, VectorXi &perm)
{
    double L = 0;
    double N = r.cols();
    for(int i=0; i<N-1; i++)  // bis zum vorletzten Element
    {
        L += (r.col(perm(i))-r.col(perm(i+1))).norm();
    }
    // Ende der Kette, verbinde mit Start
    L += (r.col(0)-r.col(N-1)).norm();
    return L;
}


/*  Initialisierung Testumgebung Simulated Annealing
 *  INPUT       r       Matrix der Ortsvektoren
 *              perm    Permutation
 *              delta   Abstand der Ortsvektoren
 *  OUTPUT      r und perm werden ueberschrieben
 */
void sa_init(MatrixXd &r, VectorXi &perm, double delta)
{
    // benoetigte Hilfsvariablen
    int steps;      // Anzahl Ortsvektoren fuer den jeweiligen Schritt
    int index;      // index der aktuellen Spalte von r

    // Erstelle Matrix mit Ortsvektoren
    unsigned int N = sa_number_of_pos_vecs(delta);
    r = MatrixXd::Zero(2, N);

    // Eckvektoren
    r.col(1) << 1,0;
    r.col(2) << 1,2;
    r.col(3) << 3,2;
    r.col(4) << 3,0;
    r.col(5) << 4,0;
    r.col(6) << 4,3;
    r.col(7) << 0,3;
    index = 8;

    // Starte unten links und gehe gegen den Uhrzeigersinn durch
    // nach rechts
    steps = int(1./delta - 1);
    r.block(0, index, 1, steps) = VectorXd::LinSpaced(steps, 0+delta, 1-delta).transpose();
    index = 8+steps;

    // nach oben
    steps = int(2./delta - 1);
    r.block(0, index, 1, steps) = VectorXd::Ones(steps).transpose();
    r.block(1, index, 1, steps) = VectorXd::LinSpaced(steps, 0+delta, 2-delta).transpose();
    index += steps;

    // nach rechts
    r.block(0, index, 1, steps) = VectorXd::LinSpaced(steps, 1+delta, 3-delta).transpose();
    r.block(1, index, 1, steps) = 2.*VectorXd::Ones(steps).transpose();
    index += steps;

    // nach unten
    r.block(0, index, 1, steps) = 3.*VectorXd::Ones(steps).transpose();
    r.block(1, index, 1, steps) = VectorXd::LinSpaced(steps, 2-delta, 0+delta).transpose();
    index += steps;

    // nach rechts
    steps = int(1./delta - 1);
    r.block(0, index, 1, steps) = VectorXd::LinSpaced(steps, 3+delta, 4-delta).transpose();
    index += steps;

    // nach oben
    steps = int(3./delta - 1);
    r.block(0, index, 1, steps) = 4.*VectorXd::Ones(steps).transpose();
    r.block(1, index, 1, steps) = VectorXd::LinSpaced(steps, 0+delta, 3-delta).transpose();
    index += steps;

    // nach links
    steps = int(4./delta - 1);
    r.block(0, index, 1, steps) = VectorXd::LinSpaced(steps, 4-delta, 0+delta).transpose();
    r.block(1, index, 1, steps) = 3.*VectorXd::Ones(steps).transpose();
    index += steps;

    // nach unten
    steps = int(3./delta - 1);
    r.block(1, index, 1, steps) = VectorXd::LinSpaced(steps, 3-delta, 0+delta).transpose();

    // bestimme initiale Permutation
    perm = VectorXi::LinSpaced(N, 0, N-1);
    std::random_shuffle(perm.data(), perm.data()+perm.size());
}


/*  Speichere die aktuelle Konfiguration
 *  INPUT       r       Matrix der Ortsvektoren
 *              file    file, in welchen gespeichert wird
 */
void save_data(MatrixXd &r, VectorXi &perm, ofstream &file)
{
    file << "x y perm" << endl;
    file.precision(10);
    for(int n=0; n<r.cols(); n++)
    {
        file << r(0,n) << " " << r(1,n) << " " << perm(n);
        file << endl;
    }
}


/*  Stimulated Annealing / Travelling Salesman Problem
 *  INPUT       r       Matrix der Ortsvektoren
 *              perm    Permutation
 *              Tstart  Initiale Temperatur
 *              Tend    ungefaehre Endtemperatur
 *              d       Daempfungsfaktor, mit welchem T multipliziert wird
 *              S       Anzahl angebotener Vertauschungen vor Verringerung von T
 *              file    file, in welchen optimale Konfiguration gespeichert wird
 */
void stimulated_annealing(
        MatrixXd &r,
        VectorXi &perm,
        double Tstart,
        double Tend,
        double d,
        unsigned int S,
        ofstream &file)
{
    // Benoetigte Hilfsvariablen
    VectorXi opti_perm = perm;  // beste ermittelte Permutation
    VectorXi sugg_perm = perm;  // jeweils vorgeschlagene Permutation
    double T = Tstart; // Temperatur des Systems
    double weg_alt = weglaenge(r, perm);     // Weglaengendifferenz des vorherigen Schrittes
    double weg_sugg;    // vorgeschlagene Weglaenge
    double weg_diff;   // Weglaengendifferenz
    double weg_opti = weg_alt;    // bisherige optimale Weglaenge
    unsigned int counter = 0;  // Anzahl Vertauschungen insgesamt

    // Ziehe Ortsvektor-index mit ziehe_index()
    std::default_random_engine gen_index;
    std::uniform_int_distribution<int> indices(0,r.cols()-1);
    auto ziehe_index = std::bind(indices, gen_index);
    // Ziehe Zufallszahl in [0,1]
    std::default_random_engine gen_number;
    std::uniform_real_distribution<double> numbers(0.,1.);
    auto ziehe_zahl = std::bind(numbers, gen_number);

    auto start = std::chrono::system_clock::now();
    while(T >= Tend)
    {
        for(int i=0; i<S; i++)
        {
            sugg_perm = perm;

            // Waehle Vertauschung
            sugg_perm.row(ziehe_index()).swap(sugg_perm.row(ziehe_index()));
            weg_sugg = weglaenge(r, sugg_perm);

            // Schlage Vertauschung vor
            // differenz = weglaenge(r, sugg_perm) - weglaenge(r, perm);
            weg_diff = weg_sugg - weg_alt;
            if(weg_diff < 0)
            {
                perm = sugg_perm;
                weg_alt = weg_sugg;

                if(weg_sugg < weg_opti)    // Ist Neuer Weg auch das Optimum?
                {
                    opti_perm = perm;
                    weg_opti = weg_sugg;
                }
            }
            else    // Weglaenge durch Vertauschung verlaengert
            {
                if(ziehe_zahl() < exp(-weg_diff/T))
                {
                    perm = sugg_perm;
                    weg_alt = weg_sugg;
                }
            }
        }
        counter += S;
        T *= d; // Absenkung der Temperatur
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    cout << "Benoetigte Zeit: " << elapsed_seconds.count() << "s" << endl;
    cout << "Anzahl Angebotener Vertauschungen: " << counter << endl;
    cout << "Beste gefundene Weglaenge: " << weg_opti << endl;

    // Speichere die beste gefundene Permutation
    save_data(r, opti_perm, file);
}


int main() {
    cout << "Beginn des Programms!\n" << endl;

    // Einstellbare Parameter
    double Tstart = 10.;         // initiale Temperatur
    double Tend = 1e-2;           // ungefaehre Endtemperatur
    double d = 0.99;             // Daempfungsfaktor
    unsigned int S = 100;             // angebotene Vertauschungen
    double delta = 0.2;         // Abstand der Ortsvektoren

    cout << "Anzahl an Ortsvektoren:    " << sa_number_of_pos_vecs(delta) << endl;

    // Initialisierung der Testumgebung
    ofstream stream;
    MatrixXd r;
    VectorXi perm;
    sa_init(r, perm, delta);
    stream.open("build/init.txt");
    save_data(r, perm, stream);
    stream.close();

    // d=0.9    S=100
    cout << endl << "d=0,9      S=1e2" << endl;
    stream.open("build/d9S1e2.txt");
    stimulated_annealing(r, perm, Tstart, Tend, 0.9, 100, stream);
    stream.close();

    // d=0.9    S=1000
    cout << endl << "d=0,9      S=1e3" << endl;
    stream.open("build/d9S1e3.txt");
    stimulated_annealing(r, perm, Tstart, Tend, 0.9, 1000, stream);
    stream.close();

    // d=0.9    S=10000
    cout << endl << "d=0,9      S=1e4" << endl;
    stream.open("build/d9S1e4.txt");
    stimulated_annealing(r, perm, Tstart, Tend, 0.9, 10000, stream);
    stream.close();

    // d=0.99   S=10
    cout << endl << "d=0,99     S=1e1" << endl;
    stream.open("build/d99S1e1.txt");
    stimulated_annealing(r, perm, Tstart, Tend, 0.99, 10, stream);
    stream.close();

    // d=0.99   S=100
    cout << endl << "d=0,99     S=1e2" << endl;
    stream.open("build/d99S1e2.txt");
    stimulated_annealing(r, perm, Tstart, Tend, 0.99, 100, stream);
    stream.close();

    // d=0.99   S=1000
    cout << endl << "d=0,99     S=1e3" << endl;
    stream.open("build/d99S1e3.txt");
    stimulated_annealing(r, perm, Tstart, Tend, 0.99, 1000, stream);
    stream.close();

    // d=0.999  S=10
    cout << endl << "d=0,999    S=1e1" << endl;
    stream.open("build/d999S1e1.txt");
    stimulated_annealing(r, perm, Tstart, Tend, 0.999, 10, stream);
    stream.close();

    // d=0.999  S=100
    cout << endl << "d=0,999    S=1e2" << endl;
    stream.open("build/d999S1e2.txt");
    stimulated_annealing(r, perm, Tstart, Tend, 0.999, 100, stream);
    stream.close();

    // Test der Funktion sa_number_of_pos_vecs
    // cout << "Hier sollte 90 stehen: " << sa_number_of_pos_vecs(0.2) << endl;

    // Test der Berechnung der Weglaenge
    // cout << endl << "Test der Weglaengenberechnung" << endl;
    // MatrixXd test(2,3);
    // test << 0, 1, 1, -3, 0, 2;
    // cout << "Matrix:" << endl << test << endl;
    // VectorXi test_perm(3);
    // test_perm << 0, 1, 2;
    // cout << "Permutation: " << test_perm.transpose() << endl;
    // cout << "Wert der Funktion (sollte 10.26129717 sein): " << weglaenge(test, test_perm) << endl;

    cout << "\nEnde des Programms!" << endl;
    return 0;
}
