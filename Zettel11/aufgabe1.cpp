#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <math.h>  // sqrt()
#include <random>

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
double weglaenge(MatrixXd &r, VectorXi &perm)
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
    double T = Tstart; // Temperatur des Systems
    double ort_1; // erster vorgeschlagener Ortsvektor zum tauschen
    double ort_2; // zweiter vorgeschlagener Ortsvektor zum tauschen
    double zwischenspeicher;   // zwischenspeicher zum vertauschen
    VectorXi sugg_perm = perm;  // jeweils vorgeschlagene Permutation
    double differenz;   // Weglaengendifferenz

    // Ziehe Ortsvektor-index mit ziehe_index()
    std::default_random_engine gen_index;
    std::uniform_int_distribution<int> indices(0,r.cols()-1);
    auto ziehe_index = std::bind(indices, gen_index);
    // Ziehe Zufallszahl in [0,1]
    std::default_random_engine gen_number;
    std::uniform_real_distribution<double> numbers(0.,1.);
    auto ziehe_zahl = std::bind(numbers, gen_number);

    while(T >= Tend)
    {
        for(int i=0; i<S; i++)  // Schlage MC-Schritt vor
        {
            sugg_perm = perm;

            // Waehle Vertauschung
            ort_1 = ziehe_index();
            ort_2 = ziehe_index();
            zwischenspeicher = sugg_perm(ort_1);
            sugg_perm(ort_1) = sugg_perm(ort_2);
            sugg_perm(ort_2) = zwischenspeicher;
            differenz = weglaenge(r, sugg_perm) - weglaenge(r, perm);

            // Pruefe Akzeptanz der Vertauschung
            if (differenz < 0)
            {
                perm = sugg_perm;
                opti_perm = perm;
            }
            else    // Weglaenge durch Vertauschung verlaengert
            {
                if(ziehe_zahl() < exp(-differenz/T))
                {
                    perm = sugg_perm;
                }
            }
            // cout << weglaenge(r, perm) << endl;
        }
        T *= d; // Absenkung der Temperatur
        cout << T << endl;
    }

    cout << endl << weglaenge(r, opti_perm) << endl;
    cout << perm.transpose() << endl;

    // Speichere die beste gefundene Permutation
    save_data(r, opti_perm, file);
}


int main() {
    cout << "Beginn des Programms!\n" << endl;

    // Einstellbare Parameter
    double Tstart = 10.;         // initiale Temperatur
    double Tend = 1e-2;           // ungefaehre Endtemperatur
    double d = 0.9;             // Daempfungsfaktor
    unsigned int S = 10000;             // angebotene Vertauschungen
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

    // random_device rd;
    // mt19937 generator(rd());
    // uniform_int_distribution<int> distribution(1,6);
    // auto randomNumber = distribution(generator);

    stream.open("build/test.txt");
    stimulated_annealing(r, perm, Tstart, Tend, d, S, stream);
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
