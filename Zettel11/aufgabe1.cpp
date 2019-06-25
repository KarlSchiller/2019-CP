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


/* Initialisierung Testumgebung Simulated Annealing
 * INPUT        r       Matrix der Ortsvektoren
 *              delta   Abstand der Ortsvektoren
 */
void sa_init(MatrixXd &r, double delta)
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

    // Starte unten links und gehe gegen den Uhrzeigersinn durch
    steps = int(1./delta - 1);
    r.block(0, 8, 1, steps) = VectorXd::LinSpaced(steps, 0+delta, 1-delta).transpose();
    index = 8+steps;

    steps = int(2./delta - 1);
    r.block(0, index, 1, steps) = VectorXd::Ones(steps).transpose();
    r.block(1, index, 1, steps) = VectorXd::LinSpaced(steps, 0+delta, 2-delta).transpose();
    index += steps;
    cout << r << endl;
    // TODO: hier weiter

    // bestimme initiale Permutation
    // TODO
}


/*  Speichere die aktuelle Konfiguration
 *  INPUT       r       Matrix der Ortsvektoren
 *              file    file, in welchen gespeichert wird
 */
void save_data(MatrixXd &r, ofstream &file)
{
    file << "x y" << endl;
    file.precision(10);
    for(int n=0; n<r.cols(); n++)
    {
        file << r(0,n) << " " << r(1,n);
        file << endl;
    }
}


int main() {
    cout << "Beginn des Programms!\n" << endl;

    // Einstellbare Parameter
    double Tstart = 0.;         // initiale Temperatur
    double Tend = 0.;           // ungefaehre Endtemperatur
    double d = 0.1;             // Daempfungsfaktor
    double S = 100;             // angebotene Vertauschungen
    double delta = 0.2;         // Abstand der Ortsvektoren

    // Initialisierung der Testumgebung
    ofstream stream;
    MatrixXd r;
    sa_init(r, delta);
    stream.open("build/init.txt");
    save_data(r, stream);
    stream.close();

    // random_device rd;
    // mt19937 generator(rd());
    // uniform_int_distribution<int> distribution(1,6);
    // auto randomNumber = distribution(generator);

    // Test der Funktion sa_number_of_pos_vecs
    // cout << "Hier sollte 90 stehen: " << sa_number_of_pos_vecs(0.2) << endl;

    cout << "\nEnde des Programms!" << endl;
    return 0;
}
