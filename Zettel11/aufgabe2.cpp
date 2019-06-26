#include <iostream>
#include <fstream>
#include <Eigen/Dense>
// #include <math.h>  // sqrt, round
// #include <complex>
#include <iomanip> // setprecision
#include <random>

using namespace std;
using namespace Eigen;


// double h(random_device &rd){
//   mt19937 generator(rd());
//   uniform_real_distribution<double> distribution(-5,5);
//   auto randomNumber = distribution(generator);
//   return randomNumber;
// }

// inline double hamilton(int sigma, double magnet){
//   return -sigma*magnet;
// }

/* Funktion für MC-Simulation für eine belibige Anzahl an Schritten

  alter_Schritt Startwert
  magnet        Magnetfeldstärke
  generator     Zufallszahlgenerator
  verteilung    Verteilung, aus der gesampled wird
  file          output file
  Schritten     Anzahl der Schritte in der MC-Simulation
*/
void mc(int alter_Schritt, double magnet, mt19937 &generator,
       uniform_real_distribution<double> &verteilung, ofstream &file,
       double schritte){
  double zwischen = 0;
  int neuer_Schritt;
  double count_pos = 0, count_neg = 0;
  bool accept = false;
  auto p = verteilung(generator);


  // mt19937 generator(rd());
  for (int i = 0; i < schritte; i++){
    p = verteilung(generator);
    neuer_Schritt = -alter_Schritt;

    // MC-Move anbieten
      zwischen = -neuer_Schritt*magnet - -alter_Schritt*magnet;
      if (zwischen < 0){
        // cout << "Zwischen < 0" << endl;
        accept = true;
      }
      else{
        zwischen = exp(-zwischen);
        // zwischen = 1;
        if (p < zwischen){
          accept = true;
        }
        else{
          accept = false;
        }
      }

    // accept überprüfen
    if (accept == true){
      if (neuer_Schritt > 0){
        // Falls +1, counter für +1 hochsetzen
        count_pos ++;
      }
      // Falls -1, counter für -1 verkleinern
      else count_neg --;
      alter_Schritt = neuer_Schritt;
    }
    else{
      if (alter_Schritt > 0){
        // Falls +1, counter für +1 hochsetzen
        count_pos ++;
      }
      // Falls -1, counter für -1 verkleinern
      else count_neg --;
    }
  }
  // Differenz in ein file schreiben
  file << count_pos + count_neg;
}

int main() {
  cout << "Beginn des Programms!\n" << endl;
  random_device rd;
  double alter_Schritt = 1;
  ofstream file;
  double schritte = 1e5;

  mt19937 generator(rd());
  uniform_real_distribution<double> distribution(0,1);

  VectorXd spins, magnetfeld, magnetfeld100;
  // Initialisieren des Magnetfeld-Vektors
  magnetfeld = VectorXd::LinSpaced(1e4, -5, 5);
  magnetfeld100 = VectorXd::LinSpaced(1e2, -5, 5);

  file.open("build/aufg2_100.txt", ios::trunc);
  for (int j = 0; j < magnetfeld100.size(); j++){
      mc(alter_Schritt, magnetfeld100(j), generator, distribution, file, schritte);
      file << endl;
      alter_Schritt = 1;
  }
  file.close();

  file.open("build/aufg2.txt", ios::trunc);
  for (int j = 0; j < magnetfeld.size(); j++){
      mc(alter_Schritt, magnetfeld(j), generator, distribution, file, schritte);
      file << endl;
      if(j%1000==0)
      {
          cout << j << endl;
      }
      alter_Schritt = 1;
  }
  file.close();

  cout << "\nEnde des Programms!" << endl;
  return 0;
}
