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

inline double hamilton(int sigma, double magnet){
  return -sigma*magnet;
}

// Funktion für einen MC-Schritt, gibt 1 oder -1 zurück
void mc(int alter_Schritt, double magnet, mt19937 &generator,
       uniform_real_distribution<double> &verteilung, ofstream &file){
  double zwischen = 0;
  int neuer_Schritt;
  double count_pos = 0, count_neg = 0;
  bool accept = false;
  auto p = verteilung(generator);


  // mt19937 generator(rd());
  for (int i = 0; i < 1e5; i++){
    p = verteilung(generator);
    neuer_Schritt = -alter_Schritt;

    // MC-Move anbieten
      // zwischen = hamilton(neuer_Schritt, magnet) - hamilton(alter_Schritt, magnet);
      zwischen = -neuer_Schritt*magnet - -alter_Schritt*magnet;
      if (zwischen < 0){
        // cout << "Zwischen < 0" << endl;
        accept = true;
      }
      else{
        // cout << "else" << endl;
        // p = verteilung(generator);
        // cout << "p: " << p << endl;
        // hamilton(neuer_Schritt, magnet) - hamilton(alter_Schritt, magnet
        zwischen = exp(-(zwischen));
        // zwischen = 1;
        if (p < zwischen){
          // cout << "p < zwischen" << endl;
          accept = true;
        }
        else{
          // cout << "else" << endl;
          accept = false;
        }
      }

    // cout << "accept: " << accept << endl;
    // accept überprüfen
    if (accept == true){
      // return neuer_Schritt;
      // if (neuer_Schritt > 0){
      //   count_pos ++;
      // }
      // else count_neg --;
      file << neuer_Schritt << " ";
      alter_Schritt = neuer_Schritt;
    }
    else{
      // if (alter_Schritt > 0){
      //   count_pos ++;
      // }
      // else count_neg --;
      // cout << "alter_Schritt" << endl;
      // return alter_Schritt;
      file << alter_Schritt << " ";
    }
  }

}

int main() {
  cout << "Beginn des Programms!\n" << endl;
  random_device rd;
  double schritte = 1e5, alter_Schritt = 1;
  ofstream file;

  mt19937 generator(rd());
  uniform_real_distribution<double> distribution(0,1);

  VectorXd spins, magnetfeld;
  // Initialisieren des Magnetfeld-Vektors
  magnetfeld = VectorXd::LinSpaced(1e2, -5, 5);

  file.open("build/aufg2.txt", ios::trunc);
  // file << alter_Schritt << " ";
  for (int j = 0; j < magnetfeld.size(); j++){
      mc(alter_Schritt, magnetfeld(j), generator, distribution, file);
      file << endl;
      cout << j << endl;
      alter_Schritt = 1;
  }
  file.close();

  cout << "\nEnde des Programms!" << endl;
  return 0;
}
