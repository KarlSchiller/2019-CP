#include <iostream>
#include <fstream>
#include <Eigen/Dense>
// #include <math.h>  // sqrt, round
// #include <complex>
#include <iomanip> // setprecision
#include <random>

using namespace std;
using namespace Eigen;


double h(random_device &rd){
  mt19937 generator(rd());
  uniform_real_distribution<double> distribution(-5,5);
  auto randomNumber = distribution(generator);
  return randomNumber;
}

double hamilton(int sigma, double magnet){
  return -sigma*magnet;
}

// Funktion für einen MC-Schritt, gibt 1 oder -1 zurück
int mc(random_device &rd, int alter_Schritt, double magnet){
  double zwischen = 0;
  int neuer_Schritt = -alter_Schritt;
  bool accept = false;


  mt19937 generator(rd());
  uniform_int_distribution<int> verteilung(0,1);
  auto p = verteilung(generator);

  // MC-Move anbieten
    zwischen = hamilton(neuer_Schritt, magnet) - hamilton(alter_Schritt, magnet);
    if (zwischen < 0){
      // cout << "Zwischen < 0" << endl;
      accept = true;
    }
    else{
      // cout << "else" << endl;
      p = verteilung(generator);
      // cout << "p: " << p << endl;
      zwischen = exp(-(hamilton(neuer_Schritt, magnet) - hamilton(alter_Schritt, magnet)));
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
    return neuer_Schritt;
  }
  else{
    // cout << "alter_Schritt" << endl;
    return alter_Schritt;
  }
}

int main() {
  cout << "Beginn des Programms!\n" << endl;
  random_device rd;
  double schritte = 1e5, alter_Schritt = 1, neuer_Schritt;
  ofstream file;

  VectorXd spins, magnetfeld;
  // Initialisieren des Magnetfeld-Vektors
  magnetfeld = VectorXd::LinSpaced(100, -5, 5);

  file.open("build/aufg2.txt", ios::trunc);
  // file << alter_Schritt << " ";
  for (int j = 0; j < magnetfeld.size(); j++){
    for (int i = 0; i < schritte; i++){
      neuer_Schritt = mc(rd, alter_Schritt, magnetfeld(j));
      file << neuer_Schritt << " ";
      alter_Schritt = neuer_Schritt;
    }
    file << endl;
    cout << j << endl;
    alter_Schritt = 1;
    // file << alter_Schritt << " ";
  }
  file.close();

  // cout << mc(rd, -1, -2) << endl;
  cout << "\nEnde des Programms!" << endl;
  return 0;
}
