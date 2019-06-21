#include <iostream>
#include <fstream>
#include <Eigen/Dense>
// #include <math.h>  // sqrt, round
// #include <complex>
#include <iomanip> // setprecision

using namespace std;
using namespace Eigen;

int main() {
  cout << "Beginn des Programms!\n" << endl;
  random_device rd;
  mt19937 generator(rd());
  uniform_int_distribution<int> distribution(1,6);
  for (int i = 0; i<10; i++){
    auto randomNumber = distribution(generator);
    cout << randomNumber << endl;
  }
  cout << "\nEnde des Programms!" << endl;
  return 0;
}
