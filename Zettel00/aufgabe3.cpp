#include <iostream>
#include <math.h> // provides exp
#include <fstream> // write to file
using namespace std;


void euler_normal(int n, double y_0, double delta_t, double* y){
  y[0] = y_0;
  for(int i=1; i<=n; i++){
    y[i] = y[i-1]*(1-delta_t);
  }
}


void euler_symm(int n, double y_0, double delta_t, double* y){
  y[0] = y_0;
  y[1] = exp(-delta_t);
  for(int i=2; i<=n; i++){
    y[i] = -2*delta_t*y[i-1]+y[i-2];
  }
}


int main(){

  cout << "Beginn des Programms!" << endl;
  double y_0 = 1.0;
  double delta_t = 0.2;
  int n = 10/delta_t;
  cout << "n=" << n << endl;
  double y_normal[n];
  euler_normal(n, y_0, delta_t, y_normal);
  double y_symm[n];
  euler_symm(n, y_0, delta_t, y_symm);

  ofstream myfile;
  myfile.open("euler_vergleich.txt", ios::trunc);
  myfile << "# n*delta_t euler-normal euler-symm analytisch" << endl;
  myfile << "t;normal;symm;ana" << endl;
  for(int i=0; i<=n; i++){
    myfile << i*delta_t << "; ";
    myfile << y_normal[i] << "; ";
    myfile << y_symm[i] << "; ";
    myfile << exp(-i*delta_t) << endl;
  }
  myfile.close();

  cout << "Ende des Programms!" << endl;
  return 0;
}
