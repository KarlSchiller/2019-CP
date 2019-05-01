#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>
#include <iomanip>

using namespace std;
using namespace Eigen;

// Definition der Funktionen
double funk_1(double x){
  return exp(-x)/x;
}

double funk_2(double x){
  if(x == 0){
    return 0;
  }
  else {
    return x*sin(1/x);
  }
}

// Definition Trapez
double trapez(double (*funptr)(double), double a, double b, double N){
  double h = (b-a)/N;
  double res=0, xk=0;
  for(int i = 1; i<N; i++){
    xk = a + i*h;
    res = res + funptr(xk);
  }
  res = res*h;
  res = res + h/2*(funptr(a) + funptr(b));
  return res;
}

// Mittelpunktsregel
double mittel(double (*funptr)(double), double a, double b, double N){
  double h = (b-a)/N;
  double res=0;
  for(int i = 1; i<N+1; i++){
    res = res + funptr(a-h/2+i*h);
  }
  res *= h;
  return res;
}

// Simpsonregel
double simpson(double (*funptr)(double), double a, double b, int N){
  double h = (b-a)/N;
  double n = N/2;
  double res=0, x2k=0, x2k1 = 0;
  if(fmod(N, 2) ==0){
    for(int i = 1; i<n; i++){
      x2k = a + i*2*h;
      res += 2*funptr(x2k);
    }
    for(int i = 1; i<n+1; i++){
        x2k1 = a + (2*i-1)*h;
        res += 4*funptr(x2k1);
    }
    res += funptr(a) + funptr(b);
    res *= h/3;
  }
  return res;
}


int main()
{
  cout << "Beginn des Programms!" << endl;
  // Initialisieren der Größen
  double a_2=0.0, a_1=1.0;
  double b_2 =1.0, b_1=100.0;
  double N = 2.0;
  double res1_1 = 0, res2_1=0, res3_1 = 0, abweichung_1=0, zwischen_1;
  double res1_2 = 0, res2_2=0, res3_2 = 0, abweichung_2=0, zwischen_2;
  // Erster Durchlauf mit Simpsonsregel
  res3_1 = simpson(funk_1, a_1, b_1, N);
  zwischen_1 = res3_1;
  res3_2 = simpson(funk_2, a_2, b_2, N);
  zwischen_2 = res3_2;

  // Zweiter Durchlauf mit Simpsonregel
  // Funktion 1
  do{
    N = N*2;
    res3_1 = simpson(funk_1, a_1, b_1, N);

    abweichung_1 = 1 - res3_1/zwischen_1;

    zwischen_1 = res3_1;
  }while(abs(abweichung_1) >= 1e-4);

  // Funktion 2
  do{
    N = N*2;
    res3_2 = simpson(funk_2, a_2, b_2, N);

    abweichung_2 = 1 - res3_2/zwischen_2;

    zwischen_2 = res3_2;
  }while(abs(abweichung_2) >= 1e-4);

  // Erster Durchlauf Mittelpunktsregel
  //N = 2.0;
  res2_1 = mittel(funk_1, a_1, b_1, N);
  zwischen_1 = res2_1;
  res2_2 = mittel(funk_2, a_2, b_2, N);
  zwischen_2 = res2_2;

  // Zweiter Durchlauf Mittelpunktsregel
  // Funktion 1
  do{
    N = N*2;
    res2_1 = mittel(funk_1, a_1, b_1, N);

    abweichung_1 = 1 - res2_1/zwischen_1;

    zwischen_1 = res2_1;
  }while(abs(abweichung_1) >= 1e-4);

  // Funktion 2
  do{
    N = N*2;
    res2_2 = mittel(funk_2, a_2, b_2, N);

    abweichung_2 = 1 - res2_2/zwischen_2;

    zwischen_2 = res2_2;
  }while(abs(abweichung_2) >= 1e-4);

  // Erster Durchlauf mit Trapezregel
  //N = 2.0;
  res1_1 = trapez(funk_1, a_1, b_1, N);
  zwischen_1 = res1_1;
  res1_2 = trapez(funk_2, a_2, b_2, N);
  zwischen_2 = res1_2;

  // Zweiter Durchlauf mit Trapezregel
  // Funktion 1
  do{
    N = N*2;
    res1_1 = trapez(funk_1, a_1, b_1, N);

    abweichung_1 = 1 - res1_1/zwischen_1;

    zwischen_1 = res1_1;
  }while(abs(abweichung_1) >= 1e-4);

  // Funktion 2
  do{
    N = N*2;
    res1_2 = trapez(funk_2, a_2, b_2, N);

    abweichung_2 = 1 - res1_2/zwischen_2;

    zwischen_2 = res1_2;
  }while(abs(abweichung_2) >= 1e-4);

  // Ausgabe der Ergebnisse
  cout << "Erste Funktion: " << endl;
  cout << "Trapez: " << setprecision(8) <<  res1_1 << endl;
  cout << "Mittel: " << setprecision(8) <<  res2_1 << endl;
  cout << "Simpson: " << setprecision(8) <<  res3_1 << endl;

  cout << endl;
  cout << "Zweite Funktion: " << endl;
  cout << "Trapez: " << setprecision(8) <<  res1_2 << endl;
  cout << "Mittel: " << setprecision(8) <<  res2_2 << endl;
  cout << "Simpson: " << setprecision(8) <<  res3_2 << endl;

  cout << "Ende des Programms!" << endl;
  return 0;
}
