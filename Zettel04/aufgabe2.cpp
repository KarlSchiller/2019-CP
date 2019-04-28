#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>

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
  double a=0.0;
  double b =1.0;
  double N = 2.0;
  double res1 = 0, res2=0, res3 = 0, abweichung=0, zwischen;
  // Erster Durchlauf
  res3 = simpson(funk_2, a, b, N);
  zwischen = res3;

  // Zweiter Durchlauf
  do{
    N = N*2;
    res3 = simpson(funk_2, a, b, N);
    abweichung = 1 - res3/zwischen;
    zwischen = res3;
  }while(abs(abweichung) >= 1e-4);

  // Erster Durchlauf
  res2 = mittel(funk_2, a, b, N);
  zwischen = res2;

  // Zweiter Durchlauf
  do{
    N = N*2;
    res2 = mittel(funk_2, a, b, N);
    abweichung = 1 - res2/zwischen;
    zwischen = res2;
  }while(abs(abweichung) >= 1e-4);

  // Erster Durchlauf
  res1 = trapez(funk_2, a, b, N);
  zwischen = res1;

  // Zweiter Durchlauf
  do{
    N = N*2;
    res1 = trapez(funk_2, a, b, N);
    abweichung = 1 - res1/zwischen;
    zwischen = res1;
  }while(abs(abweichung) >= 1e-4);

  cout << "Trapez: " <<  res1 << endl;
  cout << "Mittel: " <<  res2 << endl;
  cout << "Simpson: " <<  res3 << endl;
  cout << "Ende des Programms!" << endl;
  return 0;
}
