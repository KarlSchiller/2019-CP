#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>
#include <tuple>

using namespace std;
using namespace Eigen;

tuple<VectorXd, double> lanczos(MatrixXd A, int N, double genauigkeit)
{
  // Startvektoren
  VectorXd q0 = VectorXd::Zero(N);
  VectorXd qi = VectorXd::Zero(N);
  qi(0) = 1;
  // Array mit Einträgen der Tridiagonalmatrix
  VectorXd alpha = VectorXd::Zero(N);
  VectorXd beta = VectorXd::Zero(N+1);
  beta(0) = 1;
  // Benoetigt fuer die Iterationen
  VectorXd gamma = VectorXd::Ones(N);
  VectorXd gamma1 = VectorXd::Zero(N);
  int i = 0;

  // Testzwecke
  // MatrixXd Q = MatrixXd::Zero(N,N+1);
  // Q.col(0) = qi;

  while(gamma.norm() > genauigkeit && i < N)
  {
    alpha(i) = qi.transpose() * A * qi;
    gamma = (A - alpha(i)*MatrixXd::Identity(N,N))*qi - beta(i)*q0;
    if(i == 0){
      gamma1 = gamma;
    }
    beta(i+1) = sqrt(gamma1.dot(gamma));
    // update fuer naechsten Durchlauf
    q0 = qi;
    qi = gamma / beta(i+1);
    // Q.col(i+1) = qi;  // test
    i++;
  }

  // Test, ob Q.T * T * Q wieder A ergibt
  // MatrixXd test = Q.block(0,0,N,N);
  // MatrixXd T = MatrixXd::Zero(N, N);
  // T.diagonal() = alpha;
  // T.diagonal<1>() = beta.block(1,0,N-1,1);
  // T.diagonal<-1>() = beta.block(1,0,N-1,1);
  // cout << T << endl << endl;
  // cout << test.transpose()*T*test << endl;  // sollte wieder A ergeben

  VectorXd evec;
  double eval;
  return make_tuple(evec, eval);
}

int main()
{
  cout << "Aufgabe 1" << endl;

  const int N = 10;
  MatrixXd test = MatrixXd::Ones(N,N);
  test(0, 0) = 3;
  test(1, 1) = 9;
  // Grundzustand und zugehörige Energie
  VectorXd evec;
  double eval;
  tie(evec, eval) = lanczos(test, N, 0.001);


  // old
  // const int N = 10;
  // VectorXd m = VectorXd::LinSpaced(N,1,N);
  // VectorXd k = N*VectorXd::Ones(N-1)-VectorXd::LinSpaced(N-1,1,N-1);
  // VectorXd eivals = federkette(m, k, N);
  // cout << eivals << endl;

  return 0;
}
