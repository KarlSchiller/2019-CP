#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>

using namespace std;
using namespace Eigen;

//MatrixXd init_federmatrix(VectorXd m, VectorXd k, const int N)
//{
//  // build diagonal entries of Matrix as vector diag
//  VectorXd diag = m;
//  diag.head(N-1) += k;
//  diag.tail(N-1) += k;
//  MatrixXd A = diag.asDiagonal();
//
//  A.diagonal<1>() = -k;
//  A.diagonal<-1>() = -k;
//  return A;
//}
//
//VectorXd federkette(VectorXd m, VectorXd k, const int N)
//{
//  MatrixXd A = init_federmatrix(m, k, N);
//  // compute eigenvalues
//  VectorXd eivals = A.eigenvalues().real();
//  std::sort(eivals.data(), eivals.data()+eivals.size());
//  // a function to compute the piece-wise sqare root is missing :(
//  for (int i=0; i<N; i++){
//    eivals(i) = sqrt(eivals(i));
//  }
//  return eivals;
//}
//
//int main()
//{
//  cout << "Aufgabe 1" << endl;
//
//  const int N = 10;
//  VectorXd m = VectorXd::LinSpaced(N,1,N);
//  VectorXd k = N*VectorXd::Ones(N-1)-VectorXd::LinSpaced(N-1,1,N-1);
//  VectorXd eivals = federkette(m, k, N);
//  cout << eivals << endl;
//
//  return 0;
//}
