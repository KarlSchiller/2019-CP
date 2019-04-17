#include <iostream>
#include <Eigen/Dense>

// using Eigen::Matrix3i;
// using Eigen::PartialPivLU;
using namespace std;
using namespace Eigen;

int main() {
  // Teil a)
  // stelle Matrix auf
  MatrixXd A(10, 2);
  A << 0, 1, 2.5, 1, -6.3, 1, 4, 1, -3.2, 1, 5.3, 1, 10.1, 1, 9.5, 1, -5.4, 1, 12.7, 1;
  // Vektor y
  VectorXd y(10);
  y << 4, 4.3, -3.9, 6.5, 0.7, 8.6, 13, 9.9, -3.6, 15.1;

  // Teil b)
  MatrixXd P = A.transpose()*A;
  cout << "Matrix P = A.T*A" << endl << P << endl;
  VectorXd b = A.transpose()*y;
  cout << "Vektor b' = A.T*b" << endl << b << endl;

  // Teil c)
  // LU-Zerlegung unter Verwendung der Bibliothek Eigen zur Berechnung der Inversen von P
  MatrixXd copy = P;
  // PartialPivLU<Ref<MatrixXd>> luc(LU);
  PartialPivLU<Ref<MatrixXd>> LU(copy);
  MatrixXd L(2, 2);
  L = MatrixXd::Identity(2,2);
  L.block<2,2>(0,0).triangularView<Eigen::StrictlyLower>()= LU.matrixLU();
  cout << "L" << endl << L << endl;
  MatrixXd U(2, 2);
  U = LU.matrixLU().triangularView<Upper>();
  cout << "U" << endl << U << endl;
  MatrixXd Perm(2, 2);
  Perm = LU.permutationP();
  cout << "P" << endl << Perm << endl;
  // Berechne Inverse von P
  MatrixXd I(2, 2);
  // I = LU.inverse(); Der einfache Weg, wäre der erlaubt?...
  I = U.inverse()*L.inverse()*Perm.transpose();
  cout << "P^(-1)" << endl << I << endl;

  // Berechne Schätzer für Parameter m und n
  VectorXd a(2);
  a = I*b;
  cout << "Regressierte m und n" << endl << a << endl;
}
