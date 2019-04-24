#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "Dateien/service.cpp"

using namespace std;
using namespace Eigen;

// Approximation der Singulärwerte
VectorXd sort_vec(VectorXd &vector, int k){
    VectorXd approx();
        //for (int i = 0; i < k; i++){
        //    approx(i) = vector(i);
        //}
    return vector;
}

// Approximation der Matrizen
MatrixXd sort_mat(MatrixXd &M, int k) {
    MatrixXd approx(512, 512);
    approx = M;
    //for (int i = 0; i < k; i++){
    //    for (int j = 0; j < k; j++){
    //        approx(i, j) = M(i, j);
    //    }
    //}
    int zeros = 512 - k;
    approx.rightCols(zeros).setZero();
    return approx;
}

int main()
{
    cout << "Beginn des Programms!" << endl;
    // Einlesen der Daten
    MatrixXd M;
    int k[3] = {10, 20, 50};
    loadData(M, "Dateien/Bild", 512, 512);
    cout << M.rows() << "x" << M.cols() << endl;
    // SVD durchführen
    BDCSVD<MatrixXd> svd(M, ComputeFullU|ComputeFullV);
    VectorXd sing = svd.singularValues();
    MatrixXd U = svd.matrixU();
    MatrixXd V = svd.matrixV();

    // Initialisieren der Variablen
    VectorXd approx;
    MatrixXd approxW, approxU, approxV, approxA;
    ofstream file;
    string filename;

    // Durchführen für die verschiedenen k-Werte
    for (int l=0; l<3; l++){
      // Approximation der Singulärwerte und umschreiben in eine Diagonalmatrix
      approx = sort_vec(sing, k[l]);
      approxW = approx.asDiagonal();

      // Approximation der U- und der V-Matrix
      approxU = sort_mat(U, k[l]);
      approxV = sort_mat(V, k[l]);

      // Transformieren der Matrix
      approxA = approxU*approxW*approxV.transpose();

      // Auslesen in eine txt-Datei
      filename = "build/bild_"+to_string(k[l])+".txt";
      file.open(filename, ios::trunc);
      //file << "# Array" << endl;
      for (int i = 0; i < 512; i++){
          file << i << ";";
      }
      file << endl;
      for (int i=0; i < 512; i++){
          for (int j = 0; j < 512; j++){
              file << approxA(i, j) << "; ";
          }
          file << endl;
      }
      file.close();
    }
    cout << "Ende des Programms!" << endl;
    return 0;
}
