#include <iostream>
#include <fstream>
#include <vector>
#include <eigen/Dense>
#include "Dateien/service.cpp"

using namespace std;
using namespace Eigen;


MatrixXd reduce_mat(MatrixXd &M, int k) {
    MatrixXd approx(k, k);
    for (int i = 0; i < k; i++){
        for (int j = 0; j < k; j++){
            approx(i, j) = M(i, j);
        }
    }
    return approx;
}


int main()
{
    cout << "Aufgabe 2" << endl;
    // Teil a)
    // Einlesen
    MatrixXd A;
    cout << "\tImportiere Trainingsdaten" << endl;
    int test = loadData(A, "Dateien/Training", 112*92, 360);
    if(test==0){
      cout << "\tFehler beim Einlesen der Daten!" << endl;
      return 1;
    }
    // cout << A.rows() << "x" << A.cols() << endl;

    // SVD durchführen
    cout << "\tFuehre SVD durch" << endl;
    BDCSVD<MatrixXd> svd(A, ComputeFullV|EigenvaluesOnly);
    VectorXd sing = svd.singularValues();
    MatrixXd V = svd.matrixV();

    // Speichere Sprektrum zum plotten im build-Ordner
    cout << "\tSpeichern der Eigenwerte" << endl;
    ofstream evfile;
    evfile.open("build/aufg3-eigenvalues.txt", ios::trunc);
    evfile << "evalues" << endl;
    for (int i=0; i<360; i++){
      evfile << sing(i) << endl;
    }
    evfile.close();

    // Nehme die ersten k Eigenwerte
    int k = 200;
    MatrixXd reducedV = reduce_mat(V, k);
    // TODO: Transformiere erstes Bild in A und plotte es als heatmap


    // Teil b)

    // BDCSVD<MatrixXd> svd(M, ComputeFullU|ComputeFullV);
    // VectorXd sing = svd.singularValues();
    // VectorXd approx = sort_vec(sing, k);
    // MatrixXd approxW = approx.asDiagonal();

    // MatrixXd U = svd.matrixU();
    // MatrixXd V = svd.matrixV();
    // MatrixXd approxU = sort_mat(U, k);
    // MatrixXd approxV = sort_mat(V, k);

    // MatrixXd approxA = approxU*approxW*approxV.transpose();
    return 0;
}
