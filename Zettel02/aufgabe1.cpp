#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Vorimplementierte Methoden
int loadData(Eigen::MatrixXd &mat, std::string filename, size_t row, size_t col)
{

    std::ifstream input;
    input.open(filename, std::fstream::in | std::fstream::binary);
    if (input.is_open())
    {
        std::vector<unsigned char> vec((std::istreambuf_iterator<char>(input)),
                std::istreambuf_iterator<char>() );
        mat= Eigen::Map<Eigen::Matrix<unsigned char, Eigen::Dynamic,Eigen::Dynamic>>(vec.data(),row,col).template cast<double>();
        input.close();
        return 1;
    }
    return 0;

}

int storeData(Eigen::MatrixXd &mat, std::string filename)
{
    Eigen::Matrix<char, Eigen::Dynamic, Eigen::Dynamic> pic= mat.template cast<char>();
    std::ofstream output;
    output.open(filename,std::fstream::trunc);
    if (output.is_open())
    {
        output.write(pic.data(),pic.size());
        output.close();
        return 1;
    }
    return 0;


}

// Approximation der Singulärwerte
VectorXd sort_vec(VectorXd &vector, int k){
    VectorXd approx(k);
        for (int i = 0; i < k; i++){
            approx(i) = vector(i);
        }
    return approx;
}

// Approximation der Matrizen
MatrixXd sort_mat(MatrixXd &M, int k) {
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
    cout << "Beginn des Programms!" << endl;
    // Einlesen der Daten
    MatrixXd M;
    const int k = 50;
    loadData(M, "../Dateien/Bild", 512, 512);
    cout << M.rows() << "x" << M.cols() << endl;
    // SVD durchführen
    BDCSVD<MatrixXd> svd(M, ComputeFullU|ComputeFullV);
    VectorXd sing = svd.singularValues();

    // Approximation der Singulärwerte und umschreiben in eine Diagonalmatrix
    VectorXd approx = sort_vec(sing, k);
    MatrixXd approxW = approx.asDiagonal();

    // Approximation der U- und der V-Matrix
    MatrixXd U = svd.matrixU();
    MatrixXd V = svd.matrixV();
    MatrixXd approxU = sort_mat(U, k);
    MatrixXd approxV = sort_mat(V, k);

    // Transformieren der Matrix
    MatrixXd approxA = approxU*approxW*approxV.transpose();

    // Auslesen in eine txt-Datei
    ofstream file;
    file.open("bild_50.txt", ios::trunc);
    //file << "# Array" << endl;
    for (int i = 0; i < k; i++){
        file << i << ";";
    }
    file << endl;
    for (int i=0; i < k; i++){
        for (int j = 0; j < k; j++){
            file << approxA(i, j) << "; ";
        }
        file << endl;
    }
    file.close();
    cout << "Ende des Programms!" << endl;
    return 0;
}
