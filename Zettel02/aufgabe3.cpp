#include <iostream>
#include <eigen/Dense>
#include <eigen/Eigenvalues>
#include <tuple>

using namespace std;
using namespace Eigen;

// Maximales Nebendiagonalelement finden
tuple<int, int, double> finde_neben(MatrixXd &M){
    // Dimension herausfinden
    int N = M.cols();
    int zeile, spalte;
    double max = 0;
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i != j){
                if(abs(M(i, j)) > max){
                    max = abs(M(i, j));
                    zeile = i;
                    spalte = j;
                }
            }
        }
    }
    return make_tuple(zeile, spalte, max);
}

// c und s berechnen
double berechne_c(double &t){
    return 1/sqrt(1+t*t);
}

double berechne_s(double &t){
    return berechne_c(t)*t;
}

double berechne_t(double &Omega){
    if(Omega > 0){
        return 1/(Omega + sqrt(Omega*Omega + 1));
    }
    if(Omega < 0){
        return 1/(Omega - sqrt(Omega*Omega + 1));
    }
}

double berechne_off(MatrixXd &A){
    int N = A.cols();
    double summe = 0;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if(i != j){
                summe = summe + abs(A(i, j))*abs(A(i, j));
            }
        }
    }
    return summe;
}

MatrixXd trafo_matrix(MatrixXd &A, int &zeile, int &spalte){
    const int N = A.cols();
    MatrixXd M(N, N);
    double omega = (A(spalte, spalte) - A(zeile, zeile))/(2*A(zeile, spalte));
    cout << endl;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if ((i == zeile && j == zeile) || (i == spalte && j== spalte)){
                double t = berechne_t(omega);
                M(i, j) = berechne_c(t);
            }
            else if (i == zeile && j == spalte){
                double t = berechne_t(omega);
                M(i, j) = berechne_s(t);
            }
            else if (i == spalte && j == zeile){
                double t = berechne_t(omega);
                M(i, j) = (-1)*berechne_s(t);
            }
            else if (i == j){
                M(i, j) = 1;
            }
            else{
                M(i, j) = 0;
            }
        }
    }
    return M;
}

int main()
{
    cout << "Beginn des Programms!" << endl;
    // Aufgabenteil a)

    // Initialisieren der Matrix
    MatrixXd A(4,4);
    A << 1, -2, 2, 4, -2, 3, -1, 0, 2, -1, 6, 3, 4, 0, 3, 5;

    // Berechnen der Eigenwerte mithilfe von Eigen
    //VectorXcd eigenv = A.eigenvalues();
    //cout << "Eigenwerte: " << endl << eigenv << endl;

    // Aufgabenteil b)
    // Größtes Nebendiagonalelement bestimmen
    int zeile_max, spalte_max, lol = 0;
    double max, off;
    MatrixXd A_Trafo, z;
    do {
        cout << "Iteration: " << lol << endl;
        // Maximales Nebendiagonalelement finden
        tie(zeile_max, spalte_max, max) = finde_neben(A);
        //cout << "Max finden ausgeführt!" << endl;
        cout << "Zeile: " << zeile_max << " Spalte: " << spalte_max << endl;

        // Trafo-Matrix berechnen
        z = trafo_matrix(A, zeile_max, spalte_max);
        //cout << "Trafo-Matrix: " << endl << z << endl;
        //cout << "Trafo-Matrix ausgeführt" << endl;

        // Trafo durchführen
        A_Trafo = z.transpose()*A*z;
        //cout << A_Trafo << endl;
        //cout << "Trafo durchgeführt" << endl;

        // Berechne Off
        off = berechne_off(A) - 2*abs(A(zeile_max, spalte_max))*abs(A(zeile_max, spalte_max));
        //cout << "Off: " << off << endl;
        // Ändern
        A = A_Trafo;
        lol++;
        //cout << "A am Ende des Durchgangs: " << endl << A << endl;
    } while(off >= 1e-6);
    cout << "Finale Diagonalmatrix: " << endl << A << endl;
    cout << endl;
    cout << "Diagonalelemente: " << endl << A.diagonal() << endl;


    cout << "Ende des Programms!" << endl;
    return 0;
}
