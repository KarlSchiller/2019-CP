#include <iostream>
#include <../Eigen/Dense>

// using Eigen::Matrix3i;
// using Eigen::PartialPivLU;
using namespace std;
using namespace Eigen;

int main()
{
    cout << "Beginn des Programms" << endl;
    int N = 4;
    int A[4][4] = {
        {1, 5, -4, 2},
        {-4, -17, 14, -5},
        {2, 16, -10, 7},
        {6, 33, -22, 13}
    };
    int L[N][N]{};
    int U[N][N]{};
    for (int i = 0; i < N; i++){
        L[i][i] = 1;
    }
    for (int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            U[i][j] = A[i][j];
        }
    }
// intitialisieren beendet
    int rij = 0;
    for (int i = 0; i < N; i++){
        for (int j = i + 1; j < N; j++){
            rij = U[j][i]/U[i][i];
            for (int k = 0; k < N; k++){
                if (A[i][i] != 0){
                   // cout << j  << " " << rij << " ";
                   // cout << endl;
                   // cout << U[j][k]-rij*U[i][k] << " ";
                    U[j][k] = (U[j][k])-(U[i][k]*rij);
                    L[j][i] = rij; // L ist richtig
                }
                else{
                    cout << "Pivotisierung notwendig!";
                    break;
                }
            }
        }
    }
    cout << endl;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            cout << U[i][j] << " ";
        }
        cout << "\n";
    }
    cout << endl;

    // Aufgabenteil b)
    // initialisieren
    int Ab[N][N] = {
        {1, 5, -4, 2},
        {1, 5, -22, 13},
        {-4, 17, 14, 5},
        {2, 16, -10, 7}
    };
    int P[N][N]{};
    for (int i = 0; i < N; i++){
        P[i][i] = 1;
    }
    int Lb[N][N]{};
    int Ub[N][N]{};
    for (int i = 0; i < N; i++){
        Lb[i][i] = 1;
    }
    for (int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            Ub[i][j] = Ab[i][j];
        }
    }
    // Initialisierung beendet
    for (int i = 0; i < N; i++){
        for (int j = i + 1; j < N; j++){
            // Betragsgrößtes element finden
                int max = abs(Ub[i+1][0]);
                int max_i = i+1;
                for (int m = i+1; m < N; m++){
                    for(int k = 0; k <N; k++){
                        if(abs(Ub[m+1][k]) > max){
                            max_i = m;
                            max = abs(Ub[m+1][k]);
                        }
                    }
                }
                cout << max << " ";
                // cout << " \n";
                // zeilen tauschen
                for (int l = 0; l < N; l++){
                    int c = Ub[i][l];
                    Ub[i][l] = Ub[max_i][l];
                    Ub[max_i][l] = c;
                }
                for (int l = 0; l < N; l++){
                    int c = P[i][l];
                    P[i][l] = P[max_i][l];
                    P[max_i][l] = c;
                }
                for (int l = 0; l < i; l++){
                    int c = Lb[i][l];
                    Lb[i][l] = Lb[max_i][l];
                    Lb[max_i][l] = c;
                }
                rij = Ub[j][i]/Ub[i][i];
            for (int k = 0; k < N; k++){
                   // cout << Ub[j][k]-rij*Ub[i][k] << " ";
                    Ub[j][k] = (Ub[j][k])-(Ub[i][k]*rij);
                    Lb[j][i] = rij; // L ist richtig
                }
            }
           // cout << "\n";
        }
    cout << endl;



    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            cout << Ub[i][j] << " ";
        }
        cout << "\n";
    }
    cout << endl;

    // Aufgabenteil c)
    // Matrix A initialisieren
    Matrix3i Ac;
    Ac(0, 0) = 1;
    Ac(0, 1) = 2;
    Ac(0, 2) = 4;
    Ac(1, 0) = -2;
    Ac(1, 1) = 1;
    Ac(1, 2) = 0;
    Ac(2, 0) = 4;
    Ac(2, 1) = 2;
    Ac(2, 2) = -4;
    cout << Ac << endl;
    Matrix3i Lol = Ac.PartialPivLU<Matrix3i>().matrixLU();


    cout << "Ende des Programms" <<endl;
    return 0;
}
