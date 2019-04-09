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
    //for (int i = 0; i < N; i++){
    //    for (int j = 0; j < N; j++){
    //        cout << U[i][j] << " ";
    //    }
    //    cout << "\n";
    //}
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
                int max = abs(Ub[i+1][i]);
                int max_i = i+1;
                for (int m = i+1; m < N; m++){
                    for(int k = 0; k <N; k++){
                        if(abs(Ub[m+1][k]) > max){
                            max_i = m;
                            max = abs(Ub[m+1][k]);
                        }
                    }
                }
               // cout << max << " ";
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
        // cout << endl;



    // for (int i = 0; i < N; i++){
    //     for (int j = 0; j < N; j++){
    //         cout << Lb[i][j] << " ";
    //     }
    //     cout << "\n";
    // }
    // cout << endl;

    // Aufgabenteil c)
    // Matrix A initialisieren
    MatrixXd Ac(3,3);
    Ac << 1, 2, 4, -2, 1, 0, 4, 2, -4;
    // cout << Ac << endl;
    MatrixXd M = Ac;
    PartialPivLU<Ref<MatrixXd> > lu(M);
    cout << M << endl;
    cout << endl;

    // Definition von U und L, gewonnen aus M
    MatrixXd Uc(3,3);
    Uc << 4, 2, -4, 0, 2, -2, 0, 0, 6.5;
    MatrixXd Lc(3,3);
    Lc << 1, 0, 0, -0.5, 1, 0, 0.25, 0.75, 1;
    // Lc*Uc ergibt das Pivotisierte A
    cout << Lc*Uc << endl;
    cout << endl;
    // Daher kenne ich die Pivotisierungsmatrix :)
    MatrixXd Pc(3, 3);
    Pc << 0, 0, 1, 0, 1, 0, 1, 0, 0;
    PartialPivLU<Ref<MatrixXd> > luc(Pc);
    MatrixXd Pc1 = luc.solve(Lc*Uc);
    cout << Pc1 << endl;
    cout << endl;
    // Also wenn man die Pivotisierungsmatrix kennt, kann man mit solve wieder A herausbekommen
    MatrixXd At = Ac.transpose();
    MatrixXd Ac1 = Lc*Uc;
    MatrixXd Ac2 = Ac1.transpose();
    luc.compute(At);
    MatrixXd Loesung = luc.solve(Ac2);
    cout << Loesung << endl;




    cout << "Ende des Programms" <<endl;
    return 0;
}
