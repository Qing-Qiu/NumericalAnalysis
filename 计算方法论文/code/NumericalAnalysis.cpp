#include<bits/stdc++.h>

using namespace std;

const int N = 405;
const double eps = 1e-6;
int n = 19;
double mat[N][N << 1];
double A[N] = {0, 306.35, 311.2, 324.9, 343.4, 361.61, 390.76, 421.36,
               458.75, 494.9, 522.78, 547.22, 588.7, 647.18,
               712.34, 778.32, 835.31, 888.1, 923.16, 939.48,
               988.60, 1073.62, 1164.29},
        B[2 + 1][N], B_T[N][2 + 1], C[N], c[2 + 1][1 + 1];
double Y[N][1 + 1], P[N], F[N], G[N], inv[3][3];
double tmp[3][N];
double a, b;

void get_inv(int n) { //高斯消元
    for (int i = 1, r; i <= n; i++) {
        r = i;
        for (int j = i + 1; j <= n; j++)
            if (mat[j][i] > mat[r][i]) r = j;
        if (r != i) swap(mat[i], mat[r]);
        if (abs(mat[i][i]) < eps) {
            puts("No Solution");
            return;
        }
        double tmp = mat[i][i];
        for (int k = 1; k <= n; k++) {
            if (k == i) continue;
            double p = mat[k][i] / tmp;
            for (int j = i; j <= (n << 1); j++)
                mat[k][j] = mat[k][j] - p * mat[i][j];
        }
        for (int j = 1; j <= (n << 1); j++)
            mat[i][j] = mat[i][j] / tmp;
    }
    for (int i = 1; i <= n; i++)
        for (int j = n + 1; j <= (n << 1); j++)
            inv[i][j - n] = mat[i][j];
}

int main() {
    P[1] = A[1];
    for (int i = 2; i <= n; i++)
        P[i] = P[i - 1] + A[i];
    for (int i = 1; i < n; i++)
        C[i] = (P[i] + P[i + 1]) / 2;
    for (int i = 1; i < n; i++) {
        B[1][i] = B_T[i][1] = -C[i];
        B[2][i] = B_T[i][2] = 1;
    }
    for (int i = 1; i < n; i++)
        Y[i][1] = A[i + 1];
    for (int i = 1; i <= 2; i++)
        for (int j = 1; j <= 2; j++)
            for (int k = 1; k < n; k++)
                mat[i][j] += B[i][k] * B_T[k][j];
    for (int i = 1; i <= 2; i++)
        for (int j = 1; j <= 2; j++) {
            mat[i][i + 2] = 1;
        }
    get_inv(2);
    for (int i = 1; i <= 2; i++)
        for (int j = 1; j < n; j++)
            for (int k = 1; k <= 2; k++)
                tmp[i][j] += inv[i][k] * B[k][j];
    for (int i = 1; i <= 2; i++)
        for (int k = 1; k < n; k++)
            c[i][1] += tmp[i][k] * Y[k][1];
    a = c[1][1], b = c[2][1];
    F[1] = A[1];
    for (int i = 2; i <= n + 3; i++)
        F[i] = (A[1] - b / a) / exp(a * (i - 1)) + b / a;
    G[1] = A[1];
    for (int i = 2; i <= n + 3; i++)
        G[i] = F[i] - F[i - 1];
    for (int i = 1; i <= n + 3; i++)
        cout << 1979 + i << ' ' << G[i] << ' ' << (G[i] - A[i]) / A[i] * 100 << endl;
    return 0;
}