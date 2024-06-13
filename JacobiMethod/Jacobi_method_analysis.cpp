#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
using namespace std;
#define EPS 1e-5
const double PI = acos(-1);
const int MAX_STEP = 5000;
int n=4;
vector<vector<double>>A, P, orMtx;
vector<double>Bp, Bq, errVec, offDiag;
// Compute sum of the abs of upper off-diagonal entries
void count_offDiag(int k) {
	double ret = 0;
	for (int i = 0;i < n;i++)
		for (int j = i + 1;j < n;j++)
			ret += fabs(A[i][j]);
	offDiag[k] = ret;
}

void make_matrix() {
	vector<double>tp(n, 0);
	A.clear(), P.clear();
	A.resize(n, tp);
	P.resize(n, tp);
	for (int i = 0;i < n;i++)
		for (int j = i;j < n;j++)
			A[i][j] = A[j][i] = rand() % 20;
}
void update_P_mtx(vector<vector<double>>& P, int p, int q, double c, double s, int n) {
	for (int k = 0;k < n;k++) {
		Bp[k] = c * P[k][p] + s * P[k][q];
		Bq[k] = (-s) * P[k][p] + c * P[k][q];
	}
	for (int k = 0;k < n;k++) {
		P[k][p] = Bp[k];
		P[k][q] = Bq[k];
	}
}
void update_A_mtx(vector<vector<double>>& A, int p, int q, double c, double s, int n) {
	for (int k = 0;k < n;k++) {
		if (k != p && k != q) Bp[k] = c * A[k][p] + s * A[k][q];
		if (k != p && k != q) Bq[k] = -s * A[k][p] + c * A[k][q];
	}
	Bp[p] = c * c * A[p][p] + 2.0 * s * c * A[p][q] + s * s * A[q][q];
	Bq[q] = s * s * A[p][p] - 2.0 * s * c * A[p][q] + c * c * A[q][q];
	Bp[q] = Bq[p] = (c * c - s * s) * A[p][q] + s * c * (A[q][q] - A[p][p]);
	for (int k = 0;k < n;k++) {
		A[p][k] = A[k][p] = Bp[k];
		A[q][k] = A[k][q] = Bq[k];
	}
}
void make_identity_mtx(vector<vector<double>>& A, int N) {
	for (int i = 0;i < N;i++) {
		for (int j = 0;j < N;j++) {
			if (i != j) A[i][j] = 0.0;
			else A[i][j] = 1.0;
		}
	}
}
double max_off_diag_entry(vector<vector<double>>& A, int* p, int* q, int N) {
	double max = A[0][1];
	*p = 0;
	*q = 1;
	for (int i = 0;i < N;i++) {
		for (int j = i + 1;j < N;j++)
			if (fabs(A[i][j]) > fabs(max)) {
				*p = i;
				*q = j;
				max = A[i][j];
			}
	}
	return(max);
}
int jacobian_method(vector<vector<double>>& A, vector<vector<double>>& P, int n) {
	int p, q;
	make_identity_mtx(P, n);
	Bp.clear(), Bq.clear();
	Bp.resize(n, 0);
	Bq.resize(n, 0);
	double A_pq = max_off_diag_entry(A, &p, &q, n);
	int k = 0;
	while (fabs(A_pq) > EPS && k < MAX_STEP) {
		double a, b, theta;
		a = 2.0 * A_pq;
		b = A[p][p] - A[q][q];
		if (b == 0.0) {
			if (a > 0.0) theta = PI / 2.0;
			else theta = 3.0 * PI / 2.0;
		}
		else
			theta = atan(a / b);
		theta = theta / 2.0;
		double c = cos(theta), s = sin(theta);
		update_P_mtx(P, p, q, c, s, n);
		update_A_mtx(A, p, q, c, s, n);
		A_pq = max_off_diag_entry(A, &p, &q, n);
		k++;
		count_offDiag(k);
	}
	return (k);
}
double  inner_product(vector<double>& a, vector<double>& b, int n) {
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += a[i] * b[i];
	return (sum);
}
void  retrieve_eigen_values(vector<double>& v) {
	for (int i = 0; i < n; i++) v[i] = A[i][i];
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			double temp = P[i][j];
			P[i][j] = P[j][i];
			P[j][i] = temp;
		}
	}
	vector<double>tp(n, 0);
	orMtx.resize(n, tp);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			orMtx[i][j] = inner_product(P[i], P[j], n);
}
double vec_norm(vector<double>& a, int n) {
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += a[i] * a[i];
	sum = sqrt(sum);
	return(sum);
}
void normalize_vec(vector<double>& x, int n) {
	double norm = vec_norm(x, n);
	if (norm == 0) return;
	for (int i = 0; i < n; i++)
		x[i] = x[i] / norm;
}
void mtx_vec_mult(vector<double>& a, vector<vector<double>> A, vector<double>& b, int n) {
	for (int i = 0; i < n; i++) {
		a[i] = 0.0;
		for (int j = 0; j < n; j++)
			a[i] += A[i][j] * b[j];
	}
}
double norm(vector<double>& v) {// 2 norm
	double ret = 0;
	for (int i = 0; i < n; i++)ret += v[i] * v[i];
	return sqrt(ret);
}
void comp_err_vec(vector<double>& v, int n) {
	vector<double>a(n, 0), b(n, 0), t(n, 0);
	errVec.clear();
	errVec.resize(n, 0);
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < n; k++) t[k] = P[i][k];
		mtx_vec_mult(a, A, t, n);
		for (int k = 0; k < n; k++) b[k] = a[k] - v[i] * t[k];
		errVec[i] = vec_norm(b, n);
	}
}
int main() {
	//freopen("D:/user/Desktop/學期/數值分析/hw6/data_delta10.in", "r", stdin); 
	//freopen("D:/user/Desktop/學期/數值分析/hw6/data_delta10.out", "w", stdout); 
	srand(time(NULL));
	//for (n = 4;n <= 40;n++) {
	n = 4;
		make_matrix();
		offDiag.clear();
		offDiag.resize(MAX_STEP + 1, 0);
		count_offDiag(0);
		int iter = jacobian_method(A, P, n);
		cout << iter << " \n";
		//for (int i = 0;i <= iter;i++)cout << offDiag[i]/ offDiag[0] << '\n';
	//}
	return 0;
}