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
#define EPS 1e-6
const double PI = acos(-1);
const int MAX_STEP = 50;
int n = 4;
vector<vector<double>>A,P,orMtx;
vector<double>Bp, Bq, errVec;
void print_matrix(vector<vector<double>>&M,string s) {
	cout << s << ":\n";
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++)
			cout << setw(10) << fixed << setprecision(6) << M[i][j] << ' ';
		cout << '\n';
	}
	cout << '\n';
}
void diagonal_dominant() {
	for (int i = 0; i < n; i++) {
		A[i][i] = 5.0;
		for (int j = 0; j < n; j++) {
			if (i != j)
				A[i][i] += fabs(A[i][j]);
		}
	}
}
void make_matrix() {
	vector<double>tp(n, 0);
	A.resize(n, tp);
	P.resize(n, tp);
	for (int i = 0;i < n;i++)for (int j = i;j < n;j++)A[i][j] = A[j][i] = rand() % 20;
	diagonal_dominant();
}
void update_P_mtx(vector<vector<double>>& P, int p, int q, double c, double s, int n){
	//Compute the new columns
	for (int k = 0;k < n;k++) {
		Bp[k] = c * P[k][p] + s * P[k][q];
		Bq[k] = (-s) * P[k][p] + c * P[k][q];
	}
	//Copy into the P[][] mtx.
	for (int k = 0;k < n;k++) {
		P[k][p] = Bp[k];
		P[k][q] = Bq[k];
	}
}
void update_A_mtx(vector<vector<double>>& A, int p, int q, double c, double s, int n){
	//Compute new values of  the p- and q-col (row)
	for (int k = 0;k < n;k++) {
		if (k != p && k != q) Bp[k] = c * A[k][p] + s * A[k][q];
		if (k != p && k != q) Bq[k] = -s * A[k][p] + c * A[k][q];
	}
	//Compute the new A[p][p] & A[q][q].
	Bp[p] = c * c * A[p][p] + 2.0 * s * c * A[p][q] + s * s * A[q][q];
	Bq[q] = s * s * A[p][p] - 2.0 * s * c * A[p][q] + c * c * A[q][q];
	//Update A[p][q] amd A[q][p]
	Bp[q] = Bq[p] = (c * c - s * s) * A[p][q] + s * c * (A[q][q] - A[p][p]);
	//enter the new col. and rows.
	for (int k = 0;k < n;k++) {
		A[p][k] = A[k][p] = Bp[k];
		A[q][k] = A[k][q] = Bq[k];
	}
	// A[p][q] = A[q][p] = 0.0;
}
void make_identity_mtx(vector<vector<double>>& A, int N){
	for (int i = 0;i < N;i++) {
		for (int j = 0;j < N;j++) {
			if (i != j) A[i][j] = 0.0;
			else A[i][j] = 1.0;
		}
	}
}
double max_off_diag_entry(vector<vector<double>>& A, int* p, int* q, int N){
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
int jacobian_method(vector<vector<double>>& A, vector<vector<double>>& P, int n){
	int p, q;
	// Make P = I.
	make_identity_mtx(P, n);
	// Create the temporary space hodoubleing the p & q columns
	Bp.resize(n, 0);
	Bq.resize(n, 0);

	// Find the max off-diagonal entry.
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
		//Update the p & q col. & rows of P[][].
		update_P_mtx(P, p, q, c, s, n);
		// A = R'*A*R.
		update_A_mtx(A, p, q, c, s, n);
		// Find the max off-diagonal entry for next iteration.
		A_pq = max_off_diag_entry(A, &p, &q, n);
		//Print matrix A[][] for error checking.
		k++;
		cout << "k = " << k << ", theta = " << theta << '\n';
		print_matrix(A,"A");
	}
	return (k);
}
double  inner_product(vector<double>& a, vector<double>& b, int n){
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += a[i] * b[i];
	return (sum);
}
// Retrieve eigen value and test the orthogonality of the eigen-vectors
void  retrieve_eigen_values(vector<double>& v){
	// Retrieve eigen value
	for (int i = 0; i < n; i++) v[i] = A[i][i];
	// Transpose the eigen-vector matrix.
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			double temp = P[i][j];
			P[i][j] = P[j][i];
			P[j][i] = temp;
		}
	}
	// Test the orthogonality of the eigen-vectors.
	vector<double>tp(n, 0);
	orMtx.resize(n, tp);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			orMtx[i][j] = inner_product(P[i], P[j], n);
	cout << "\nOrthogonality between eigen-vectors:\n";
	print_matrix(orMtx, "orMtx");
}
double vec_norm(vector<double>& a, int n){
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += a[i] * a[i];
	sum = sqrt(sum);
	return(sum);
}
void normalize_vec(vector<double>& x, int n){
	double norm = vec_norm(x, n);
	if (norm == 0) return;
	for (int i = 0; i < n; i++)
		x[i] = x[i] / norm;
}
void mtx_vec_mult(vector<double>& a,vector<vector<double>> A, vector<double>& b, int n){
	for (int i = 0; i < n; i++) {
		a[i] = 0.0;
		for (int j = 0; j < n; j++)
			a[i] += A[i][j] * b[j];
	}
}
/*
	Procedure to compute the norms of || A * v - u * v ||
	where A is the matrix, v the eigen vector,and u the eigen value.
	Compute this norm for all eigen values and their eigen vectors.
	Keep them in errVec[].
*/
double norm(vector<double>& v) {// 2 norm
	double ret = 0;
	for (int i = 0; i < n; i++)ret += v[i] * v[i];
	return sqrt(ret);
}
void comp_err_vec(vector<double>& v, int n){
	vector<double>a(n, 0), b(n, 0), t(n, 0);
	errVec.resize(n, 0);
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < n; k++) t[k] = P[i][k];
		// a = A*v.
		mtx_vec_mult(a, A, t, n);
		// b = a - u*v
		for (int k = 0; k < n; k++) b[k] = a[k] - v[i] * t[k];
		// err = ||b||
		errVec[i] = vec_norm(b, n);
	}
	cout << "Residual of ||A*v - u*v|| = ";
	for (int i = 0; i < n; i++)
		cout << errVec[i] << ' ';
	cout << "\n";

	// Compute 2-norm of ||A*v - u*v||
	cout << "2-norm of ||A*v - u*v|| = " << norm(errVec) << "\n\n";
}
int main() {
	srand(time(NULL));
	make_matrix();
	print_matrix(A,"A");
	int iter = jacobian_method(A, P, n);
	cout << "Need " << iter << " times iterations\n";

	// Retrieve eigen value and test the orthogonality of the eigen-vectors
	vector<double> eigenVec(n, 0);
	retrieve_eigen_values(eigenVec);
	cout << "EigenValues = ";
	for (int i = 0; i < n; i++)
		cout << eigenVec[i] << ' ';
	cout << '\n';

	// Find maximum and minimun eigen value(abs)
	double max = 0, min = 0;
	for (int i = 0; i < n; i++) {
		if (fabs(eigenVec[i]) > fabs(eigenVec[max])) max = i;
		if (fabs(eigenVec[i]) < fabs(eigenVec[min])) min = i;
	}
	// Find eigne vectors
	cout << "The eigen-vectors are: \n";
	for (int i = 0; i < n; i++) {
		cout << "x[" << i << "] = ";
		for (int j = 0; j < n; j++)
			cout << setw(10) << P[i][j] << ' ';
		cout << "\n";
	}
	cout << '\n';

	// Compute the condition number
	comp_err_vec(eigenVec, n);
	double condNum = fabs(eigenVec[max] / eigenVec[min]);
	cout << "max eigenvalue = " << eigenVec[max] << '\n';
	cout << "min eigenvalue = " << eigenVec[min] << '\n';
	cout << "Condition number = " << condNum << "\n\n";
	return 0;
}