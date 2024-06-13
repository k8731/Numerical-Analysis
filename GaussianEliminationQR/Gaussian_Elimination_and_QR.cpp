#include<iostream>
#include<cmath>
#include<algorithm>
#include<iomanip>
#include<memory.h>
#include<vector>
#include<string>
#include <cstdio>
#include <cstdlib>
#include <ctime>
using namespace std;
typedef long double ld;
#define EPS0 1e-6
#define EPS1 1e-7
#define EPS2 1e-8
#define EPS3 1e-9
ld epsArray[4] = { EPS0,EPS1,EPS2,EPS3 };
const int N = 21;// # of sample point
const int n = 8;
vector<ld>xs(N, 0), ys(N, 0);// sample points
vector<vector<ld>>B, A, At,pA;// p: perturbation matrix
vector<ld>d(n), c(n), y(N),py;
void print_vector(vector<ld>&v,int n,string s,bool endline) {
	cout << s ;
	if (!endline)cout << '\n';
	for (int i = 0; i < n; i++) {
		if (endline && (i % 7 == 0))cout << '\n';
		cout << setw(12) << fixed << setprecision(16) << v[i] << ' ';
	}
	cout << '\n';
}
void print_matrix(vector<vector<ld>>& v, int n,int m, string s) {
	cout << s << '\n';
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)cout << setw(16) << fixed << setprecision(12) << v[i][j] << ' ';
		cout << '\n';
	}
	cout << '\n';
}
void Horner() {
	cout << "Create sample points:\n";
	// counting x
	for (int i = 0;i < N;i++)xs[i] = 2.0 + i * 0.1;
	print_vector(xs, N, "Xi:",1);
	// counting y
	for (int i = 0;i < N;i++) {
		ld cur = 1;
		for (int j = 6; j >= 0; j--)
			cur = cur * xs[i] + 1;
		ys[i] = cur;
	}
	print_vector(ys, N, "Yi:", 1);
	cout << "\n\n";
}
void back_substitute(vector<vector<ld>>& B,vector<ld>&d){
	cout << "Back substitute:\n";
	int   i, j;
	for (i = n - 1;i >= 0;i--) {
		c[i] = d[i] / B[i][i];
		for (j = i - 1;j >= 0;j--)
			d[j] = d[j] - B[j][i] * c[i];
	}
}
// procedure to add perturbation
void make_perturbation() {
	cout << "Add perturbation:\n";
	// setting perturbation matrix
	vector<ld>ta(n, 0.0);
	pA.resize(N, ta);
	py.resize(N, 0.0);

	// add perturbation
	/*
	srand(time(NULL));
	int rs = rand() % 50;// 0~50個perturbation
	cout << "have " << rs << " perturbations.\n";
	for (int i = 0; i < rs; i++) {
		int r = rand() % N;
		int pos = rand() % 4;
		int fac = rand() % 10;
		py[r] += epsArray[pos] * fac;
	}
	*/
	srand(time(NULL));
	int rs = rand() % 50;// 0~50個perturbation
	cout << "have " << rs << " perturbations.\n";
	for (int i = 0; i < rs; i++) {
		int rx = rand() % N;
		int ry = rand() % n;
		int pos = rand() % 4;
		int fac = rand() % 10;
		pA[rx][ry] += epsArray[pos] * fac;
		//py[r] += epsArray[pos] * fac;
	}

	// setting initial matrix with perturbation
	A = pA;
	for (int i = 0; i < (int)ys.size(); i++)ys[i] += py[i];
	// print data
	print_matrix(pA, N, n, "A_p:");
	print_vector(py, N, "y_p:", 1);
	cout << "\n\n";
}
void build_matrix() {
	make_perturbation();
	cout << "Build matrix:\n";
	vector<ld>ta(n, 0.0);
	vector<ld>tb(N, 0.0);
	B.resize(n, ta);
	At.resize(n, tb);

	for (int i = 0;i < N;i++) {
		ld cur = 1;
		for (int j = 0;j < n;j++) {
			A[i][j] += cur;
			cur *= xs[i];
		}
		y[i] += ys[i];
	}
}
void Gauss_elimination() {
	cout << "Gauss elimination:\n";
	int     p, i, j, k;
	ld  maxEntry, t, r;
	//int n = n;
	for (i = 0; i < n - 1; i++) {
		// Partial pivoting
		maxEntry = fabs(B[i][i]);
		p = i;
		for (k = i; k < n; k++)
			if (fabs(B[k][i]) > maxEntry) {
				p = k;
				maxEntry = fabs(B[k][i]);
			}
		if (p != i) {
			for (j = i; j < n; j++) {
				t = B[p][j];
				B[p][j] = B[i][j];
				B[i][j] = t;
			}
			t = d[p];
			d[p] = d[i];
			d[i] = t;
		}
		//Forward elimination.
		for (k = i + 1; k < n; k++) {
			if (B[k][i] == 0.0) continue;

			r = B[k][i] / B[i][i];
			for (j = i; j < n; j++)
				B[k][j] = B[k][j] - r * B[i][j];
			d[k] = d[k] - r * d[i];
		}
	}
}
void precondition() {
	cout << "Precondition:\n";
	// make At
	for (int i = 0;i < N;i++)
		for (int j = 0;j < n;j++) At[j][i] = A[i][j];
	// make B = At*A
	for (int i = 0;i < n;i++)
		for (int j = 0;j < n;j++)
			for (int k = 0;k < N;k++)
				B[i][j] += At[i][k] * A[k][j];
	// make d = At*y
	for (int i = 0;i < n;i++)
		for (int j = 0;j < N;j++)
			d[i] += At[i][j] * y[j];
}
void verify() {
	cout << "Verify:\n";
	ld newY[N];
	for (int i = 0;i < N;i++) newY[i]=0;
	for (int i = 0;i < N;i++)
		for (int j = 0;j < n;j++)
			newY[i] += A[i][j] * c[j];
	cout << "new Y:\n";
	for (int i = 0;i < N;i++)cout << newY[i] << ' ';
	cout << '\n';
	cout << "origin Y:\n";
	for (int i = 0;i < N;i++)cout << ys[i] << ' ';
	cout << '\n';
	for (int i = 0;i < N;i++)cout << fabs(ys[i]-newY[i]) << ' ';
	cout << '\n';
}
// Compute inner product <a, b>, where a and b are n-dimensional vectors.
ld  inner_product(vector<ld> a, vector<ld> b, int n){
	ld sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += a[i] * b[i];
	return (sum);
}
// Compute a = A*b, where a and b are vectors and A an n x n matrix.
void mtx_vec_mult(vector<ld> &a, vector<vector<ld>> A, vector<ld> b, int n){
	for (int i = 0; i < n; i++) {
		a[i] = 0.0;
		for (int j = 0; j < n; j++)
			a[i] += A[i][j] * b[j];
	}
}

// Normalize a vector. If the vector norm is 0, do nothing.
void normalize_vec(vector<ld>& a, int n){
	ld sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += a[i] * a[i];
	if (sum == 0.0) return;
	sum = sqrt(sum);
	for (int i = 0; i < n; i++)
		a[i] = a[i] / sum;
}

ld norm() {// 1 norm
	ld ret = 0;
	for (int i = 0; i < n; i++)ret += fabs(c[i] - 1.0);
	return ret;
}
// Create a Householder vector for a matrix. This vector
void create_H_vec(vector<vector<ld>>& A, int i, int j, vector<ld>& v, int n1,int n2){
	for (int k = 0; k < i; k++) v[k] = 0.0;
	ld norm2 = 0.0;
	for (int k = i; k < n2; k++) norm2 += A[k][j] * A[k][j];
	norm2 = sqrt(norm2);

	if (A[i][j] >= 0.0) v[i] = A[i][j] + norm2;
	else v[i] = A[i][j] - norm2;

	for (int k = i + 1; k < n2; k++) v[k] = A[k][j];
	print_vector(v, n2, "Householder vector:", n2 != n1);
}
// 2: n1=n2=n 3:n1=n n2=N
void QR_reflect(vector<vector<ld>> &A, vector<ld>& b, int n1,int n2){
	int     i, j, k;
	double vv, vx;
	vector<ld>v(n2);
	// Eliminate each column to make A[][] into U[][].
	// 若是j<n1-1 over-constraint會解不到最後一條 方便起見取j<n1
	for (j = 0; j < n1; j++) {
		// Create vector v[] = A.j + ||A.j||*e1.
		create_H_vec(A, j, j, v, n1,n2);
		vv = inner_product(v, v, n2);
		// Update each column A.k, j<=k<=n-1.
		for (k = j; k < n1; k++) {
			// Compute vx=<v, A.k>
			vx = 0.0;
			for (i = j; i < n2; i++) vx += A[i][k] * v[i];
			// A.k = A.k - 2(<v, A.k>/<v,v>)v ;
			for (i = j; i < n2; i++)
				A[i][k] = A[i][k] - 2.0 * (vx / vv) * v[i];
		}
		// Update b. b = b - 2(<v,b>/<v,v>)v ;
		vx = inner_product(v, b, n2);
		for (i = j; i < n2; i++)
			b[i] = b[i] - 2.0 * (vx / vv) * v[i];
		print_matrix(A, n2, n1, "B:");
		print_vector(b, n2, "d:", n2 != n1);
		cout << "\n-----------------------------------------------------------------";
		cout << "----------------------------------------------------------------------\n";
	}
}
void solve1() {
	Horner();
	build_matrix();
	print_matrix(A, N, n, "A:");
	print_vector(y, N, "y:",1);
	cout << "\n\n";

	precondition();
	print_matrix(B, n, n, "B:");
	print_vector(d, n, "d:", 0);
	cout << "\n\n";

	Gauss_elimination();
	print_matrix(B, n, n, "B:");
	print_vector(d, n, "d:", 0);
	back_substitute(B, d);
	print_vector(c, n, "c:", 0);
	cout << "\n\n";
	// verify
	//verify();
}
void solve2() {
	Horner();
	build_matrix();
	print_matrix(A, N, n, "A:");
	print_vector(y, N, "y:", 1);
	cout << "\n\n";

	precondition();
	print_matrix(B, n, n, "B:");
	print_vector(d, n, "d:", 0);
	cout << "\n\n";

	QR_reflect(B, d, n, n);
	back_substitute(B,d);

	print_vector(c, n, "c:", 0);
	cout << "\n\n";
	// verify
	//verify();
}
void solve3() {
	Horner();
	build_matrix();
	//print_matrix(A, N, n, "A:");
	//print_vector(y, N, "y:", 1);
	cout << "\n\n";

	B = A;
	d = y;
	QR_reflect(B, d, n,N);
	back_substitute(B, d);

	//print_vector(c, n, "c:", 0);
	//cout << "\n\n";
	// verify
	//verify();
}

int main() {
	/*
	solve3();
	cout << "1 norm error is: " << fixed << setprecision(8) << norm();
	*/
	return 0;
}