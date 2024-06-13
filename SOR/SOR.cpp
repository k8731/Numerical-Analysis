#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<algorithm>
#include<cmath>
#include<iomanip>
#include<vector>
#define EPS 1e-6
#define INF 0x7fffffff
typedef long double ld;
using namespace std;
vector<vector<ld> >A;
vector<ld>b, x;
int n;
const int N = 1000;
void create_matrix() {
	x.resize(N, 0);
	A.resize(N, x);
	for (int i = 0; i < N; i++) {
		for (int j = max(0, i - 2); j <= min(N - 1, i + 2); j++)A[i][j] = -1;
		A[i][i] = 4;
	}
}
ld norm() {// 1 norm
	ld error = 0;
	// count residual
	for (int i = 0; i < n; i++) {
		ld ret = 0;
		for (int j = 0; j < n; j++) {// e[i]
			ret += A[i][j] * x[j];
		}
		//error = max(error, abs(ret - b[i]));
		error += abs(ret - b[i]);
	}
	return error;
}
void initial_guess() {
	for (int i = 0; i < n; i++)x[i] = (rand() % 100)/1000.0;
	//for (int i = 0; i < n; i++)cout << x[i] << ' ';
}
int SOR(ld w){
	ld  err = INF;
	// Guess an initial solution.
	for (int i = 0;i < n;i++) x[i] = 0.0;
	//initial_guess();
	int k = 0; // number of iterations
	while (err > EPS) {// Iterate until being converged
		// Compute X(k+1)[i]
		k++;// Increase num of iterations
		for (int i = 0;i < n;i++) {
			ld sum = b[i];
			for (int j = 0;j < n;j++)
				sum = sum - A[i][j] * x[j];// the residual
			x[i] = x[i] + w * sum / A[i][i];// the relaxation
		}
		// Compute the residual
		err = norm();// Compute the norm of the residual
		/*
		cout << k << " iteration ";
		for (int i = 0; i < n; i++)cout << fixed << setprecision(8) << x[i] << ' ';
		cout << '\n';
		*/
	}
	cout <<"need "<< k << " iterations\n";
	return k;
}

int main() {
	//freopen("D:/user/Desktop/學期/數值分析/hw5/data4.in", "r", stdin); 
	//freopen("D:/user/Desktop/學期/數值分析/hw5/data4.out", "w", stdout); 
	n = 40;
	//srand(time(NULL));
	create_matrix();
	n = 20;
	b.clear(), b.resize(n, 0);
	b[0] = b[n - 1] = 2;
	b[1] = b[n - 2] = 1;
	
	for (n = 20; n <= 100; n += 10) {
		cout << n << '\n';
		b.clear(), b.resize(n, 0);
		b[0] = b[n - 1] = 2;
		b[1] = b[n - 2] = 1;
		int best = 1000000000;
		double ome = 1.0;
		for (ld w = 1.8; w <= 1.9+EPS; w += 0.1) {
			cout << w << ' ';
			int it = SOR(w);
			if (it < best) {
				best = it;
				ome = w;
			}
		}
		cout << "best w: " << ome << " \n";
	}
	return 0;
}