#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<cmath>
using namespace std;
typedef long double ld;
ld factor[5][3] = { 0 };
vector<vector<ld>>coef(5),root(5);
ld func(ld x, ld y) {
	return sin(2 * x) * sin(2 * x) - cos(3 * y) * cos(3 * y) + 1.0;
}
ld func_mapping(ld x, ld y,double a,double b,double c,double d) {
	ld s = pow(sin(2.0 * ((b - a) * x + b + a) / 2.0), 2.0);
	ld t = pow(cos(3.0 * ((d - c) * y + d + c) / 2.0), 2.0);
	return (s-t + 1.0) * (b - a) * (d - c) / 4.0;
}
ld func_mapping2(ld x, ld y, int a, int b, int c, int d) {
	return (8 * pow(double((b-a) * x + a+b)/2.0, 4.0) + 6 * pow(double((d-c) * y + d+c)/2.0, 4.0) + 4) *(b - a) * (d - c) / 4.0;
}
void pre_create() {
	coef[2].push_back(1);
	coef[2].push_back(1);
	coef[3].push_back(0.55555555555555);
	coef[3].push_back(0.88888888888888);
	coef[3].push_back(0.55555555555555);
	coef[4].push_back(0.34785484513745);
	coef[4].push_back(0.65214515486255);
	coef[4].push_back(0.65214515486255);
	coef[4].push_back(0.34785484513745);
	root[2].push_back(-0.57735026918963);
	root[2].push_back(0.57735026918963);
	root[3].push_back(-0.77459666924148);
	root[3].push_back(0.0);
	root[3].push_back(0.77459666924148);
	root[4].push_back(-0.86113631159405);
	root[4].push_back(-0.33998104358486);
	root[4].push_back(0.33998104358486);
	root[4].push_back(0.86113631159405);
}
ld count(int cell,int n) {
	ld ans = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (double x = 0; x < 4; x += 4.0 / cell)
				for (double y = 0; y < 4; y += 4.0 / cell)
					ans += coef[n][i] * coef[n][j] * func_mapping(root[n][i], root[n][j], x, x + 4.0 / cell, y, y + 4.0 / cell);
	return ans ;
}
ld absolute_error(ld real, ld comp) {
	return fabs(real - comp);
}
ld relative_error(ld real, ld comp) {
	return fabs(real - comp)/fabs(real);
}
void print_data(int cell, int n) {
	ld real_integral = 16.0 - sin(24.0) / 3.0 - sin(16.0) / 2.0;
	ld data = count(cell, n);
	int t = (int)(log(cell+1)/log(4.0));
	cout << fixed << setprecision(12) << "cell = " << cell << ", N = " << n << ": " << data <<' ';
	cout << "absolute error = "  << absolute_error(data, real_integral) << ' ';
	cout << "relative error = " << relative_error(data, real_integral) << " \n";
	factor[t][n - 2] = absolute_error(data, real_integral);
}
void factor_analysis() {
	ld real_integral = 16.0 - sin(24.0) / 3.0 - sin(16.0) / 2.0;
	
	for (int cell = 0; cell < 5; cell++) {
		cout << "cell " << (int)pow(4, cell) << ": ";
		for (int n = 1; n <= 2; n++) {
			//if (abs(factor[cell][n - 1] - factor[cell][n]) < 1e-9)cout << "inf ";
			cout << factor[cell][n-1] / factor[cell][n] << ' ';
			//cout << factor[cell][n - 2] << ' ';
		}
		cout << '\n';
			
	}
	for (int n = 0; n <= 2; n++) {
		cout << "N " << n+2 << ": ";
		for (int cell = 1; cell < 5; cell++) {
			//if (abs(factor[cell-1][n] - factor[cell][n]) < 1e-9)cout << "inf ";
			cout << factor[cell-1][n] / factor[cell][n] << ' ';
			//cout << factor[cell][n - 2] << ' ';
		}
		cout << '\n';

	}
	
}
ld trun_count(double a,double b,int n) {
	ld ret = pow(b - a, 2*n+1);
	for (int i = n; i >= 1; i--) {
		ret *= pow(i, 4);
	}
	ret /= (2 * n + 1);
	for (int i = 2*n; i >= 1; i--) {
		ret /= pow(i, 3);
	}

	return ret;
}
int main() {
	//freopen("D:/user/Desktop/學期/數值分析/hw3/data.in", "r", stdin); 
	//freopen("D:/user/Desktop/學期/數值分析/hw3/data.out", "w", stdout); 
	ld real_integral = 16.0 - sin(24.0) / 3.0 - sin(16.0) / 2.0;
	//ld real_integral = 11532.8;
	cout << "actual value: " << fixed << setprecision(20) << real_integral << '\n';
	pre_create();
	for (int cell = 1; cell <= 256; cell *= 4) {
		for (int n = 2; n <= 4; n++)
			print_data(cell, n);
	}
	cout << "factor analysis:\n";
	factor_analysis();
	return 0;
}