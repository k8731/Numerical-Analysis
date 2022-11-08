#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#define EPS 1e-15
using namespace std;
typedef long double ld;
const double π = acos(-1);
ld angle[50], x[50], y[50], coefx[50], coefy[50], qx[400], qy[400], px[400], py[400];
ld error_2_group[80], error_inf_group[80];
// θ Δ π
int N = 8;// N+1個sample points 牛頓：N次多項式
ld r = 10.0;// radius
double range = 2 * π;// 0~range
void print_table_θ(int s, int t) {
	// print angle
	cout << "                   |  ";
	for (int i = s; i <= t; i++)
		cout << "  θ" << setw(2)<<i << "      ";
	cout << "\nθ[i] = i*Δθ     |";
	cout << "-----------------------------------------------------------------------------------------------------------------------------------";
	cout << "\n                   | ";
	for (int i = s; i <= t; i++) 
		cout << setw(3) << i << "π/" << setw(2) << N/2 << "    ";
	cout << "\n\n";
}
void print_table_x(int s, int t) {
	// print x
	cout << "                   |  ";
	for (int i = s; i <= t; i++)
		cout << "   X" << setw(2)<< i << "      ";
	cout << "\nX[i] = r*cos(θ[i])|";
	cout << "-----------------------------------------------------------------------------------------------------------------------------------";
	cout << "\n                   | ";
	for (int i = s; i <= t; i++) 
		cout << fixed << setprecision(5) << setw(8) << x[i] << "    ";
	cout << "\n\n";
}
void print_table_y(int s, int t) {
	// print y
	cout << "                   |  ";
	for (int i = s; i <= t; i++)
		cout << "   Y" << setw(2)<<i << "      ";
	cout << "\nY[i] = r*sin(θ[i])|";
	cout << "-----------------------------------------------------------------------------------------------------------------------------------";
	cout << "\n                   | ";
	for (int i = s; i <= t; i++)
		cout << fixed << setprecision(5) << setw(8) << y[i] << "    ";
	cout << "\n\n";
}
ld func_by_hornor_x(ld p) {
	ld f = coefx[N];
	for (int i = N - 1; i >= 0; i--)
		f = f * (p - angle[i]) + coefx[i];
	return f;
}
ld func_by_hornor_y(ld p) {
	ld f = coefy[N];
	for (int i = N - 1; i >= 0; i--)
		f = f * (p - angle[i]) + coefy[i];
	return f;
}
void print_table_8() {
	print_table_θ(0, N);
	cout << '\n';
	print_table_x(0, N);
	cout << '\n';
	print_table_y(0, N);
}
void print_table_16() {
	print_table_θ(0, 8);
	print_table_θ(9, 16);
	cout << '\n';
	print_table_x(0, 8);
	print_table_x(9, 16);
	cout << '\n';
	print_table_y(0, 8);
	print_table_y(9, 16);
}
void print_table_24() {
	print_table_θ(0, 8);
	print_table_θ(9, 16);
	print_table_θ(17, 24);
	cout << '\n';
	print_table_x(0, 8);
	print_table_x(9, 16);
	print_table_x(17, 24);
	cout << '\n';
	print_table_y(0, 8);
	print_table_y(9, 16);
	print_table_y(17, 24);
}
void print_table_40() {
	print_table_θ(0, 10);
	print_table_θ(11, 20);
	print_table_θ(21, 30);
	print_table_θ(31, 40);
	cout << '\n';
	print_table_x(0, 10);
	print_table_x(11, 20);
	print_table_x(21, 30);
	print_table_x(31, 40);
	cout << '\n';
	print_table_y(0, 10);
	print_table_y(11, 20);
	print_table_y(21, 30);
	print_table_y(31, 40);
}
void count() {
	//cout << "sample point:\n";
	for (int i = 0; i <= N; i++) {
		angle[i] = i * range / N;
		if (fabs(angle[i]) < EPS)angle[i] = 0;
		x[i] = r * cos(angle[i]);
		y[i] = r * sin(angle[i]);
		//cout << i << ' ' << fixed << setprecision(16) << x[i] << ' ' << y[i] << " \n";
	}
}
void random_sample() {
	//cout << "sample point:\n";
	angle[0] = 0;
	angle[1] = 10 * 2 * π / 360;
	angle[2] = 30 * 2 * π / 360;
	angle[3] = 50 * 2 * π / 360;
	angle[4] = 100 * 2 * π / 360;
	angle[5] = 120 * 2 * π / 360;
	angle[6] = 180 * 2 * π / 360;
	angle[7] = 200 * 2 * π / 360;
	angle[8] = 360 * 2 * π / 360;
	for (int i = 0; i <= N; i++) {
		x[i] = r * cos(angle[i]);
		y[i] = r * sin(angle[i]);
		//cout << i << ' ' << fixed << setprecision(16) << x[i] << ' ' << y[i] << " \n";
	}
}
void Divided_Difference_x() {
	for (int i = 0; i <= N; i++)coefx[i] = x[i];
	for (int t = 1; t <= N; t++) {
		for (int i = N; i >= t; i--) {
			coefx[i] = (coefx[i] - coefx[i - 1]) / (angle[i] - angle[i - t]);
		}
	}
	//cout << "x coefficients:\n";
	//for (int i = 0; i <= N; i++)cout << "c" << i << ":" << coefx[i] << '\n';
}
void Divided_Difference_y() {
	for (int i = 0; i <= N; i++)coefy[i] = y[i];
	for (int t = 1; t <= N; t++) {
		for (int i = N; i >= t; i--) {
			coefy[i] = (coefy[i] - coefy[i - 1]) / (angle[i] - angle[i - t]);
		}
	}
	//cout << "y coefficients:\n";
	//for (int i = 0; i <= N; i++)cout << "c" << i << ":" << coefy[i] << '\n';
}
ld error_2() {// 2 norm 
	ld sum = 0;
	ld ret = 0;
	for (int i = 0; i < 360; i++) {
		sum += ((px[i] - qx[i]) * (px[i] - qx[i]) + (py[i] - qy[i]) * (py[i] - qy[i]));
		ret += ((px[i] - qx[i]) * (px[i] - qx[i]) + (py[i] - qy[i]) * (py[i] - qy[i]));
		if (i&&i % 5 == 0) {
			error_2_group[i / 5 - 1] = sqrt(ret);
			ret = 0;
		}
	}
	error_2_group[360 / 5-1] = sqrt(ret);
	return sqrt(sum);
}
ld error_inf() {// inf norm 
	ld sum = 0;
	ld ret = 0;
	for (int i = 0; i < 360; i++) {
		sum = max(sum, sqrt((px[i] - qx[i]) * (px[i] - qx[i]) + (py[i] - qy[i]) * (py[i] - qy[i])));
		ret = max(ret, sqrt((px[i] - qx[i]) * (px[i] - qx[i]) + (py[i] - qy[i]) * (py[i] - qy[i])));
		if (i && i % 5 == 0) {
			error_inf_group[i / 5 - 1] = ret;
			ret = 0;
		}
	}
	error_inf_group[360 / 5 - 1] = ret;
	return sum;
}
int main() {
	//freopen("D:/user/Desktop/學期/數值分析/hw2/code/random4_inf_norm.in", "r", stdin); 
	//freopen("D:/user/Desktop/學期/數值分析/hw2/code/random4_inf_norm.out", "w", stdout); 
	//for (; N <= 40; N++) {
		//cout << N << "\n";
	//count();// 計算取樣點
	random_sample();
	//print_table_8();
	Divided_Difference_x();// 計算差分
	Divided_Difference_y();// 計算差分
	double delta_t = 2 * π / 360;
	// Newton interpolation
	//cout << "Newton interpolation:\n";
	for (int i = 0; i <= 360; i++) {
		ld ang = i * delta_t;
		qx[i] = func_by_hornor_x(ang);
		qy[i] = func_by_hornor_y(ang);
		//cout << i << fixed << setprecision(16) << ' ' << qx[i] << ' ' << qy[i] << ' ' << '\n';
	}
	//}
	// real function
	//cout << "\nreal function:\n";
	for (int i = 0; i <= 360; i++) {
		px[i] = r * cos(i * delta_t);
		py[i] = r * sin(i * delta_t);
		//cout << i << fixed << setprecision(16) << ' ' << px[i] << ' ' << py[i] << ' ' << '\n';
	}
	// 計算誤差
	//error_2();
	//error_inf();
	//cout << "error 2:" << fixed << setprecision(20)<<error_2() << '\n';
	//cout << "error inf:" << fixed << setprecision(20)<< error_inf() << '\n';
	// 誤差最大和最小的區域分析(五個點/四個間隔為一單位)
	//cout << "error 2:\n";
	//for (int i = 0; i <= 360 / 5 - 1; i++)cout <<fixed<<setprecision(20)<< i * 5 << " - " << i * 5 + 4 << " : " << error_2_group[i] << " .\n";
	//cout << "error inf:\n";
	//for (int i = 0; i <= 360 / 5 - 1; i++)cout <<fixed<<setprecision(20)<< i * 5 << " - " << i * 5 + 4 << " : " << error_inf_group[i] << " .\n";
	return 0;
}