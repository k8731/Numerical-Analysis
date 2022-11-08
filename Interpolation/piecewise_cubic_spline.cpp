#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#define EPS 1e-15
using namespace std;
typedef long double ld;
const double π = acos(-1);
struct term {
	ld coef;//coefficient
	int expo;// exponent
	term() { coef = expo = 0; }
	term(int e, ld c) { expo = e, coef = c; }
};
void print(string s, vector<term>& p,int op) {// op=1:x op=2:y
	cout << s << " = ";
	for (int i = p.size() - 1;i >= 0;i--) {
		// print coef
		if (p[i].coef < 0 || i == p.size() - 1) // no '+'
			cout << p[i].coef;
		else
			cout << '+' << p[i].coef;
		// print x or y
		if (p[i].expo != 0) {
			if (op == 1)cout << 'x';
			else cout << 'y';
		}
		// print expo
		if (p[i].expo <= 1)continue;
		else cout << '^' << p[i].expo;

	}
	cout << "\n";
}
ld func(vector<term>& p, ld an) {// 計算函數值
	ld ret = 0;
	for (int i = 0;i < p.size();i++) {
		if (p[i].expo == 0)ret += p[i].coef;// constant
		else {
			ld x_ = 1;
			for (int j = 0;j < p[i].expo;j++)x_ *= an;// x^(p[i].expo)
			ret += p[i].coef * x_;
		}
		
	}
	return ret;
}
vector<term>xf[6],yf[6];
ld angle[10], x[10], y[10], qx[400], qy[400], px[400], py[400];
ld coefx[6][4], coefy[6][4];// coef[i]:pi開始的多項式
ld error_2_group[80], error_inf_group[80];
// θ Δ π
int N = 8;// N+1個sample points
ld r = 10.0;// radius
double range = 2.0 * π;// 0~range
void count() {
	//cout << "sample point:\n";
	for (int i = 0; i <= N; i++) {
		angle[i] = ld(i) * range / ld(N);
		if (fabs(angle[i]) < EPS)angle[i] = 0;
		x[i] = r * cos(angle[i]);
		y[i] = r * sin(angle[i]);
		//cout << angle[i] << " .\n";
	}
}
void Divided_Difference_x(int p0) {
	for (int i = 0; i <= 3; i++)
		coefx[p0][i] = x[p0 + i];
	for (int t = 1; t <= 3; t++) {
		for (int i = 3; i >= t; i--) {
			coefx[p0][i] = (coefx[p0][i] - coefx[p0][i - 1]) / (angle[p0 + i] - angle[p0 + i - t]);
			//if (fabs(coefx[p0][i]) < EPS)coefx[p0][i] = 0;
		}
	}
}
void Divided_Difference_y(int p0) {
	for (int i = 0; i <= 3; i++)coefy[p0][i] = y[p0 + i];
	for (int t = 1; t <= 3; t++) {
		for (int i = 3; i >= t; i--) {
			coefy[p0][i] = (coefy[p0][i] - coefy[p0][i - 1]) / (angle[p0 + i] - angle[p0 + i - t]);
			//if (fabs(coefy[p0][i]) < EPS)coefy[p0][i] = 0;
		}
	}
}
ld func_by_hornor_x(int s,ld p) {
	ld f = coefx[s][3];
	for (int i = 3 - 1; i >= 0; i--)
		f = f * (p - angle[s+i]) + coefx[s][i];
	return f;
}
ld func_by_hornor_y(int s, ld p) {
	ld f = coefy[s][3];
	for (int i = 3 - 1; i >= 0; i--)
		f = f * (p - angle[s+i]) + coefy[s][i];
	return f;
}
ld error_2() {// 2 norm 
	ld sum = 0;
	ld ret = 0;
	for (int i = 0; i < 300; i++) {
		sum += ((px[i] - qx[i]) * (px[i] - qx[i]) + (py[i] - qy[i]) * (py[i] - qy[i]));
		ret += ((px[i] - qx[i]) * (px[i] - qx[i]) + (py[i] - qy[i]) * (py[i] - qy[i]));
		if (i && i % 5 == 0) {
			error_2_group[i / 5 - 1] = sqrt(ret);
			ret = 0;
		}
	}
	error_2_group[300 / 5 - 1] = sqrt(ret);
	return sqrt(sum);
}
ld error_inf() {// inf norm 
	ld sum = 0;
	ld ret = 0;
	for (int i = 0; i < 300; i++) {
		sum = max(sum, sqrt((px[i] - qx[i]) * (px[i] - qx[i]) + (py[i] - qy[i]) * (py[i] - qy[i])));
		ret = max(ret, sqrt((px[i] - qx[i]) * (px[i] - qx[i]) + (py[i] - qy[i]) * (py[i] - qy[i])));
		if (i && i % 5 == 0) {
			error_inf_group[i / 5 - 1] = ret;
			ret = 0;
		}
	}
	error_inf_group[300 / 5 - 1] = ret;
	return sum;
}
ld derivative_x(int p0,ld t) {
	// f'(x)=c3[3x^2-2x(x0+x1+x2)+(x0x1+x0x2+x1x2)]+c2[2x-(x0+x1))+c1
	ld c3 = coefx[p0][3], c2 = coefx[p0][2], c1 = coefx[p0][1], c0 = coefx[p0][0];
	ld x0 = angle[p0 + 0], x1 = angle[p0 + 1], x2 = angle[p0 + 2], x3 = angle[p0 + 3];
	return c3 * ((t - x1) * (t - x2) + (t - x0) * (t - x1) + (t - x0) * (t - x2)) + c2 * (2 * t - (x0 + x1)) + c1;
}
ld derivative_y(int p0, ld t) {
	// f'(x)=c3[3x^2-2x(x0+x1+x2)+(x0x1+x0x2+x1x2)]+c2[2x-(x0+x1))+c1
	ld c3 = coefy[p0][3], c2 = coefy[p0][2], c1 = coefy[p0][1], c0 = coefy[p0][0];
	ld x0 = angle[p0 + 0], x1 = angle[p0 + 1], x2 = angle[p0 + 2], x3 = angle[p0 + 3];
	return c3 * ((t - x1) * (t - x2) + (t - x0) * (t - x1) + (t - x0) * (t - x2)) + c2 * (2 * t - (x0 + x1)) + c1;
}
int main() {
	//freopen("D:/user/Desktop/學期/數值分析/hw2/piecewise.in", "r", stdin); 
	//freopen("D:/user/Desktop/學期/數值分析/hw2/piecewise.out", "w", stdout); 
	//freopen("D:/user/Desktop/學期/數值分析/hw2/piecewise_p.in", "r", stdin); 
	//freopen("D:/user/Desktop/學期/數值分析/hw2/piecewise_p.out", "w", stdout); 
	//freopen("D:/user/Desktop/學期/數值分析/hw2/piecewise_N=8_2norm.in", "r", stdin);
	//freopen("D:/user/Desktop/學期/數值分析/hw2/piecewise_N=8_2norm.out", "w", stdout);
	//freopen("D:/user/Desktop/學期/數值分析/hw2/piecewise_N=8_inf_norm.in", "r", stdin);
	//freopen("D:/user/Desktop/學期/數值分析/hw2/piecewise_N=8_inf_norm.out", "w", stdout);
	count();
	//cout << "Piecewise method:\n";

	for (int i = 0;i <= 5;i++) {
		Divided_Difference_x(i);
		Divided_Difference_y(i);
		int k = 0;
		double j;
		for (k=0,j = (i + 1) * 45.0; k <= 50; j += 0.9, k++) {// 270度->300個點
			double ang = j * 2 * π / 360;
			qx[k + i * 50] = func_by_hornor_x(i, ang);
			qy[k + i * 50] = func_by_hornor_y(i, ang);
			px[k + i * 50] = r * cos(ang);
			py[k + i * 50] = r * sin(ang);
		}
		
		/*
		cout << "x: ";
		for (int j = 0;j <= 3;j++)cout << fixed << setprecision(16) << coefx[i][j] << ' ';
		cout << "\ny: ";
		for (int j = 0;j <= 3;j++)cout << fixed << setprecision(16) << coefy[i][j] << ' ';
		cout << '\n';
		*/
		
	}
	// 計算誤差
	//error_2();
	//error_inf();
	//cout << "error 2:" << fixed << setprecision(20) << error_2() << '\n';
	//cout << "error inf:" << fixed << setprecision(20) << error_inf() << '\n';
	// 誤差最大和最小的區域分析(五個點/四個間隔為一單位)
	//cout << "error 2:\n";
	//for (int i = 0; i <= 300 / 5 - 1; i++)cout << fixed << setprecision(20) << i * 5 << " - " << i * 5 + 4 << " : " << error_2_group[i] << " .\n";
	//cout << "error inf:\n";
	//for (int i = 0; i <= 300 / 5 - 1; i++)cout << fixed << setprecision(20) << i * 5 << " - " << i * 5 + 4 << " : " << error_inf_group[i] << " .\n";
	// 測試是否C0-continuity
	cout << "p1" << ":\n";
	cout << fixed << setprecision(12) << "X(θ" << 1 << ") = " << func_by_hornor_x(0, angle[1]) << " Y(θ" << 1 << ") = " << func_by_hornor_y(0, angle[1]) << '\n';
	cout << fixed << setprecision(12) << "X(θ" << 1 << ") = " << func_by_hornor_x(1, angle[1]) << " Y(θ" << 1 << ") = " << func_by_hornor_y(1, angle[1]) << '\n';
	cout << "p2" << ":\n";
	cout << fixed << setprecision(12) << "X(θ" << 2 << ") = " << func_by_hornor_x(0, angle[2]) << " Y(θ" << 2 << ") = " << func_by_hornor_y(0, angle[2]) << '\n';
	cout << fixed << setprecision(12) << "X(θ" << 2 << ") = " << func_by_hornor_x(2, angle[2]) << " Y(θ" << 2 << ") = " << func_by_hornor_y(2, angle[2]) << '\n';
	for (int ps = 3; ps <= 5; ps++) {
		cout << "p" << ps << ":\n";
		cout << fixed << setprecision(12) << "X(θ" << ps << ") = " << func_by_hornor_x(ps - 3, angle[ps]) << " Y(θ" << ps << ") = " << func_by_hornor_y(ps - 3, angle[ps]) << '\n';
		cout << fixed << setprecision(12) << "X(θ" << ps << ") = " << func_by_hornor_x(ps, angle[ps]) << " Y(θ" << ps << ") = " << func_by_hornor_y(ps, angle[ps]) << '\n';
	}
	cout << "p6" << ":\n";
	cout << fixed << setprecision(12) << "X(θ" << 6 << ") = " << func_by_hornor_x(3, angle[6]) << " Y(θ" << 6 << ") = " << func_by_hornor_y(3, angle[6]) << '\n';
	cout << fixed << setprecision(12) << "X(θ" << 6 << ") = " << func_by_hornor_x(5, angle[6]) << " Y(θ" << 6 << ") = " << func_by_hornor_y(5, angle[6]) << '\n';
	cout << "p7" << ":\n";
	cout << fixed << setprecision(12) << "X(θ" << 7 << ") = " << func_by_hornor_x(4, angle[7]) << " Y(θ" << 7 << ") = " << func_by_hornor_y(4, angle[7]) << '\n';
	cout << fixed << setprecision(12) << "X(θ" << 7 << ") = " << func_by_hornor_x(5, angle[7]) << " Y(θ" << 7 << ") = " << func_by_hornor_y(5, angle[7]) << '\n';
	
	// 測試是否C1-continuity
	cout << "p1"<< ":\n";
	cout << fixed << setprecision(12) << "X'(θ" << 1 << ") = " << derivative_x(0, angle[1]) << " Y'(θ" << 1 << ") = " << derivative_y(0, angle[1]) << '\n';
	cout << fixed << setprecision(12) << "X'(θ" << 1 << ") = " << derivative_x(1, angle[1]) << " Y'(θ" << 1 << ") = " << derivative_y(1, angle[1]) << '\n';
	cout << "p2" << ":\n";
	cout << fixed << setprecision(12) << "X'(θ" << 2 << ") = " << derivative_x(0, angle[2]) << " Y'(θ" << 2 << ") = " << derivative_y(0, angle[2]) << '\n';
	cout << fixed << setprecision(12) << "X'(θ" << 2 << ") = " << derivative_x(2, angle[2]) << " Y'(θ" << 2 << ") = " << derivative_y(2, angle[2]) << '\n';
	for (int ps = 3; ps <= 5; ps++) {
		cout << "p" << ps << ":\n";
		cout << fixed << setprecision(12) << "X'(θ" << ps << ") = " << derivative_x(ps - 3, angle[ps]) << " Y'(θ" << ps << ") = " << derivative_y(ps - 3, angle[ps]) << '\n';
		cout << fixed << setprecision(12) << "X'(θ" << ps << ") = " << derivative_x(ps, angle[ps]) << " Y'(θ" << ps << ") = " << derivative_y(ps, angle[ps]) << '\n';
	}
	cout << "p6" << ":\n";
	cout << fixed << setprecision(12) << "X'(θ" << 6 << ") = " << derivative_x(3, angle[6]) << " Y'(θ" << 6 << ") = " << derivative_y(3, angle[6]) << '\n';
	cout << fixed << setprecision(12) << "X'(θ" << 6 << ") = " << derivative_x(5, angle[6]) << " Y'(θ" << 6 << ") = " << derivative_y(5, angle[6]) << '\n';
	cout << "p7" << ":\n";
	cout << fixed << setprecision(12) << "X'(θ" << 7 << ") = " << derivative_x(4, angle[7]) << " Y'(θ" << 7 << ") = " << derivative_y(4, angle[7]) << '\n';
	cout << fixed << setprecision(12) << "X'(θ" << 7 << ") = " << derivative_x(5, angle[7]) << " Y'(θ" << 7 << ") = " << derivative_y(5, angle[7]) << '\n';
	return 0;
}