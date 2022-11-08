#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <memory.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iomanip>
#define EPS 1e-6
#define ld long double
using namespace std;
struct term {
    ld coef;//coefficient
    int vari;// constant:0,x:1,y:2
    int expo;// exponent
    term() { coef = vari = expo = 0; }
    term(int v, int e, ld c) { vari = v, expo = e, coef = c; }
};
void print(string s,vector<term>&p) {
    cout << s << " = ";
    for (int i = p.size() - 1;i >= 0;i--) {
        // print coef
        if (p[i].coef < 0 || i == p.size() - 1) // no '+'
            cout << p[i].coef;
        else
            cout << '+'<< p[i].coef;
        // print x or y
        if (p[i].vari == 0) {}
        else if (p[i].vari == 1)cout << 'x';
        else cout << 'y';
        // print expo
        if (p[i].expo == 1)continue;
        else cout << '^' << p[i].expo;
    }
    cout << "\n";
}
vector<term>partial_differential(vector<term>& p,int v) {
    vector<term>newp;
    for (int i = 0;i < p.size();i++) {
        if (p[i].vari == 0)continue;// �`�ƷL��=0
        if (p[i].vari != v)continue;// �����`�ƷL��
        if (p[i].expo == 1)// �@����ƷL��=�`��
            newp.push_back(term(0, 1, p[i].coef));
        else // �G����ƥH�W
            newp.push_back(term(p[i].vari, p[i].expo-1, p[i].coef*p[i].expo));
    }
    return newp;
}
ld func(vector<term>& p, ld x, ld y) {// �p���ƭ�
    ld ret = 0;
    for (int i = 0;i < p.size();i++) {
        if (p[i].vari == 0)ret += p[i].coef;// constant
        else if (p[i].vari == 1) {
            ld x_ = 1;
            for (int j = 0;j < p[i].expo;j++)x_ *= x;// x^(p[i].expo)
            ret += p[i].coef * x_;
        }
        else {
            ld y_ = 1;
            for (int j = 0;j < p[i].expo;j++)y_ *= y;// y^(p[i].expo)
            ret += p[i].coef * y_;
        }
    }
    return ret;
}

ld error_1(ld e0, ld e1) {// 1 norm 
    return abs(e0) + abs(e1);
}
ld error_2(ld e0,ld e1) {// 2 norm 
    return sqrt(e0 * e0 + e1 * e1);
}
ld error_inf(ld e0, ld e1) {// inf norm 
    return max(abs(e0), abs(e1));
}
ld C(ld fx,ld fy,ld gx,ld gy) {// condition number
    // count eigenvalues
    ld eigen1 = (fx + gy + sqrt((fx + gy) * (fx + gy) - 4 * (fx * gy - gx * fy))) / 2.0;
    ld eigen2 = (fx + gy - sqrt((fx + gy) * (fx + gy) - 4 * (fx * gy - gx * fy))) / 2.0;
    eigen1 = abs(eigen1), eigen2 = abs(eigen2);
    // C=abs(eigen_max)/abs(eigen_min)
    if (eigen1 > eigen2)return eigen1 / eigen2;
    else return eigen2 / eigen1;
}
vector<term>f, g, fx, fy, gx, gy;
double c;
int Newton_method_2D(ld x0,ld y0) {
    ld x = x0, y = y0;
    ld h,k,h0,k0;
    int cycle = 0;
    cout << "start!!\n";
    cout << "initial point: [" << x0 << ',' << y0 << "]\n"; 
    cout << " i           x                    y                   error                   C" << '\n';
    //cout <<" i           x                y                f(x,y)            g(x,y)         error           C"<<'\n';
    //cout << " i           x                y                f(x,y)            g(x,y)         error" << '\n';
    //cout <<" i           x                y               fx               fy               gx              gy"<<'\n';
    cout << fixed << setw(2) << cycle << "  ";
    cout << fixed << setprecision(12) << setw(20) << x << ' ' << setw(20) << y << '\n';
    //cout << fixed << setprecision(12) << setw(16) << x << ' ' << setw(16) << y << ' ' << setw(16) << func(f, x, y) << ' ' << setw(16) << func(g, x, y) << '\n';
    //cout << fixed << setprecision(12) << setw(16) << x << ' ' << setw(16) << y << ' ';
    //cout << fixed << setprecision(12) << setw(16) << func(fx, x, y) << ' ' << setw(16) << func(fy, x, y) << ' ' << setw(16) << func(gx, x, y) << ' ' << setw(16) << func(gy, x, y) << '\n';
    //cout << fixed << setprecision(10) << C(func(fx, x, y), func(fy, x, y), func(gx, x, y), func(gy, x, y)) << '\n';
    //c = C(func(fx, x, y), func(fy, x, y), func(gx, x, y), func(gy, x, y));
    do{
        // det=fy*gx-fx*gy
        // h=(1/det)*(gy*f-fy*g)
        // k=(1/det)*(-gx*f+fx*g)
        cycle++;
        ld det = func(fy, x, y) * func(gx, x, y) - func(fx, x, y) * func(gy, x, y);
        if (abs(det) <= 1e-12) {// det=0 �L�k��
            cout << "divergence!!!\n\n";
            return -1;// �o��
        } 
        h = 1 / det * (func(gy, x, y) * func(f, x, y) - func(fy, x, y) * func(g, x, y));
        k = 1 / det * (-func(gx, x, y) * func(f, x, y) + func(fx, x, y) * func(g, x, y));
        x += h, y += k;
        //cout << fixed << setw(2)<<cycle << "  ";
        cout << fixed << setw(2) << cycle << "  ";
        cout << fixed << setprecision(12) << setw(20) << x << ' ' << setw(20) << y << ' ';
        //cout << fixed << setprecision(12) << setw(16) << x << ' ' << setw(16) << y << ' ' << setw(16) << func(f, x, y) << ' ' << setw(16) << func(g, x, y) << ' ';
        //cout << fixed << setprecision(12) << setw(16) << x << ' ' << setw(16) << y << ' ';
        //cout << fixed << setprecision(12) << setw(16) << func(fx, x, y) << ' ' << setw(16) << func(fy, x, y) << ' ' << setw(16) << func(gx, x, y) << ' ' << setw(16) << func(gy, x, y) <<' ';
        cout<< fixed << setprecision(12) << setw(20)<< error_1(h, k)<<' ';
        //cout << '\n';
        cout << fixed << setprecision(12) << setw(20)<<C(func(fx, x, y), func(fy, x, y), func(gx, x, y), func(gy, x, y)) << '\n';
        //if (cycle == 1)c = C(func(fx, x, y), func(fy, x, y), func(gx, x, y), func(gy, x, y));
    } while (error_2(h, k) >= EPS);
    cout << "finished!!\n";
    cout << "used " << cycle << " times.\n\n";
    return cycle;
}
int test_oval(int x, int y) {// �b���Wor����or�~��
    int ret = 4 * x * x + 9 * y * y;
    if (ret > 36)return 1;// ���~
    else if (ret == 36)return 0;// ���W
    else return -1;// ��ꤺ
}
int test_parabola(int x, int y) {// �b�ߪ��u����or�~��
    int ret = x * x - 2;
    if (ret > y)return 1;// �ߪ��u�~
    else if (ret == y)return 0;// �ߪ��u��
    else return -1;// �ߪ��u��
}
struct point {
    double x, y;
    int iter;
    point() { x = y = iter = 0; }
    point(double x_, double y_, int it) { x = x_, y = y_, iter = it; }
};
int main() {
    
    /*
    freopen("file.in", "r", stdin); // Ū file.in ��
    freopen("file.out", "w", stdout); // �g�J file.out ��
    */
    
    /*
    freopen("D:/user/Desktop/�Ǵ�/�ƭȤ��R/hw1/animation4.in", "r", stdin); // Ū animation.in ��
    freopen("D:/user/Desktop/�Ǵ�/�ƭȤ��R/hw1/animation4.out", "w", stdout); // �g�J animation.out ��
    */

    // f(x,y)=-1+(1/9)x^2+(1/4)y^2
    f.push_back(term(0,1,-1));
    f.push_back(term(2, 2, 1 / ld(4)));
    f.push_back(term(1, 2, 1 / ld(9)));
    // g(x,y)=-2+x^2-y
    g.push_back(term(0, 1, -2));
    g.push_back(term(2, 1, -1));
    g.push_back(term(1, 2, 1));
    print("f(x,y)",f);
    print("g(x,y)",g);
    // ���L��f�My
    fx = partial_differential(f, 1);
    fy = partial_differential(f, 2);
    gx = partial_differential(g, 1);
    gy = partial_differential(g, 2);
    print("fx", fx);
    print("fy", fy);
    print("gx", gx);
    print("gy", gy);
    cout << '\n';
    // x=-4~+4 y=-3~+3��������Ʈ�l�I
    vector<point>div;// �|�o�����I
    vector<vector<point>>group(5);// �����Igroup[g]:group=g��iteration����
    vector<double>speed(5);
    vector<point>c_small, c_big,c_nan;
    vector<point>ps;

    /*
    for (double x = -4;x <= 4;x += 1) {
        for (double y = -3;y <= 3;y += 1) {
            int iteration = Newton_method_2D(x, y);
            if (iteration == -1)div.push_back(point( x,y,0 ));
            else {
                int oval = test_oval(x, y);
                int para = test_parabola(x, y);
                if (oval == -1 && para == -1)group[0].push_back(point(x, y, iteration));
                else if (oval == -1 && para == 1)group[1].push_back(point(x, y, iteration));
                else if (oval == 1 && para == -1)group[2].push_back(point(x, y, iteration));
                else if (oval == 1 && para == 1)group[3].push_back(point(x, y, iteration));
                else group[4].push_back(point(x, y, iteration));
                if (c <= 1.9&&c>=0)c_small.push_back(point(x, y, iteration));
                else if(c>=0)c_big.push_back(point(x, y, iteration));
                else c_nan.push_back(point(x, y, iteration));
                ps.push_back(point(x, y, iteration));
            }
        }
    }
    */
    // condition number analysis
    /*
    cout << "small\n";
    for (point p : c_small) {
        cout << p.x << ' ' << p.y << ' ' << p.iter << " \n";
    }
    cout << "big\n";
    for (point p : c_big) {
        cout << p.x << ' ' << p.y << ' ' << p.iter << " \n";
    }
    cout << "nan\n";
    for (point p : c_nan) {
        cout << p.x << ' ' << p.y << ' ' << p.iter << " \n";
    }
    */
    /*
    
    for (int i = 0;i < 5;i++) {
        cout << i << '\n';
        for (point p : group[i])cout << p.x << ' ' << p.y << '\n';
        cout << '\n';
    }
    cout << "div:\n";
    for (point p : div)cout << p.x << ' ' << p.y << '\n';
    cout << '\n';
    */
    // choose �S���I
    //Newton_method_2D(0.000001, 2);
    // ���R�ϰ��I���ĳt��
    /*
    for (int i = 0;i < 5;i++) {
        cout << i <<" group" <<": ";
        double cnt = 0;
        for (point p : group[i])cnt += p.iter;
        speed[i] = cnt / group[i].size();
        cout << speed[i] << '\n';
    }
    */
    Newton_method_2D(0.00001,0.00001);
    return 0;
}