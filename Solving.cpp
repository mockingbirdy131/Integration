#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <random>
#include <cmath>
using namespace std;

const double a = 0.5;       
const double b = 1;   
const int K = 1;

double func(double x);
double abs_int(double a, double b);
void print(double value, int size, int space, ofstream &fout);
double rectangles(double a, double b);
double trapezoid(double a, double b);
double Simpson(double rect, double trap);
double Nut_Kot(double a, double b);
double Gauss(double a, double b);
void deltas (double k1, double k2, double abs, ofstream &fout, const string name);

double func(double x){
    return pow(x, 7);
}
double abs_int(double a, double b){
    return (pow(b, 8) - pow(a, 8))/8.0; 
}

void print(double value, int size, int space, ofstream &fout){
    fout << scientific << setprecision(size) << setw(space) << right << value;
}

double rectangles(double a, double b){
    return (b-a) * func((b+a)/2.0);
}

double trapezoid(double a, double b){
    return (func(a) + func(b)) * (b-a)/2.0;
}

double Simpson(double rect, double trap){
    return 2.0*rect/3.0 + trap/3.0;
}

double Nut_Kot(double a, double b){
    int i;
    double f[5] = {0.0};                 // 5 узлов
    double h = (b-a)/4.0;
    for (i = 0; i < 5; i ++)
        f[i] = func(a+i*h);
    return 2.0*h/45.0*(7.0*f[0] + 32.0*f[1] + 12.0*f[2] + 32.0*f[3] + 7.0*f[4]);
}

double Gauss(double a, double b){
    int i;
    double f[3] = {0.0};                 // 3 узла
    double h = (b-a)/2.0;
    double x = sqrt(0.6), c1 = 5.0/9.0, c2 = 8.0/9.0;
    f[0] = func((1-x)*h+a);
    f[1] = func(h+a);
    f[2] = func((1+x)*h+a);
    return h*(c1*f[0] + c2*f[1] + c1*f[2]);
}

void deltas(double k1, double k2, double abs, ofstream &fout, const string name){
    int i;
    double tmp;
    double otn1, otn2;

    otn1 = fabs(abs-k1)/fabs(abs);
    otn2 = fabs(abs-k2)/fabs(abs);
    
    fout << name << "   ";
    print(otn1, 3, 9, fout);
    fout << "   "; 
    print(otn2, 3, 9, fout);
    fout << "\n";
}

int main() {
    int i;
    double x1, x2, h, tmp1, tmp2, abs;
    // K интервалов
    double rect1 = 0.0, trap1 = 0.0, simp1 = 0.0, n_k1 = 0.0, g1 = 0.0;
    h = (b-a)/K;                     
    for (i = 0; i < K; i ++){
        x1 = a+i*h;
        x2 = a+(i+1)*h;
        tmp1 = rectangles(x1, x2); rect1 += tmp1;
        tmp2 = trapezoid(x1, x2); trap1 += tmp2;
        simp1 += Simpson(tmp1, tmp2);
        n_k1 += Nut_Kot(x1, x2);
        g1 += Gauss(x1, x2);
    }
    // 2K интервалов
    double rect2 = 0.0, trap2 = 0.0, simp2 = 0.0, n_k2 = 0.0, g2 = 0.0;
    h = (b-a)/(2*K);                     
    for (i = 0; i < 2*K; i ++){
        x1 = a+i*h;
        x2 = a+(i+1)*h;
        tmp1 = rectangles(x1, x2); rect2 += tmp1;
        tmp2 = trapezoid(x1, x2); trap2 += tmp2;
        simp2 += Simpson(tmp1, tmp2);
        n_k2 += Nut_Kot(x1, x2);
        g2 += Gauss(x1, x2);
    }
    // Погрешности
    ofstream fout ("delta.txt");
    fout << setw(12) << right << "K";
    fout << setw(13) << right << "2K\n";
    abs = abs_int(a, b);
    deltas(rect1, rect2, abs, fout, "rect");
    deltas(trap1, trap2, abs, fout, "trap");
    deltas(simp1, simp2, abs, fout, "simp");
    deltas(n_k1, n_k2, abs, fout, "nu_k");
    deltas(g1, g2, abs, fout, "gaus");

    fout.close();
    return 0;
}
