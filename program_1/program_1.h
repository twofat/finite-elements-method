#include<stdio.h>
#include<math.h>
#include<functional>
#include <tr1/functional>
#include <algorithm>
#include <iostream>

#define Inf 10000000000
#define accurate_num 100000
using namespace std;


FILE *outestimete[5];
FILE *outreal[5];
int global_k;
int global_i;
int global_j;
double *global_x;
int global_n;
double *global_test = new double [102];
double length;

void printvector(double *vector, int n){
	int i;
	for(i=1; i<=n; i++){
		printf("%lf ",vector[i]);
	}
	printf("\n");
}
double I(double x){
	if(x<=length&&x>=0){
		return 1;
	}
	return 0;
}

double integrate(double (* k)(double), double (*p)(double), double(*q)(double), double a, double b, int N)
{
	double h = (b-a)/N;
	double i;
	double temp;
	double result=0;
	double t1 = 0.5773502692;
	double t2 = -0.5773502692;
	for(i=0; i<(b-a-h/4); i=i+h){
	    temp = (2*a+2*i+h+h*t1)/2;
	    result += h/2 * (*k)(temp)*(*p)(temp)*(*q)(temp);
	    temp = (2*a+2*i+h+h*t2)/2;
	    result += h/2 * (*k)(temp)*(*p)(temp)*(*q)(temp);
    }
	return result;
}

double* persue(double *a, double *b, double *c, double *f, int n){
	double  *u = new double[n+2];
	double  *v = new double[n+2];
	double  *y = new double[n+2];
	double  *x = new double[n+2];
	int i;

	for(i=1; i<=n; i++)
	{
		u[i] = a[i] - c[i] * v[i-1];
		v[i] = b[i] / u[i];
		y[i] = (f[i] - c[i] * y[i-1]) / u[i];
	}
	for(i=n; i>=1; i--)
	{
		x[i] = y[i] - v[i] * x[i+1];
	}
	delete[] u;
	delete[] v;
	delete[] y;
	return x;
}

void init(int n){
	int i;
	double h = length/(double)n;
	global_x = new double[n+2];
	for(i=0; i<=n+1; i++){
		global_x[i] = h * (double)i;
	}
}

double phii(double x){
	double h = length/(double)global_n;
	if(x <= global_x[global_i] && x>global_x[global_i-1]){
		return (x-global_x[global_i-1])/h;
	}
	if(x < global_x[global_i+1] && x>global_x[global_i]){
		return (global_x[global_i+1]-x)/h;
	}
	return 0;
}

double phij(double x){
	double h = length/(double)global_n;
	if(x <= global_x[global_j] && x>global_x[global_j-1]){
		return (x-global_x[global_j-1])/h;
	}
	if(x < global_x[global_j+1] && x>global_x[global_j]){
		return (global_x[global_j+1]-x)/h;
	}
	return 0;
}

double dphii(double x){
	double h = length/(double)global_n;
	if(x <= global_x[global_i] && x>global_x[global_i-1]){
		return 1/h;
	}
	if(x < global_x[global_i+1] && x>global_x[global_i]){
		return -1/h;
	}
	return 0;
}
double dphij(double x){
	double h = length/(double)global_n;
	if(x <= global_x[global_j] && x>global_x[global_j-1]){
		return 1/h;
	}
	if(x < global_x[global_j+1] && x>global_x[global_j]){
		return -1/h;
	}
	return 0;
}
double k(double x){
	if(x>1){
		return 0;
	}
	if(x<0.5){
		return 1;
	}
	if(x>=0.5){
		return 0.5;
	}
}
double u(double x){
	double temp;
	if(x<0.5){
		temp = x + 0.5*(1-exp(x))/sqrt(exp(1));
		return temp;
	}
	if(x>=0.5){
		temp = x + 0.5*(1-exp(0.5))/sqrt(exp(1));
		return temp;
	}
}

double uD(double x){
	double temp;
	if(x<0.5){
		temp = 1 + 0.5*(-exp(x))/sqrt(exp(1));
		return temp;
	}
	if(x>=0.5){
		return 1;
	}
}


double estimate(double x, double* result, int n){
	double temp=0;
	for(global_i=1; global_i<=n; global_i++){
		temp += phii(x) * result[global_i];
	}
	return temp;
}

double* problem1(){
	int n = global_n;
	int i;
	double max = -1;
	double *a = new double [n+2];
	double *b = new double [n+2];
	double *c = new double [n+2];
	double *f = new double [n+2];
	double *result = new double[n+2];
	double zk;
	length = 1;
	init(global_n);

	for(global_i = 1; global_i <= n; global_i++){
		a[global_i] = integrate(k, dphii, dphii, global_x[global_i-1], global_x[global_i+1], accurate_num)
			+ integrate(I, dphii, phii, global_x[global_i-1], global_x[global_i+1], accurate_num) ;
	}
	for(global_i = 1; global_i <= n-1; global_i++){
		global_j = global_i + 1;
		b[global_i] = integrate(k, dphij, dphii, global_x[global_i], global_x[global_j], accurate_num)
			+ integrate(I, dphij, phii, global_x[global_i], global_x[global_j], accurate_num) ;
	}
	for(global_i = 2; global_i <= n; global_i++){
		global_j = global_i - 1;
		c[global_i] = integrate(k, dphij, dphii, global_x[global_j], global_x[global_i], accurate_num)
			+ integrate(I, dphij, phii, global_x[global_j], global_x[global_i], accurate_num) ;
	}
	for(global_i = 1; global_i <= n; global_i++){
		f[global_i] = integrate(I, I, phii, global_x[global_i-1], global_x[global_i+1], accurate_num);
	}
	f[n] += 0.5;
	result = persue(a, b, c, f, n);
	return result;
}





double F(double x){
	double q = 200;
	double F = 100;
	return q/F/2*x*(length-x);
}

double k2(double x){
	if(x>length){
		return 0;
	}
	return 8.8e5;
}
double W1(double x){
	double a = 100 * length * length / (8.8e7);
	double b = 200 * length * length * length * length / 2 / (8.8e7);
	double t = x / length;
	double result = b/a *(-t * t + t - 2/a + 2/a/sinh(sqrt(a))*(sinh(sqrt(a)*t) + sinh(sqrt(a)*(1-t)))  );
	return result;
}
double W1D(double x){
	double a = 100 * length * length / (8.8e7);
	double b = 200 * length * length * length * length / 2 / (8.8e7);
	double t = x / length;
	double result = b/a/length *(-2 * t + 1 + 2/sqrt(a)/sinh(sqrt(a))*(cosh(sqrt(a)*t) - cosh(sqrt(a)*(1-t)))  );
	return result;
}
double W2(double x){
	double a = 100 * length * length / (8.8e7);
	double b = 200 * length * length * length * length / 2 / (8.8e7);
	double t = x / length;
	double A = -b/a;
	double B = b/a;
	double C = -2*b/a/a;
	double D = b/a/sqrt(a)/cosh(sqrt(a)) + 2*b/a/a/cosh(sqrt(a))/sinh(sqrt(a));
	double E = 2*b/a/a/sinh(sqrt(a));
	double result = A * t*t + B*t + C + D * sinh(sqrt(a)*t) + E*sinh(sqrt(a)*(1-t));
	return result;
}
double W2D(double x){
	double a = 100 * length * length / (8.8e7);
	double b = 200 * length * length * length * length / 2 / (8.8e7);
	double t = x / length;
	double A = -b/a;
	double B = b/a;
	double C = -2*b/a/a;
	double D = b/a/sqrt(a)/cosh(sqrt(a)) + 2*b/a/a/cosh(sqrt(a))/sinh(sqrt(a));
	double E = 2*b/a/a/sinh(sqrt(a));
	double result = 2* A * t + B + sqrt(a) * D * cosh(sqrt(a)*t) - sqrt(a)*E*cosh(sqrt(a)*(1-t));
	return result/length;
}

double * problem2(){
	int n = global_n;
	int i;
	double temp=0;
	double max = -1;
	double *a = new double [n+2];
	double *b = new double [n+2];
	double *c = new double [n+2];
	double *f = new double [n+2];
	double *result = new double[n+2];
	double zk;
	length = 50;
	init(global_n);

	for(global_i = 1; global_i <= n-1; global_i++){
		a[global_i] = integrate(k2, dphii, dphii, global_x[global_i-1], global_x[global_i+1], accurate_num)
			+ integrate(I, phii, phii, global_x[global_i-1], global_x[global_i+1], accurate_num) ;
	}
	for(global_i = 1; global_i <= n-2; global_i++){
		global_j = global_i + 1;
		b[global_i] = integrate(k2, dphij, dphii, global_x[global_i], global_x[global_j], accurate_num)
			+ integrate(I, phij, phii, global_x[global_i], global_x[global_j], accurate_num) ;
}
	for(global_i = 2; global_i <= n-1; global_i++){
		global_j = global_i - 1;
		c[global_i] = integrate(k2, dphij, dphii, global_x[global_j], global_x[global_i], accurate_num)
			+integrate(I, phij, phii, global_x[global_j], global_x[global_i], accurate_num) ;
}
	for(global_i = 1; global_i <= n-1; global_i++){
		f[global_i] = integrate(I, F, phii, global_x[global_i-1], global_x[global_i+1], accurate_num);
	}
//	printvector(a, n-1);
//	printvector(b, n-1);
//	printvector(c, n-1);
//	printvector(f, n-1);
	result = persue(a, b, c, f, n-1);
	return result;
}



double* problem3(){
	int n = global_n;
	int i;
	double temp=0;
	double max = -1;
	double *a = new double [n+2];
	double *b = new double [n+2];
	double *c = new double [n+2];
	double *f = new double [n+2];
	double *result = new double[n+2];
	length = 50;
	init(global_n);

	for(global_i = 1; global_i <= n; global_i++){
		a[global_i] = integrate(k2, dphii, dphii, global_x[global_i-1], global_x[global_i+1], accurate_num)
			+ integrate(I, phii, phii, global_x[global_i-1], global_x[global_i+1], accurate_num) ;
	}
	for(global_i = 1; global_i <= n-1; global_i++){
		global_j = global_i + 1;
		b[global_i] = integrate(k2, dphij, dphii, global_x[global_i], global_x[global_j], accurate_num)
			+ integrate(I, phij, phii, global_x[global_i], global_x[global_j], accurate_num) ;
	}
	for(global_i = 2; global_i <= n; global_i++){
		global_j = global_i - 1;
		c[global_i] = integrate(k2, dphij, dphii, global_x[global_j], global_x[global_i], accurate_num)
			+integrate(I, phij, phii, global_x[global_j], global_x[global_i], accurate_num) ;
	}
	for(global_i = 1; global_i <= n; global_i++){
		f[global_i] = integrate(I, F, phii, global_x[global_i-1], global_x[global_i+1], accurate_num);
	}
	result = persue(a, b, c, f, n);
	return result;
}












































