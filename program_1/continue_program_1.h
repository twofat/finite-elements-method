#ifndef CONTINUE_PROGRAM_1_H_INCLUDED
#define CONTINUE_PROGRAM_1_H_INCLUDED
#include "program_1.h"

double estimateD(double x, double* result, int n){
	double temp=0;
	for(global_i=1; global_i<=n; global_i++){
		temp += dphii(x) * result[global_i];
	}
	return temp;
}

double maxnorm(double *result, int n, double(*u)(double)){
    double temp;
    int i;
	double max = -1;
    for(i=0; i<=100; i++){
        global_test[i] = 0.01 * (double)i * length;
    }
    for(i=0; i<=global_n; i++){
//		printf("estimate: %lf \t real: %lf\n", estimate(global_test[i], result, n), u(global_test[i]));
		temp = fabs(u(global_x[i]) - estimate(global_x[i], result, n));
		if(temp>max){
			max = temp;
		}
	}
	return max;
}
double l2norm(double *result, int n, double(*u)(double)){
    double i;
    double a=0;
    double b=length;
	double temp;
	double tempvalue;
	double re=0;
	double t1 = 0.5773502692;
	double t2 = -0.5773502692;
	double h = (b-a)/accurate_num;
	for(i=0; i<(b-a-h/4); i=i+h){
	    temp = (2*a+2*i+h+h*t1)/2;
	    tempvalue = (u(temp) - estimate(temp, result, n)) * (u(temp) - estimate(temp, result, n));;
	    re += h/2 * tempvalue;
	    temp = (2*a+2*i+h+h*t2)/2;
	    re += h/2 * tempvalue;
    }
    return sqrt(re);
}

double h2norm(double *result, int n, double(*u)(double), double (*uD)(double)){
    double i;
    double a=0;
    double b=length;
	double temp;
	double tempvalue;
	double re1 = 0;
	double re2 = 0;
    double re;
	double t1 = 0.5773502692;
	double t2 = -0.5773502692;
	double h = (b-a)/accurate_num/5;
	for(i=0; i<(b-a-h/4); i=i+h){
	    temp = (2*a+2*i+h+h*t1)/2;
	    tempvalue = (u(temp) - estimate(temp, result, n)) * (u(temp) - estimate(temp, result, n));
	    re1 += h/2 * tempvalue;
	    temp = (2*a+2*i+h+h*t2)/2;
	    tempvalue = (u(temp) - estimate(temp, result, n)) * (u(temp) - estimate(temp, result, n));
	    re1 += h/2 * tempvalue;
    }
    for(i=0; i<(b-a-h/4); i=i+h){
	    temp = (2*a+2*i+h+h*t1)/2;
	    tempvalue = (uD(temp) - estimateD(temp, result, n)) * (uD(temp) - estimateD(temp, result, n));
	    re2 += h/2 * tempvalue;
	    temp = (2*a+2*i+h+h*t2)/2;
	    tempvalue = (uD(temp) - estimateD(temp, result, n)) * (uD(temp) - estimateD(temp, result, n));
	    re2 += h/2 * tempvalue;
    }
    re = re1 + re2;
    return sqrt(re);
}
void printtofile(double *result, int n, double (*u)(double)){
    int i;
    double h = length/1000;
    fprintf(outestimete[global_k], "estimate\n");
    for(i=0;i<=global_n;i++){
        fprintf(outestimete[global_k], "%16.9lf\n", estimate(global_x[i], result, n));
        fprintf(outreal[global_k], "%16.9lf\n", u(global_x[i]));
    }
    fprintf(outreal[global_k], "real\n");
    if(global_k>1) return;
    for(i=0;i<=1000;i++){
        fprintf(outreal[0], "%16.9lf\n", u(h*i));
    }
}







#endif // CONTINUE_PROGRAM_1_H_INCLUDED
