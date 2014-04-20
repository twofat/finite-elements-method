#include "continue_program_1.h"
#include <string>

int main(){
	int i;
	double * result;
	FILE *output;
    output = fopen("norm.txt","w");
    outestimete[1] = fopen("outestimate20.txt","w");
    outestimete[2] = fopen("outestimate40.txt","w");
    outestimete[3] = fopen("outestimate80.txt","w");
    outestimete[4] = fopen("outestimate160.txt","w");
    output = fopen("norm.txt","w");
    outreal[0] = fopen("outreal.txt","w");
    outreal[1] = fopen("outreal20.txt","w");
    outreal[2] = fopen("outreal40.txt","w");
    outreal[3] = fopen("outreal0.txt","w");
    outreal[4] = fopen("outreal160.txt","w");
/*
    length=1;
    for(global_n=20, global_k=1; global_n<=160; global_n = global_n *2){
		result = problem1();
		fprintf(output, "%.6e\t%d\n",maxnorm(result, global_n, u),global_n);
		fprintf(output, "%.6e\t%d\n",l2norm(result, global_n, u),global_n);
		fprintf(output, "%.6e\t%d\n",h2norm(result, global_n, u, uD),global_n);
		printtofile(result, global_n, u);
		global_k++;
	}



	length = 50;
	for(global_n=20, global_k=1; global_n<=160; global_n = global_n *2){
		result = problem2();
		fprintf(output, "%.6e\t%d\n",maxnorm(result, global_n-1, W1),global_n);
		fprintf(output, "%.6e\t%d\n",l2norm(result, global_n-1, W1),global_n);
		fprintf(output, "%.6e\t%d\n",h2norm(result, global_n-1, W1, W1D),global_n);
		printtofile(result, global_n-1, W1);
		global_k++;
	}


*/
	length = 50;
	for(global_n=20, global_k=1; global_n<=160; global_n = global_n *2){
		result = problem3();
		fprintf(output, "%.6e\t%d\n",maxnorm(result, global_n, W2),global_n);
		fprintf(output, "%.6e\t%d\n",l2norm(result, global_n, W2),global_n);
		fprintf(output, "%.6e\t%d\n",h2norm(result, global_n, W2, W2D),global_n);
		printtofile(result, global_n, W2);
		global_k++;
	}

	return 0;
}








