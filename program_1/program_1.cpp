#include "program_1.h"


int main(){
	int i;
/*	for(i=0; i<=100; i++){
		global_test[i] = 0.01 * (double)i;
	}

	for(global_n=20; global_n<=160; global_n = global_n *2){
		problem1();
	}
	
*/
	for(i=0; i<=100; i++){
		global_test[i] = 0.01 * (double)i * 50;
	}
	for(global_n=20; global_n<=20; global_n = global_n *2){
		problem3();
	}
	return 0;
}






