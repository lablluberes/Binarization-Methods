#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

//Get mean of values in array
double mean(double* array, int i, int j){
	
	double sum = 0;
	//printf("i = %d, j = %d\n", i, j);
	
	for(int k = i; k < j; k++){
		sum = sum + array[k];	
	}
	
	double m = sum/(j-i);
	return m;
	
}


//Count the amount of values in the set
int countFile(FILE *f){
	
	int counter = 0;
	double num;
	while(fscanf(f, "%le,", &num)==1){
		counter++;
	}
	return counter;
}

//Read numbers and store them
//Please reset file reading pointer before
//calling this
double* readFile(FILE *f, int n){
	
	double* data = malloc(n*sizeof(double));
	double num;
	
	for(int i = 0; i < n; i++){
		fscanf(f, "%le, ", &num);
		data[i] = num;
	}
	
	return data;
	
}


//SSTOT
double getSSTOT(double* x, double xmean, int n){
	
	double sum = 0;
	
	for(int i = 0; i < n; i++)
		sum = sum + (x[i] - xmean) * (x[i] - xmean);
	
	
	return sum;
}


double onestepSSE(double* x, double left, double right, int i, int n){
	
	double sum = 0;
	
	for(int k = 0; k < n; k++){
		if(k < i+1)
			sum = sum + (x[k] - left) * (x[k] - left);
		else
			sum = sum + (x[k] - right) * (x[k] - right);
		
	}
	
	return sum;
	
}


double twostepSSE(double* x, double left, double right,int i, int j, int n){
	
	double sum = 0;
	//printf("in twostep\n");
	
	for(int k = 0; k < n; k++){
		if(k < i+1 || k > j)
			sum = sum + (x[k] - left) * (x[k] - left);
		else
			sum = sum + (x[k] - right) * (x[k] - right);
		
	}
	
	return sum;
	
}

//STEPMINER

int* stepminer(double* x, int n){
	
	//distributing threads
	//int nt = omp_get_num_threads();
	//int chunkSize = n/nt;
	
	double xmean = mean(x, 0, n);
	double sstot = getSSTOT(x, xmean, n);
	
	double min1, min2;
	min1 = sstot;
	min2 = sstot;
	
	double maxF = 0;
	
	//stuff i need to calculate
	double left, left1, right, sse, ssr, F, msr, mse;
	int m, step, greaterThan;
	
	//stuff i need to store indexes
	int i1, i2, j2;
	
	double wtime = omp_get_wtime();
	//parallel them separately
	//onestep

	for(int i = 0; i < n-1; i++){
		
		left1 = mean(x, 0, i+1);
		
		for(int j = i+1; j < n-1; j++){
				

			right = mean(x, i+1, j+1);

			left = ((left1 * (i+1)) + (mean(x, j+1, n) * (n-j-1)))/(n-j+i);

			sse = twostepSSE(x, left, right, i, j, n);
					
					
			if(sse < min2){
				min2 = sse;
				ssr = sstot - sse;
				if(n > 4)
				m = 4;
				else
				m = 2;
			
			
				msr = ssr/(m-1);
				mse = sse/(n-m);
				F = msr/mse;
					
				
					if(F > maxF){
						maxF = F;
						i2 = i;
						j2 = j;
						greaterThan = (left > right) ? 1 : 0;
						step = 0;
					}
				
			}
		}
		
		right = mean(x, i+1, n);
		
		sse = onestepSSE(x, left1, right, i, n);
		
		if(sse < min1){
			
			min1 = sse;
			
			ssr = sstot - sse;
			if(n > 4)
			m = 3;
			else
			m = 2;
		
			msr = ssr/(m-1);
			mse = sse/(n-m);
			F = msr/mse;
			
			
				if(F > maxF){
					maxF = F;
					i1 = i;
					greaterThan = (left > right) ? 1 : 0;
					step = 1;
				}
			
		}
		
	}
	printf("step: %d\n", step);
	wtime = omp_get_wtime() - wtime;
	printf("time elapsed is %f\n", wtime);
	//printf("twostep done\n");
	
	int* u = malloc(n*sizeof(int));
	
	if(step == 1){
		for(int k = 0; k < n; k++){
			if(k < i1+1)
				u[k] = (greaterThan == 1) ? 1 : 0;
			else
				u[k] = (greaterThan == 1) ? 0 : 1;
		}	
	}
		
	else{
		for(int k = 0; k < n; k++){
			if(k < i2+1 || k >= j2+1)
				u[k] = (greaterThan == 1) ? 1 : 0;
			else
				u[k] = (greaterThan == 1) ? 0 : 1;
		}	
	}

	
	//return final vector later
	return u;
	
}

int main(int argc, char **argv){


	FILE *f; //*w;
	
	//open file to read
	
	f = fopen(argv[1], "r");
	//w = fopen("stepminermodc.txt", "a");

	int n = countFile(f);
	printf("size: %d\n", n);
	//for testing
	//int n = 6;
	double* x;
	int* u;

	//reset pointer to start of file
	fseek(f, 0, SEEK_SET);
	
	//for(int a = 0; a < 418; a++){
		
		x = readFile(f, n);
		
		u = stepminer(x, n);
		
		for(int i = 0; i < n; i++)
		{
			//fprintf(w, " %d", u[i]);	
			printf("%d ", u[i]);
		}
		
		//fprintf(w, "\n");	
		printf("\n");
	//}
	fclose(f);
	
	return 0;

}