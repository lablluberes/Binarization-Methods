#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>



//count rows in file

int getRows(FILE *f){
	
	//read rows
	int rows = 0;
	char buffer[1024];
	//printf("rows:%d\n",rows);
    while (fgets(buffer, 1024, f)){		
		rows++;	
	} 
	
	
	return rows;
	
}


//count columns in file

int getCols(FILE *f){
	
	int cols = 0;
	
	char line[1024];
	
	fgets(line, 1024, f);
	
	char *scan = line;
	
	//count cols
	
	double dummy;
	int offset = 0;
	
	while(sscanf(scan, "%le,%n", &dummy, &offset) == 1)
	{
		scan += offset;
		cols++;
	}
	
	return cols;
	
	
}


//read file data


double** readFile(FILE *f, int rows, int cols){
	
	double num; 
	double** matrix = NULL;
	
	//allocate rows
	matrix = malloc(rows*sizeof(double));
	
	//allocate columns
	for(int i = 0; i < rows; i++)
		matrix[i] = malloc(cols*sizeof(double));


	//read data
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			fscanf(f, "%le,", &num);
			matrix[i][j] = num;
			//printf("%le ",matrix[i][j]);
		}
		//printf("\n");
	}


	return matrix;
}


//allocate empty result matrix

int** allocMat(int rows, int cols){
	
	int** matrix = NULL;
	
	//allocate rows
	matrix = malloc(rows*sizeof(int));
	
	//allocate columns
	for(int i = 0; i < rows; i++)
		matrix[i] = malloc(cols*sizeof(int));

	
	return matrix;
}


//for double
void freeMat(double** M, int rows){
	
	
	for(int i = 0; i < rows; i++)
		free(M[i]);
	
	free(M);
	
}


//for int
void freeMatInt(int** M, int rows){
	
	
	for(int i = 0; i < rows; i++)
		free(M[i]);
	
	free(M);
	
}



//STEPMINER FUNCTIONS


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

//SSTOT
double getSSTOT(double* x, double xmean, int n){
	
	double sum = 0;
	
	for(int i = 0; i < n; i++)
		sum = sum + (x[i] - xmean) * (x[i] - xmean);
	
	
	return sum;
}

//SSE for onestep

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

//SSE for twostep

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
	//printf("step: %d\n", step);
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



//main


int main(int argc, char **argv)
{

	FILE *f, *w;
	
	f = fopen(argv[1], "r");

	int rows;
	int cols;
	
	
	if (f == NULL){
		printf("fake.\n");
		return 0;
	}
	
	if(argc > 2){
		w = fopen(argv[2], "a");
	}
	else{
		w = fopen("result.txt", "a");
	}
	
	//count rows
	rows = getRows(f);
	
	//reset pointer
	fseek(f, 0, SEEK_SET);
	
	//countCols
	cols = getCols(f);
	printf("size:%dx%d\n",rows,cols);
	
	//reset pointer again
	fseek(f, 0, SEEK_SET);
	
	//read nums
	double** matrix = readFile(f, rows, cols);
	
	//result matrix
	int** u = allocMat(rows,cols);
	
	//divide work
	
	
	double wtime = omp_get_wtime();
	
	
	#pragma omp parallel for
	for(int i = 0; i < rows; i++)
		u[i] = stepminer(matrix[i],cols);
	
	wtime = omp_get_wtime() - wtime;
	printf("time elapsed is %f\n", wtime);
	
	//empty init matrix
	freeMat(matrix,rows);
	
	
	//write result matrix to file
	
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			fprintf(w, " %d", u[i][j]);	
		}
		fprintf(w,"\n");
	}
	
	//empty result matrix
	//freeMatInt(u,rows);

	fclose(f);
	fclose(w);
	
	return 0;
}
