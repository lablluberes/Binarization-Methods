#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

//HOW TO USE
//programname (name of matrix file) (OPTIONAL: name of result file, defaults to result.txt)


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


//for double
void freeMat(double** M, int rows){
	
	
	for(int i = 0; i < rows; i++)
		free(M[i]);
	
	free(M);
	
}

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

//STEPMINER

double onestep(double* x, int n){
	
	//distributing threads
	//int nt = omp_get_num_threads();
	//int chunkSize = n/nt;
	
	double xmean = mean(x, 0, n);
	double sstot = getSSTOT(x, xmean, n);
	
	double min = sstot;
	
	//stuff i need to calculate
	double left, right, sse, t;
	
	
	//onestep

	for(int i = 0; i < n-1; i++){
		
		left = mean(x, 0, i+1);
		
		right = mean(x, i+1, n);
		
		sse = onestepSSE(x, left, right, i, n);
		
		if(sse < min){
			
			min = sse;
			
			t = (right + left)/2;
		}
		
	}
	return t;
	
}



//main


int main(int argc, char **argv)
{

	FILE *f, *w;
	
	f = fopen(argv[1], "r");

	int rows;
	int cols;
	
	
	if (f == NULL){
		printf("no input.\n");
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
	
	//result array of thresholds
	double* u = malloc(rows*sizeof(double));
	
	//divide work
	
	
	double wtime = omp_get_wtime();
	
	
	#pragma omp parallel for
	for(int i = 0; i < rows; i++)
		u[i] = onestep(matrix[i],cols);
	
	wtime = omp_get_wtime() - wtime;
	printf("time elapsed is %f\n", wtime);
	
	//empty init matrix
	freeMat(matrix,rows);
	
	
	//write result matrix to file
	
	for(int i = 0; i < rows; i++){
		fprintf(w, "%le\n", u[i]);
	}
	

	fclose(f);
	fclose(w);
	
	return 0;
}
