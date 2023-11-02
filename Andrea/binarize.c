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

//copied quicksort off geeks4geeks cos i dont remember it

int partition(double* u,int low,int high){
	
	double pivot = u[high];
	double temp;

	int i = low-1;

	for(int j=low;j<=high;j++){
		
		if(u[j]<pivot){
		i++;
		//swap
		temp = u[i];
		u[i] = u[j];
		u[j] = temp;
		}
	}
	//swap
	temp = u[i+1];
	u[i+1] = u[high];
	u[high] = temp;
	
	return (i+1);
}


void quickSort(double* u,int low,int high){

	if(low<high)
	{

		int pi=partition(u,low,high);
		
		//recursive
		quickSort(u,low,pi-1);
		quickSort(u,pi+1,high);
	}
}

//sort fun QUICKSORT!!!

double* sort(double* x, int n){
	
	//make deep copy of array
	double* u = malloc(n*sizeof(double));
	
	for(int i = 0; i < n; i++)
		u[i] = x[i];
	
	//sort
	quickSort(u, 0, n-1);
	
	return u;
}


//binarize fun

double binarize(double* x, int n){
	
	//sorted vector
	double* s = sort(x,n);
	
	//vector for indexes
	double* d = malloc((n-1)*sizeof(double));
	
	for(int i = 0; i < n-1; i++)
		d[i] = s[i+1] - s[i];
	
	double t = (s[n-1] - s[0])/(n-1);
	
	double min = s[n-1];
	int index = 0;
	
	for(int i = 0; i < n-1; i++){
		
		if((d[i] > t) && (d[i] < min)){
			min = d[i];
			index = i;
		}
		
		
	}
	free(d);
	double z = s[index + 1];
	return z;
}


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
		u[i] = binarize(matrix[i],cols);
	
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