#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv)
{

	//first argument is rows
	//second argument is columns
	
	
	//unsigned int randval;
	FILE *f;
	time_t t;
	char* p;
	srand((unsigned) time(&t));
	//s = fopen("/dev/random","r");

	/*generating*/

	printf("\nGenerating Random Time Series.\n");

	f = fopen("tseriestest.txt", "a");

	long i = strtol(argv[1], &p, 10);
	long j = strtol(argv[2], &p, 10);
	int x;

	for (long k = 0; k<i; k++)
	{	
		for(long l = 0; l<j; l++){
			x = rand()%100;
			if(l < j-1)
				fprintf(f, "%d,", x);
			else
				fprintf(f, "%d", x);
		}
		if(k < i-1)
			fprintf(f,"\n");
	}

	printf("Generated in tseriestest.txt :)\n");
	fclose(f);
	//fclose(s);

	return 0;

}

