#include "Headers.h"


//if after all the iterasions in checkForGoodClusters we dont get good enaf qualety  then we print this
void writeToFileNotFound() {
	FILE* f = fopen(FILE_OUT, "w");
	fprintf(f, "\n\nNever Found Good Clusters\n\n");
}

//checking if file opend corectly
void checkFile(const FILE* f)
{
	if (f == NULL)
	{
		puts("Could not open the file , please check file path");
		MPI_Finalize();
		exit(1);
	}
}

//initielise file by data from define 
void initFile(char* fileName, int numofPoints)
{
	FILE* f = fopen(INIT_FILENAME, "w");
	checkFile(f);
	/*First line of file*/
	fprintf(f, "%d %d %lf %lf %d %lf\n", NUMPOINTS, NUMCLUSTERS, T, DT, LIMIT, QUALITY_MEASURE);

	for (int i = 0; i < numofPoints; i++) {
		/* Write to file Xi , Yi , Zi*/
		fprintf(f, "%lf %lf %lf ", generatRandomPositiveNegitiveValue(MAXVALUE), generatRandomPositiveNegitiveValue(MAXVALUE), generatRandomPositiveNegitiveValue(MAXVALUE));
		/* Write to file Vxi , Vyi , Vzi*/
		fprintf(f, "%lf %lf %lf\n", generatRandomPositiveNegitiveValue(MAXVALUE/100), generatRandomPositiveNegitiveValue(MAXVALUE / 100), generatRandomPositiveNegitiveValue(MAXVALUE / 100));
	}
	fclose(f);
}


//reading point from file
Point* readPointsFromFile(const char* fileName, int* numofPoints, int* num_clusters, double* t, double* dt, int* limit, double* quality_measure)
{

	FILE* f = fopen(fileName, "r");
	checkFile(f);

	fscanf(f, "%d %d %lf %lf %d %lf", numofPoints, num_clusters, t, dt, limit, quality_measure);
	Point* allpoints = (Point*)calloc(*numofPoints, sizeof(Point));
	checkDynamicAllocation(allpoints);

	for (int i = 0; i < *numofPoints; i++)
	{
		/* Read from file Xi , Yi , Zi to allpoints[i]*/
		fscanf(f, "%lf %lf %lf", &allpoints[i].x, &allpoints[i].y, &allpoints[i].z);
		/* Read from file Vxi , Vyi , Vzi to allpoints[i]*/
		fscanf(f, "%lf %lf %lf", &allpoints[i].Vx, &allpoints[i].Vy, &allpoints[i].Vz);
	}
	fclose(f);
	return allpoints;
}

void writeToFile(int numofClusters, Cluster *clusterArray, double elapsedTime, double quality)
{
	int i;
	FILE *fp = fopen(FILE_OUT, "w");
	if (fp == NULL) {
		printf("the file is not open\n");
		exit(1);
	}
	fprintf(fp, "First occurrence at t = %lf with q = %lf\nCenters of the clusters:\n", elapsedTime, quality);
	for (i = 0; i < numofClusters; i++)
	{
		fprintf(fp, " Cluster Centre ID : [ %d ]  ( %lf , %lf,  %lf )\n", clusterArray[i].id, clusterArray[i].x, clusterArray[i].y, clusterArray[i].z);
	}
	fclose(fp);
}
