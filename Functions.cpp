#include "Headers.h"
#include "StructsHeader.h"



//generates random nambers for init file
double generatRandomPositiveNegitiveValue(int value)
{
	double ii = (-rand() % value + rand() % value);
	return ii;
}


void freeAll(Point* allpoints, Cluster* clustersArray, Point* pointsBufferForEachProcess) {
	free(allpoints);
	free(clustersArray);
	free(pointsBufferForEachProcess);
}


void checkDynamicAllocation(const void* ptr)
{
	if (!ptr)
	{
		printf("Dynamic Memory Allocation Error --> Not Enough Memory !");
		fflush(stdout);
		MPI_Finalize();
		exit(2);
	}
}

// Master Process Broadcasts to all processes : numofpoints,numclusters,t,dt,limit
void broadcastData(int* numofpoints, int* numclusters, double* t, double *dt, int *limit, double* qualitymeasure)
{
	MPI_Bcast(numofpoints, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(numclusters, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(t, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(limit, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(dt, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(qualitymeasure, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

}

//scaters points from arrPoints in master to PointsPerProcess in every proccess 
Point* scatterPoints(int sizeOfAllPoints, Point *arrPoints, int* numofPointsPerRank, Point *PointsPerProcess, int myRank, int numOfProcs, MPI_Datatype MPI_POINT_TYPE)
{
	int i, wholesize, remaindersize;
	MPI_Status status;

	wholesize = sizeOfAllPoints / numOfProcs;
	remaindersize = sizeOfAllPoints % numOfProcs;

	//master will have the remainder of points
	if (myRank == MASTER)
		*numofPointsPerRank = wholesize + remaindersize;
	else
		*numofPointsPerRank = wholesize;

	PointsPerProcess = (Point*)malloc(sizeof(Point)*(*numofPointsPerRank));
	checkDynamicAllocation(PointsPerProcess);
	if (myRank == MASTER)
	{
		//  Master sends to himself his share of points in Parallel in openMP until wholesize + remaindersize.
#pragma omp parallel for
		for (i = 0; i < *numofPointsPerRank; i++)
		{
			(PointsPerProcess)[i] = arrPoints[i];
		}
		//  Master sends to all other processes their share of points 
		for (i = 1; i < numOfProcs; i++)
		{
			MPI_Send(arrPoints + (*numofPointsPerRank + (wholesize*(i - 1))), wholesize, MPI_POINT_TYPE, i, MPI_FLAG, MPI_COMM_WORLD);
		}
	}
	else {
		// Other Processes recieve their share of points from Master.
		MPI_Recv(PointsPerProcess, wholesize, MPI_POINT_TYPE, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}
	return PointsPerProcess;
}


Cluster* initClusters(int numOfClusters, Point *pointsArray, Cluster* allClusters, MPI_Datatype MPI_CLUSTER_TYPE, int myRank)
{
	int i;
	//Master initializes all clusters' datafields. First K (numOfClusters) Points in Order of PointsArray will be assigned as initial Clusters.
	if (myRank == MASTER)
	{
#pragma omp parallel for
		for (i = 0; i < numOfClusters; i++)
		{
			allClusters[i].x = pointsArray[i].x;
			allClusters[i].y = pointsArray[i].y;
			allClusters[i].z = pointsArray[i].z;
			allClusters[i].diameter = 0;
			allClusters[i].id = i;
			allClusters[i].numPointsInCluster = 0;
		}
	}
	//Master broadcasts all first K Initial Points- Clusters' datafields to all salves. (Processes 1..n-1)
	MPI_Bcast(allClusters, numOfClusters, MPI_CLUSTER_TYPE, MASTER, MPI_COMM_WORLD);
	return allClusters;
}

void printArray(Point* allPoints, int size, int processnumber)
{
	printf("\n\n\n");
	fflush(stdout);
	for (int i = 0; i < size; i++)
	{
		printf("processnumber :  %d  allPoints[%d].x = %lf    allPoints[%d].y = %lf     allPoints[%d].z = %lf \n\n", processnumber, i, allPoints[i].x, i, allPoints[i].y, i, allPoints[i].z);
		fflush(stdout);
		printf("processnumber :  %d  allPoints[%d].Vx = %lf     allPoints[%d].Vy =  %lf      allPoints[%d].Vz =  %lf     ----->  This Point's membership to Cluster ID [%d]  \n\n\n\n", processnumber, i, allPoints[i].Vx, i, allPoints[i].Vy, i, allPoints[i].Vz, allPoints[i].idCluster);
		fflush(stdout);
	}
	printf("\n\n\n");
	fflush(stdout);
}


//in checkForGoodClusters we will gather clusters and check thear qualety
Cluster* checkForGoodClusters(double quality_measure, int* numofPoints, int* num_clusters, double t, double dt, int limit, Point* allpoints, int* numofPointsPerProcess, Point*  pointsPerProcess, int* myid, int* numprocs, MPI_Datatype  MPI_POINT_TYPE, MPI_Datatype MPI_CLUSTER_TYPE, Cluster* allClusters)
{
	double quality = 0;
	double elapsedTime = 0;
	int iterationNumber = 1;

	//looping T/dt times or until we fing good enuff qualety
	for (double k = 0; k <= t; k += dt)
	{
		//create clusters
		mpiKmeans(pointsPerProcess, *numofPointsPerProcess, *num_clusters, allClusters, limit, myid);

		//colecting  new points in to master, for qualety check
		collectPoints(*myid, allpoints, *numofPoints, *numofPointsPerProcess, pointsPerProcess, MPI_POINT_TYPE, *numprocs);

		//measher time for printing to file the time it took to find good clasters
		elapsedTime = iterationNumber*dt;

		if (*myid == MASTER)
		{
			calculateQuality(*numofPoints, *num_clusters, allpoints, allClusters, &quality);
			if (quality <= quality_measure)
			{
				printGoodClusters(quality, elapsedTime, allClusters, *num_clusters, *myid);
				writeToFile(*num_clusters, allClusters, elapsedTime, quality);
			}
			printf("\n\n------------------quality = %lf------------------------------------\n\n", quality);
			fflush(stdout);
		}

		// notify all slaves if we reached to requested quality
		MPI_Bcast(&quality, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (quality <= quality_measure)
			return allClusters;

		//moving points for anather round
		pointsLocation(*numofPointsPerProcess, dt, pointsPerProcess);
		//	printArray(pointsPerProcess, *numofPointsPerProcess, *myid);
		iterationNumber++;
	}
	//if  i am out of the for loop it means i didn't find good clasters
	writeToFileNotFound();
	return allClusters;
}

//calculate qualety of clusters by one proccess
void calculateQuality(int numberOfPoints, int numberOfClusters, Point* pointsArray, Cluster* clusterArray, double *quality)
{
	int i, j;
	double distance = 0;
	double theQual = 0;
	int denominator = numberOfClusters*(numberOfClusters - 1);

	giveClustersDiameter(numberOfPoints, numberOfClusters, pointsArray, clusterArray);

	//private(distance,j)

	// calculate the qualety from array of clusters ond they diameters
#pragma omp parallel for private(distance,j) reduction(+ : theQual)
	for (i = 0; i < numberOfClusters; i++)
	{
		for (j = 0; j < numberOfClusters; j++)
		{
			if (i != j)
			{
				distance = distanceClusters(&clusterArray[i], &clusterArray[j]);
				theQual += (clusterArray[i].diameter / distance);
			}
		}
	}
	theQual = theQual / denominator;
	//returns qualety
	*quality = theQual;
}

//calculate diameters of clusters with cuda help
void giveClustersDiameter(int numberOfPoints, int numOfClusters, Point* pointsArray, Cluster* clusterArray)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < numOfClusters; i++)
	{
		clusterArray[i].diameter = theDiameter(numberOfPoints, numOfClusters, pointsArray, clusterArray[i].id);
	}
}

double theDiameter(int numberOfPoints, int numOfClusters, Point* pointsArray, int clusterId)
{
	int i, j;
	double distance;
	double TheMaxDistance = 0;
	for (i = 0; i < numberOfPoints; i++)
	{
		if (pointsArray[i].idCluster == clusterId)
		{
			for (j = (i + 1); j < numberOfPoints; j++)
			{
				if (pointsArray[j].idCluster == clusterId)
				{
					distance = distancePoints(&pointsArray[i], &pointsArray[j]);
					if (distance > TheMaxDistance)
						TheMaxDistance = distance;
				}
			}
		}
	}
	return sqrt(TheMaxDistance);
}

// calculate distance between to points without sqrt
double distancePoints(Point* point1, Point* point2)
{
	double dist = 0;
	dist = pow((point1->x - point2->x), 2.0) + pow((point1->y - point2->y), 2.0) + pow((point1->z - point2->z), 2.0);
	return dist;
}
// calculate distance between to clusters
double distanceClusters(Cluster* c1, Cluster* c2)
{
	double dist = 0;
	dist = sqrt(pow((c1->x - c2->x), 2.0) + pow((c1->y - c2->y), 2.0) + pow((c1->z - c2->z), 2.0));
	return dist;
}

/* Returns the distance between two points - in Square */
double euclid_dist_2(Point p, Cluster c) {
	double euclid_dist_square = 0.0;
	euclid_dist_square = pow((p.x - c.x), 2.0) + pow((p.y - c.y), 2.0) + pow((p.z - c.z), 2.0);
	return euclid_dist_square;
}


/* Returns index of nearest cluster centre point to the the given point parameter.*/
int find_nearest_cluster(int numClusters, Point point, Cluster* clusters, int *myid)
{

	int   index = 0, i;
	float dist = 0, min_dist = 0;

	/* find the cluster id that has min distance to object */
	index = clusters[0].id;
	/*Calculate the minimum distance between single point and each Cluster's Middle-Centre Point*/
	min_dist = euclid_dist_2(point, clusters[0]);

	for (i = 1; i < numClusters; i++) {
		dist = euclid_dist_2(point, clusters[i]);
		/* no need square root */
		if (dist < min_dist) { /* find the min and its array index */
			min_dist = dist;
			index = clusters[i].id;
		}
	}

	return(index);
}

void mpiKmeans(Point* pointsPerProcess, int num_of_points, int numClusters, Cluster* clustersArr, int limit, int* myid)
{
	int i, loop = 0, nearestClusterID;
	int continuekflag = TRUE;

	while (continuekflag && loop < limit)
	{

		/* For subset of points in each process - Iterate over this subset of points. This is the Parallel Job, each process does it parallel*/
		continuekflag = FALSE;

#pragma omp parallel for
		for (i = 0; i < numClusters; i++) {
			clustersArr[i].numPointsInCluster = 0;
		}

		for (i = 0; i < num_of_points; i++) {
			//gives the id of nearest claster center to the point
			nearestClusterID = find_nearest_cluster(numClusters, pointsPerProcess[i], clustersArr, myid);


			if (pointsPerProcess[i].idCluster != nearestClusterID)
			{
				pointsPerProcess[i].idCluster = nearestClusterID;
				continuekflag = TRUE;
			}
			clustersArr[nearestClusterID].numPointsInCluster++;
		}
		//collecting point's centers by clasters
		collectPointsCenters(pointsPerProcess, clustersArr, num_of_points, numClusters);

		// calculate new claster centers for enather k-means iteration
		calculateAverage(numClusters, clustersArr);

		toContinueKMeans(&continuekflag);

		loop++;
	}
}

/* MPI_Allreduce will reduce the values and distribute the results to all processes. */
void calculateAverage(int numClusters, Cluster* arrClusters)
{
	int i, total = 0;
	double x = 0, y = 0, z = 0;

	for (i = 0; i < numClusters; i++)
	{

		MPI_Allreduce(&(arrClusters[i].x), &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		MPI_Allreduce(&(arrClusters[i].y), &y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		MPI_Allreduce(&(arrClusters[i].z), &z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		MPI_Allreduce(&(arrClusters[i].numPointsInCluster), &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		/* Each Process Calculate new Clusters' Centre of Gravity. Now all Processes will contain updated Cluster Centres.*/
		arrClusters[i].x = x / total;
		arrClusters[i].y = y / total;
		arrClusters[i].z = z / total;
	}
}

//Collect all Data to Master. Now Master has all points & Cluster
void collectPoints(int myRank, Point *allPoints, int numOfAllPoints, int numPointsPerProcess, Point *PointsForProcess, MPI_Datatype MPI_POINT_TYPE, int numOfProcs)
{
	int i;
	int div = numOfAllPoints / numOfProcs;
	if (myRank == MASTER)
	{

#pragma omp parallel for 
		for (i = 0; i < numPointsPerProcess; i++)
		{
			allPoints[i] = PointsForProcess[i];
		}
		MPI_Status status;

		/* Master recieves points from each other process.*/
#pragma omp parallel for 
		for (i = 1; i < numOfProcs; i++)
		{
			MPI_Recv(allPoints + (numPointsPerProcess + (div*(i - 1))), div, MPI_POINT_TYPE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}
	}
	/* Other processes send their share of points to master.*/
	else {
		MPI_Send(PointsForProcess, numPointsPerProcess, MPI_POINT_TYPE, MASTER, 0, MPI_COMM_WORLD);
	}

}

/*Input continuekflag.*/
void toContinueKMeans(int* continuekflag)
{
	int tempflag;
	// Send Sum of continueflags to all processes
	MPI_Allreduce(continuekflag, &tempflag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	/*If any of the processes' continueflag is 1 - then sum >=1. If any of the processes' continueflag = 1 --> we need to continue mpi_kmeans. */
	if (tempflag >= 1)
		*continuekflag = TRUE;
	/*Else all processes' continueflag are 0 - then sum < 1. Then we know for sure that we need to stop for each process its mpi_kmeans operation. */
	else
		*continuekflag = FALSE;
}

void printGoodClusters(double quality_found, double time, Cluster* clustersArray, int num_clusters, int id)
{
	printf("\n First Occurence t = %lf with q = %lf \n\n ", time, quality_found);
	fflush(stdout);
	printf("Centers of the clusters: \n \n");
	fflush(stdout);
	printClusters(clustersArray, num_clusters, id);
}

void printClusters(Cluster* clustersArray, int num_clusters, int myid)
{
	printf("\n\n\n\n");
	fflush(stdout);
	for (int i = 0; i < num_clusters; i++)
	{
		printf(" \n  Process [ %d ]  Cluster [ %d ] %lf  %lf  %lf \n ", myid, i + 1, clustersArray[i].x, clustersArray[i].y, clustersArray[i].z);
		fflush(stdout);
	}
	printf("\n\n\n\n");
	fflush(stdout);
}

//it is function that cold from kMeans to add all point centers to cluster center
void collectPointsCenters(Point* pointsPerProcess, Cluster* clustersArr, int num_of_points, int numClusters) {
	int i;
#pragma omp parallel for
	for (i = 0; i < numClusters; i++) {
		clustersArr[i].x = 0;
		clustersArr[i].y = 0;
		clustersArr[i].z = 0;
	}

	for (i = 0; i < num_of_points; i++) {
		clustersArr[pointsPerProcess[i].idCluster].x += pointsPerProcess[i].x;
		clustersArr[pointsPerProcess[i].idCluster].y += pointsPerProcess[i].y;
		clustersArr[pointsPerProcess[i].idCluster].z += pointsPerProcess[i].z;
	}
}