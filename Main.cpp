#include "Headers.h"
#include "StructsHeader.h"

int main(int argc, char *argv[])
{

	int  numprocs;        //number of proccesses
	int	myid;              //tred id
	int numofPoints;      //number of point that will be read from file
	int limit;              // maximum number of iterations for K-Means Algorithm.
	int num_clusters;    //number of point that will be read from file
	double t = 0;         // T in project guidelines that will be read from file
	double dt = 0;       // dt - time block  that will be read from file   
	double quality_measure = 0;  // quality measure to stop that will be read from file

	MPI_Status status;
	MPI_Comm comm;

	Point* allpoints = NULL; // whole dataset of input points
	Point* pointsBufferForEachProcess = NULL; // chunk of dataset - part of points (numofPoints/numprocesses) for each process.
	Cluster *clustersArray = NULL;// Array of Clusters.
	int numofPointsPerRank; //number of points that etch proccess will have for parallel work on k-means algoritem


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	// MPI Datatype variables
	MPI_Datatype MPI_POINT_TYPE = createMPIDataTypePoint();
	MPI_Datatype MPI_CLUSTER_TYPE = createMPIDataTypeCluster();

	//master will do  init file and will read points from file
	if (myid == MASTER) {
		initFile(INIT_FILENAME, NUMPOINTS);
		allpoints = readPointsFromFile(FILE_IN, &numofPoints, &num_clusters, &t, &dt, &limit, &quality_measure);
		printf("\n ----- Main calling after readPointsFromFile function -------  quality_measure %lf \n ", quality_measure);
		fflush(stdout);
	}

	// Master broadcasts first row in file to all other processes
	broadcastData(&numofPoints, &num_clusters, &t, &dt, &limit, &quality_measure);
	//function that teackes point that master have and restrebute them to all pccesses
	pointsBufferForEachProcess = scatterPoints(numofPoints, allpoints, &numofPointsPerRank, pointsBufferForEachProcess, myid, numprocs, MPI_POINT_TYPE);

	
	//-----------------------------------------------------parallel work begins---------------------------------------------------------------------------------------//

	// every proccess creates clasters array
	clustersArray = (Cluster*)calloc(num_clusters, sizeof(Cluster));
	//every proccess checking dynamic allocation of claster array
	checkDynamicAllocation(clustersArray);
	
	//every proccess  initialize all the Clusters 
	clustersArray = initClusters(num_clusters, allpoints, clustersArray, MPI_CLUSTER_TYPE, myid);

    //in checkForGoodClusters we will gather clusters and check thear qualety
	checkForGoodClusters(quality_measure, &numofPoints, &num_clusters, t, dt, limit, allpoints, &numofPointsPerRank, pointsBufferForEachProcess, &myid, &numprocs, MPI_POINT_TYPE, MPI_CLUSTER_TYPE, clustersArray);

	freeAll(allpoints, clustersArray, pointsBufferForEachProcess);
	MPI_Finalize();
	return 0;
}

