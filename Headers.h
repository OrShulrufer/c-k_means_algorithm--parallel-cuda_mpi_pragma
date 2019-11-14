#ifndef  HEADERS_H_
#define  HEADERS_H_
#pragma once


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <omp.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda.h"
#include "StructsHeader.h"


#define _CRT_SECURE_NO_WARNINGS
#define MASTER 0
#define MPI_FLAG 0
#define KMEANS_TERMINATE_FLAG 1
#define CLOCK_TERMINATE_FLAG 2
#define FILE_IN "D:\\FinalProject_310994710\\FinalProject_310994710\\inputValid.txt"
#define FILE_OUT "D:\\FinalProject_310994710\\FinalProject_310994710\\output.txt"
#define INIT_FILENAME "D:\\FinalProject_310994710\\FinalProject_310994710\\data.txt"
#define CONST 3
#define NUMPOINTS 50000
#define NUMCLUSTERS 30
#define LIMIT 5
#define QUALITY_MEASURE 1.0
#define T 5.0
#define DT 1.0
#define NUM_DIMENSIONS 3
#define TRUE 1
#define FALSE 0
#define MAXVALUE 10000



/* --------------------------- Function Prototypes ---------------------------*/
void printArray(Point* allPoints, int size, int processnumber);
void checkDynamicAllocation(const void* ptr);
void checkFile(const FILE* f);
void initFile(char* fileName, int numofPoints);
void writeToFile(int numofClusters, Cluster *clusterArray, double elapsedTime, double quality);
void broadcastData(int* numofpoints, int* numclusters, double* t, double *dt, int *limit, double* quality);
void calculateAverage(int numClusters, Cluster* arrClusters);
void mpiKmeans(Point* pointsPerProcess, int num_of_points, int numClusters, Cluster* clustersArr, int limit, int* myid);
void collectPoints(int myRank, Point *allPoints, int numOfAllPoints, int numPointsPerProcess, Point *PointsForProcess, MPI_Datatype MPI_POINT_TYPE, int numOfProcs);
void calculateQuality(int numberOfPoints, int numberOfClusters, Point* pointsArray, Cluster* clusterArray, double *quality);
void giveClustersDiameter(int numberOfPoints, int numOfClusters, Point* pointsArray, Cluster* clusterArray);
void toContinueKMeans(int* continuekflag);
void printGoodClusters(double quality_found, double time, Cluster* clustersArray, int num_clusters, int id);
void printClusters(Cluster* clustersArray, int num_clusters, int id);
void freeAll(Point* allpoints, Cluster* clustersArray, Point* pointsBufferForEachProcess);
void writeToFileNotFound();
void collectPointsCenters(Point* pointsPerProcess, Cluster* clustersArr, int num_of_points, int numClusters);

cudaError_t pointsLocation(int allPointsSize, double theTime, Point* pointsArray);

Cluster* checkForGoodClusters(double quality_measure, int* numofPoints, int* num_clusters, double t, double dt, int limit, Point* allpoints, int* numofPointsPerProcess, Point*  pointsPerProcess, int* myid, int* numprocs, MPI_Datatype  MPI_POINT_TYPE, MPI_Datatype MPI_CLUSTER_TYPE, Cluster* allClusters);
Cluster* initClusters(int numOfClusters, Point *pointsArray, Cluster* allClusters, MPI_Datatype MPI_CLUSTER_TYPE, int myRank);

Point* scatterPoints(int sizeOfAllPoints, Point *arrPoints, int* numofPointsPerRank, Point *PointsPerProcess, int myRank, int numOfProcs, MPI_Datatype MPI_POINT_TYPE);
Point* readPointsFromFile(const char* fileName, int* numofPoints, int* num_clusters, double* t, double* dt, int* limit, double* quality_measure);

MPI_Datatype createMPIDataTypePoint();
MPI_Datatype createMPIDataTypeCluster();

double euclid_dist_2(Point p, Cluster c);
double theDiameter(int numberOfPoints, int numOfClusters, Point* pointsArray, int clusterId);
double  generatRandomPositiveNegitiveValue(int value);
double distanceClusters(Cluster* c1, Cluster* c2);
double distancePoints(Point* point1, Point* point2);

int find_nearest_cluster(int numClusters, Point point, Cluster* clusters, int *myid);

#endif // ! HEADERS_H_
