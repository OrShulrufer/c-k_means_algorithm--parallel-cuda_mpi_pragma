#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "StructsHeader.h"

#define THREADS_RER_BLOCK 1000


__global__ void movePoints (Point*devPoints, int allPointsSize, unsigned int numofThreadsperBlock, unsigned int numofBlocks, double theTime, int flag) {
	int index = blockIdx.x * numofThreadsperBlock + threadIdx.x;
	/*Amount of Work per Thread */
	int threadwork = allPointsSize /(numofThreadsperBlock*(numofBlocks-flag));

	for (int i = index*threadwork; i < (index* threadwork) + threadwork; i++)
	{
	if (i < allPointsSize) {
		    devPoints[i].x  = devPoints[i].x + theTime*devPoints[i].Vx;
			devPoints[i].y = devPoints[i].y + theTime*devPoints[i].Vy;
			devPoints[i].z = devPoints[i].z + theTime*devPoints[i].Vz;
		}
	}
}

cudaError_t pointsLocation(int allPointsSize, double theTime, Point* pointsArray) {

	Point* devPoints = NULL;

	unsigned int numofThreadsperBlock = THREADS_RER_BLOCK;
	unsigned int numofBlocks = allPointsSize / numofThreadsperBlock;

	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		cudaFree(devPoints);
	}
	cudaStatus = cudaMalloc((void**)&devPoints, allPointsSize * sizeof(Point));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		cudaFree(devPoints);
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(devPoints, pointsArray, allPointsSize * sizeof(Point), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		cudaFree(devPoints);
	}

	int flag = 0;
	if (0 < allPointsSize % (numofThreadsperBlock)) {
	    numofBlocks += 1;
		flag = 1;
    }
		// Launch a kernel on the GPU
		movePoints << <numofBlocks, numofThreadsperBlock >> > (devPoints,  allPointsSize, numofThreadsperBlock, numofBlocks, theTime, flag);
	

		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "CalculateMaximumClusterDiameters launch failed: %s\n", cudaGetErrorString(cudaStatus));
		cudaFree(devPoints);
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		cudaFree(devPoints);
	}

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(pointsArray, devPoints, allPointsSize * sizeof(Point), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		cudaFree(devPoints);
	}

	cudaFree(devPoints);
	return cudaStatus;
}

