
#include "Headers.h"

MPI_Datatype createMPIDataTypePoint()
{
	Point point;
	MPI_Datatype MPI_Point;
	MPI_Datatype allType[7] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
		MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };
	int blocklen[7] = { 1, 1, 1 , 1 , 1, 1, 1 };
	MPI_Aint disp[7];

	disp[0] = (char *)&point.x - (char *)&point;
	disp[1] = (char *)&point.y - (char *)&point;
	disp[2] = (char *)&point.z - (char *)&point;
	disp[3] = (char *)&point.Vx - (char *)&point;
	disp[4] = (char *)&point.Vy - (char *)&point;
	disp[5] = (char *)&point.Vz - (char *)&point;
	disp[6] = (char *)&point.idCluster - (char *)&point;
	MPI_Type_create_struct(7, blocklen, disp, allType, &MPI_Point);
	MPI_Type_commit(&MPI_Point);
	return MPI_Point;
}


MPI_Datatype createMPIDataTypeCluster()
{
	Cluster cluster;
	MPI_Datatype typeForMPICluster;
	MPI_Datatype type[6] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE , MPI_INT ,MPI_INT };
	int blocklen[6] = { 1, 1, 1 , 1 , 1, 1 };
	MPI_Aint disp[6];

	disp[0] = (char *)&cluster.x - (char *)&cluster;
	disp[1] = (char *)&cluster.y - (char *)&cluster;
	disp[2] = (char *)&cluster.z - (char *)&cluster;
	disp[3] = (char *)&cluster.diameter - (char *)&cluster;
	disp[4] = (char *)&cluster.id - (char *)&cluster;
	disp[5] = (char *)&cluster.numPointsInCluster - (char *)&cluster;
	MPI_Type_create_struct(6, blocklen, disp, type, &typeForMPICluster);
	MPI_Type_commit(&typeForMPICluster);
	return typeForMPICluster;

}
