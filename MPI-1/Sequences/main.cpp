#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <iomanip> 
#include <ctime>

using namespace std;

#define ROOT 0

int main(int argc, char* argv[])
{
	double start, finish, serialtime, paralleltime, boost;
	int size(0), rank(0), len(0), serialresult(0), parallelresult(0);
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == ROOT)
	{
		cout << "Enter the size of the vector" << endl;
		cin >> len;
	}
	int *vector = new int[len + 1];
	if (rank == ROOT) 
	{
		srand(unsigned int(time(NULL)));
		if (len > 100)
		{
			for (int i = 0; i < len; i++)
				vector[i] = rand() % 100 - 50;
		}
		else
		{
			for (int i = 0; i < len; i++)
			{
				vector[i] = rand() % 100 - 50;
				cout << vector[i] << "  ";
			}
			cout << endl;
		}
		vector[len] = 0;
	}
	if (rank == ROOT)
	{
		start = MPI_Wtime();
		for (int i = 0; i < len; i++)
			if (vector[i] * vector[i + 1] < 0)
				serialresult++;
		finish = MPI_Wtime();
		serialtime = finish - start;
		cout << "Serial result " << serialresult << endl;
		cout << "Serial time " << serialtime << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == ROOT)
	{
		start = MPI_Wtime();
		int psize, lsize;
		parallelresult = 0;
		psize = len / size;
		for (int i = 1; i < size - 1; i++)
		{
			MPI_Send(&psize, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(vector + psize*i, psize + 1, MPI_INT, i, 2, MPI_COMM_WORLD);
		}
		if (size != 1)
		{
			lsize = len - psize*(size - 1);
			MPI_Send(&lsize, 1, MPI_INT, size - 1, 1, MPI_COMM_WORLD);
			MPI_Send(vector + psize*(size - 1), lsize + 1, MPI_INT, size - 1, 2, MPI_COMM_WORLD);
		}
		for (int i = 0; i < psize; i++)
			if (vector[i] * vector[i + 1] < 0)
				parallelresult++;
		int resValue;
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&resValue, 1, MPI_INT, i, 3, MPI_COMM_WORLD, &Status);
			parallelresult += resValue;
		}
	}
	if (rank != 0)
	{
		int psize, psum = 0;
		MPI_Recv(&psize, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &Status);
		int *pvector = new int[psize + 1];
		MPI_Recv(pvector, psize + 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &Status);
		for (int i = 0; i < psize; i++)
			if (pvector[i] * pvector[i + 1] < 0)
				psum++;
		MPI_Send(&psum, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
		delete[] pvector;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == ROOT)
	{
		finish = MPI_Wtime();
		paralleltime = finish - start;
		boost = serialtime / paralleltime;
		cout << "Parallel result " << parallelresult << endl;
		cout << "Parallel time " << paralleltime << endl;
		cout << "Boost is " << boost << endl;
	}
	delete[] vector;
	MPI_Finalize();
	return 0;
}
