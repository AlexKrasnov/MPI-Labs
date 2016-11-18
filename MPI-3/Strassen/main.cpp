#include <mpi.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>

using namespace std;

#define MATRIX_SIZE 8 // размер матарицы по умолчанию
#define THRESHOLD 1 // размер матриц для перехода к обычному умножению по умолчанию

void Input(double *Matrix1, double *Matrix2, int N) // формирование матриц
{
	cout << "\t\t***Multiplication dense matrices. Strassen Algorithm***\n\n";
	srand(unsigned int(time(NULL)));
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Matrix1[i*N+j] = rand() % 10 - 5;
			Matrix2[i*N+j] = rand() % 10 - 5;
		}
	}
	if (N <= 8)
	{
		cout << "Matrix1: " << endl;
		for (int i = 0; i < N; i++, cout << endl)
			for (int j = 0; j < N; j++)
				cout << Matrix1[i*N + j] << " ";
		cout << endl << "Matrix2: " << endl;
		for (int i = 0; i < N; i++, cout << endl)
			for (int j = 0; j < N; j++)
				cout << Matrix2[i*N + j] << " ";
		cout << endl;
	}
}

void Output(double *Result, int N) // печать матрицы
{
	if (N <= 8)
	{
		cout << "Result: " << endl;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
				cout << Result[i*N + j] << " ";
			cout << endl;
		}
		cout << endl;
	}
}

double *Add(double *Matrix1, double *Matrix2, int N) // сложение матриц
{
	double *result = new double[N*N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			result[i*N+j] = Matrix1[i*N+j] + Matrix2[i*N+j];
	return result;
}

double *Sub(double *Matrix1, double *Matrix2, int N) // вычитание матриц
{
	double *result = new double[N*N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			result[i*N+j] = Matrix1[i*N+j] - Matrix2[i*N+j];
	return result;
}

/* Обычный метод умножения матриц*/
double *StandardAlgorithm(double *Matrix1, double *Matrix2, int N) 
{
	double *result = new double[N*N];
	fill(result,result+N*N, 0);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			for (int k = 0; k < N; k++)
				result[i*N+j] += Matrix1[i*N+k] * Matrix2[k*N+j];
	return result;
}

//------------------------------------------------------------------------
//      Последовательная реализация
//------------------------------------------------------------------------

/* Алгоритм Штрассена умножения матриц*/
double *StrassenAlgorithm(double *Matrix1, double *Matrix2, int N, int threshold)
{
	double *result;
	/* Когда размер матриц становится достаточно малым, 
	используем обычный метод умножения матриц, 
	так как алгоритм Штрассена теряет эффективность */
	if (N <= threshold)
		result = StandardAlgorithm(Matrix1, Matrix2, N);
	else
	{
		result = new double[N*N];
		N /= 2; // Матрица разбивается на 4 блока размерами N/2 * N/2

		/* Создаём вспомогательные матрицы, выделяем память */
		double *A[4], *B[4], *C[4], *P[7];
		for (int i = 0; i < 4; i++)
		{
			A[i] = new double [N*N];
			B[i] = new double [N*N];
		}
		/* Разбиваем матрицы на 4 блока */
		for (int i =0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				//int index_new = i*N+j,index_old = 2*i*N+j, N_N = 2*N*N;
				A[0][i*N+j] = Matrix1[2*i*N+j];
				A[1][i*N+j] = Matrix1[2*i*N+j+N];
				A[2][i*N+j] = Matrix1[2*i*N+j+2*N*N];
				A[3][i*N+j] = Matrix1[2*i*N+j+2*N*N+N];

				B[0][i*N+j] = Matrix2[2*i*N+j];
				B[1][i*N+j] = Matrix2[2*i*N+j+N];
				B[2][i*N+j] = Matrix2[2*i*N+j+2*N*N];
				B[3][i*N+j] = Matrix2[2*i*N+j+2*N*N+N];
			}
		}

		/* Производим рекурсивный процесс N раз, до тех пор
		пока размер матриц не станет достаточно малым*/
		P[0] = StrassenAlgorithm(Add(A[0],A[3],N),Add(B[0],B[3],N),N,threshold);
		P[1] = StrassenAlgorithm(Add(A[2],A[3],N),B[0],N,threshold);
		P[2] = StrassenAlgorithm(A[0],Sub(B[1],B[3],N),N,threshold);
		P[3] = StrassenAlgorithm(A[3],Sub(B[2],B[0],N),N,threshold);
		P[4] = StrassenAlgorithm(Add(A[0],A[1],N),B[3],N,threshold);
		P[5] = StrassenAlgorithm(Sub(A[2],A[0],N),Add(B[0],B[1],N),N,threshold);
		P[6] = StrassenAlgorithm(Sub(A[1],A[3],N),Add(B[2],B[3],N),N,threshold);

		/* Находим результирующие блоки */
		C[0] = Sub(Add(Add(P[0],P[3],N),P[6],N),P[4],N); // P[0] + P[3] + P[6] - P[4]
		C[1] = Add(P[2],P[4],N); // P[2] + P[4];
		C[2] = Add(P[1],P[3],N); // P[1] + P[3];
		C[3] = Sub(Add(Add(P[0],P[2],N),P[5],N),P[1],N); // P[0] + P[2] + P[5] - P[1]

		/*Формируем результирующую матрицу*/
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				result[i*2*N+j] = C[0][i*N+j];
				result[i*2*N+j+N] = C[1][i*N+j];
				result[i*2*N+j+2*N*N] = C[2][i*N+j];
				result[i*2*N+j+2*N*N+N] = C[3][i*N+j];
			}

			for(int i = 0; i < 4; i++)
			{
				delete[] A[i];
				delete[] B[i];
				delete[] C[i];
			}
			for(int i = 0; i < 7; i++)
				delete[] P[i];
	}
	return result;
}

int main(int argc, char *argv[])
{
	int rank, procnum;
	int Size = MATRIX_SIZE, threshold = THRESHOLD;
	if (argc >= 2) Size = atoi(argv[1]);
	if (argc >= 3) threshold = atoi(argv[2]);
	double start_time, end_time, serial_time, parallel_time;
	double *Matrix1 = new double [Size*Size], *Matrix2 = new double [Size*Size];
	double *result_s = new double [Size*Size], *result_p = new double [Size*Size];
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		Input(Matrix1, Matrix2, Size);
		start_time = MPI_Wtime();
		result_s = StrassenAlgorithm(Matrix1,Matrix2,Size,threshold);
		end_time = MPI_Wtime();
		cout << "Serial realisation: " << endl;
		Output(result_s, Size); // результат умножения матриц обычным методом
		serial_time = end_time - start_time;
		cout << "Serial time: " << serial_time << endl;
	}



	MPI_Finalize();
	return 0;
}