#include <mpi.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>

using namespace std;

#define MATRIX_SIZE 8 // размер матарицы по умолчанию
#define THRESHOLD 64 // размер матриц для перехода к обычному умножению по умолчанию

void Input(int *Matrix1, int *Matrix2, int N) // формирование матриц
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

void Output(int *Result, int N) // печать матрицы
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

int *Add(int *Matrix1, int *Matrix2, int N) // сложение матриц
{
	int *result = new int[N*N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			result[i*N+j] = Matrix1[i*N+j] + Matrix2[i*N+j];
	return result;
}

int *Sub(int *Matrix1, int *Matrix2, int N) // вычитание матриц
{
	int *result = new int[N*N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			result[i*N+j] = Matrix1[i*N+j] - Matrix2[i*N+j];
	return result;
}

/* Обычный метод умножения матриц*/
int *StandardAlgorithm(int *Matrix1, int *Matrix2, int N) 
{
	int *result = new int[N*N];
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
int *StrassenAlgorithm(int *Matrix1, int *Matrix2, int N, int threshold)
{
	int *result;
	/* Когда размер матриц становится достаточно малым, 
	используем обычный метод умножения матриц, 
	так как алгоритм Штрассена теряет эффективность */
	if (N <= threshold)
		result = StandardAlgorithm(Matrix1, Matrix2, N);
	else
	{
		result = new int[N*N];
		N /= 2; // Матрица разбивается на 4 блока размерами N/2 * N/2

		/* Создаём вспомогательные матрицы, выделяем память */
		int *A[4], *B[4], *C[4], *P[7];
		for (int i = 0; i < 4; i++)
		{
			A[i] = new int [N*N];
			B[i] = new int [N*N];
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
	int rank, procnum, sqr, new_N, Size = MATRIX_SIZE, threshold = THRESHOLD, 
		**A, **B, **tmp, *Matrix1, *Matrix2, *result_s, *result_p;
	if (argc >= 2) Size = atoi(argv[1]);
	if (argc >= 3) threshold = atoi(argv[2]);
	double start_time, end_time, serial_time, parallel_time;
	result_s = new int [Size*Size];
	MPI_Status Status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		Matrix1 = new int [Size*Size], Matrix2 = new int [Size*Size], result_p = new int [Size*Size];
		Input(Matrix1, Matrix2, Size);
		start_time = MPI_Wtime();
		result_s = StrassenAlgorithm(Matrix1,Matrix2,Size,threshold);
		end_time = MPI_Wtime();
		cout << "Serial realisation: " << endl;
		Output(result_s, Size); // результат умножения матриц обычным методом
		serial_time = end_time - start_time;
		cout << "Serial time: " << serial_time << endl;

		//------------------------------------------------------------------------
		//      Параллельная реализация
		//------------------------------------------------------------------------

		start_time = MPI_Wtime();
		MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&threshold, 1, MPI_INT, 0, MPI_COMM_WORLD);
		sqr = (int)sqrt((double)procnum), new_N = Size/sqr;
		A = new int*[procnum], B = new int*[procnum];
		for(int i = 0; i < procnum; i++)
		{
			A[i] = new int [new_N*new_N];
			B[i] = new int [new_N*new_N];
		}
		for(int i = 0; i < Size; i++)
		{
			for(int j = 0; j < Size; j++)
			{
				A[sqr*(i/new_N)+j/new_N][(i%new_N)*new_N+(j%new_N)] = Matrix1[i*Size+j];
				B[sqr*(i/new_N)+j/new_N][(i%new_N)*new_N+(j%new_N)] = Matrix2[i*Size+j];
			}
		}
		/*Рассылка данных другим процессам*/
		MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(&threshold, 1, MPI_INT, 0, MPI_COMM_WORLD);
		for(int i = 1; i < procnum; i++)
		{
			int coef_A = sqr*(i / sqr), coef_B = i % sqr;
			for(int j = 0; j < sqr; j++)
			{
				MPI_Send(A[coef_A], new_N*new_N, MPI_INT, i , 0, MPI_COMM_WORLD); 
				MPI_Send(B[coef_B], new_N*new_N, MPI_INT, i , 0, MPI_COMM_WORLD); 
				coef_A++;
				coef_B += sqr;
			}
		}

		/*Swap указателей для однообразных вычислений*/
		for(int i = 0; i < sqr; i++)
		{
			int* С = B[i];
			B[i] = B[i*sqr];
			B[i*sqr] = С;
		}
	}
	/*Прием данных от процесса-root и формировка нужных данных*/
	if(rank != 0)
	{
		MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&threshold, 1, MPI_INT, 0, MPI_COMM_WORLD);
		sqr = (int)sqrt((double)procnum), new_N = Size/sqr;
		A = new int*[sqr], B = new int*[sqr];
		for(int i = 0; i < sqr; i++)
		{
			A[i] = new int [new_N*new_N];
			B[i] = new int [new_N*new_N];
		}
		for(int i = 0; i < sqr; i++)
		{
			MPI_Recv(A[i], new_N*new_N, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
			MPI_Recv(B[i], new_N*new_N, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
		}
	}
	/*вычисление каждым процессом своего куска матрицы*/
	tmp = new int*[sqr+1];
	for(int i = 0; i < sqr; i++)
		tmp[i+1] = StandardAlgorithm(A[i], B[i], new_N);
	if (procnum / 4 == 0)
		tmp[0] = tmp[1];
	if(procnum == 4)
		tmp[0] = Add(tmp[1], tmp[2], new_N);
	if(procnum == 16)
		tmp[0] = Add(Add(Add(tmp[1], tmp[2], new_N),tmp[3],new_N),tmp[4],new_N);
	/*Освобождение вспомогательной памяти*/
	if(procnum == 0)
	{
		for(int i = 0; i < procnum; i++)
		{
			delete[] A[i];
			delete[] B[i];
		}
	}
	else
	{
		for(int i = 0; i < sqr; i++)
		{
			delete[] A[i];
			delete[] B[i];
		}
	}
	for(int i = 1; i < sqr+1; i++)
		delete[] tmp[i];
	delete[] A;
	delete[] B;

	/*Отправка результата на 0 процесс*/
	if (rank != 0)	
	{
		MPI_Send(tmp[0], new_N*new_N, MPI_INT, 0, 0, MPI_COMM_WORLD);
		delete[] tmp[0];
	}
	if (rank == 0)
	{
		int coef = (int)sqrt((double)procnum);
		/*ззаписываем результат совей работы*/
		for(int i = 0; i < new_N; i++)
			for(int j = 0; j < new_N; j++)
				result_p[coef*i*new_N+j] = tmp[0][i*new_N+j];
		for(int k = 1; k < procnum; k++)
		{
			/*принимаем и записываем результаты работы других процессов*/
			MPI_Recv(tmp[0], new_N*new_N, MPI_INT, k, 0, MPI_COMM_WORLD, &Status);
			for(int i = 0; i < new_N; i++)
				for(int j = 0; j < new_N; j++)
					result_p[(k/coef)*new_N*Size+(k%coef)*new_N+coef*i*new_N+j] = tmp[0][i*new_N+j];
		}
	}
	if (rank == 0)
	{
		end_time = MPI_Wtime();
		cout << "Parallel realisation: " << endl;
		Output(result_p, Size); // результат умножения матриц обычным методом
		parallel_time = end_time - start_time;
		cout << "Parallel time: " << parallel_time << endl;
		bool flag = true;
		for (int i = 0; (i < Size * Size) && (flag == true); i++)
			if (result_s[i] != result_p[i])
				flag = false;
		if(flag)
			printf("Results are equal\n");
		else
			printf("Results are not equal\n");
		printf("Boost: %f\n", serial_time / parallel_time);
	}

	MPI_Finalize();
	return 0;
}